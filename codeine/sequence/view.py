import random

from typing import Dict, Generator, List, Optional, Sequence, Union

from codeine.sequence.display import format_banned_sequences, format_count, format_restrictions
from codeine.sequence.graph import CodonGraph, CodonNode
from codeine.utils.sampling import Sampler, Seedable


CodonRestriction = Union[str, Sequence[str]]


class CodonGraphView:
    """
    View of a codon graph. The view allows optional temporary constraints such as pinned codons
    and banned nucleotide sequences to be added without affecting the underlying codon graph.

    It is on this object that most operations (counting, sampling, enumeration....) take place.
    """

    def __init__(self,
                 graph: CodonGraph,
                 banned_sequences: Optional[Sequence[str]] = None,
                 seed: Seedable = None,
                 rng: Optional[random.Random] = None,
                 ) -> None:
        """
        Constructor for the CodonGraphView

        Parameters
        ----------
        graph
            The underlying codon graph.
        banned_sequences
            Nucleotide sequences that are forbidden in this view.
        seed
            Seed used to initialise a random number generator, if not providing an RNG.
        rng
            Random number generator used by the view for sampling.
        """
        if seed is not None and rng is not None:
            raise ValueError('Provide either seed or rng, not both.')

        self._rng = rng if rng is not None else random.Random(seed)

        self.graph = graph
        self.pinned_codons: Dict[int, List[str]] = {}
        self.banned_sequences: List[str] = self._validate_banned_sequences(banned_sequences)

        self.valid_paths_by_choice = {}
        self.weight_mass_by_choice = {}

        self.n_valid_sequences = None
        self.samplers = {}

        self.compile()

    @property
    def aa_seq(self):
        """
        The amino acid sequence on the underlying graph.

        Returns
        -------
        The aa seq.
        """
        return self.graph.aa_seq

    def __getitem__(self, index: int) -> str:
        """
        Return the valid sequence at a given index.

        Parameters
        ----------
        index
            Zero-based sequence index.

        Returns
        -------
        str
            The indexed valid DNA sequence.
        """
        return self.sequence_at(index)

    def __iter__(self) -> Generator[str, None, None]:
        """
        Iterate over all valid sequences in this graph view.

        Yields
        ----------
        All valid sequences in the graph, in order.
        """
        yield from self.enumerate()

    def __contains__(self, seq: str) -> bool:
        """
        Does the given seq exist in this space?

        Returns
        ----------
        True if and only if this is a valid sequence in this space.
        """
        return self.contains(seq)

    def __repr__(self) -> str:
        molecule = 'RNA' if self.graph.tt.rna else 'DNA'

        lines = [
            f'{type(self).__name__}',
            '',
            f'Translation table: {self.graph.tt.table_id} ({self.graph.tt.name})',
            f'Molecule type: {molecule}',
            '',
            f'Amino acid sequence ({len(self.aa_seq)} aa)',
            f'{self.aa_seq}',
            ''
        ]

        if self.graph.codon_restrictions:
            lines += [
                'Codon restrictions:',
                *format_restrictions(
                    self.graph.codon_restrictions,
                    label='restricted positions',
                ),
                '',
                ]

        if self.banned_sequences:
            lines += [
                'Banned sequences:',
                *format_banned_sequences(
                    self.banned_sequences,
                ),
                '',
                ]

        if self.pinned_codons:
            lines += [
                'Temporary pins:',
                *format_restrictions(
                    self.pinned_codons,
                    label='pinned positions',
                ),
                '',
                ]

        lines.append(f'Num. valid coding sequences: {format_count(self.n_valid_sequences)}')

        return '\n'.join(lines)

    def pin_codons(self, pinned_codons: Dict[int, CodonRestriction]) -> None:
        """
        Pin (temporarily fix) a codon in this codon graph view

        Parameters
        ----------
        pinned_codons
            A dict specifying which codons to pin, by pos: codon.
        """
        pinned_codons = self.graph.validate_codon_restrictions(pinned_codons)
        self.pinned_codons.update(pinned_codons)
        self.compile()

    def unpin_codons(self, positions: Sequence[int]) -> None:
        """
        Unpin codon nodes by pos.

        Parameters
        ----------
        positions
            A list of positions to unpin.
        """
        for pos in positions:
            if pos not in self.graph.codon_nodes_by_pos:
                raise ValueError(f'Pinned codon position {pos} is out of range.')

            self.pinned_codons.pop(pos, None)

        self.compile()

    def set_pinned_codons(self, pinned_codons: Dict[int, CodonRestriction]) -> None:
        """
        Pin a specified group codons, leaving all others unpinned.

        Parameters
        ----------
        pinned_codons:
            A dict specifying which codons to pin, by pos: codon
        """
        pinned_codons = self.graph.validate_codon_restrictions(pinned_codons)
        self.pinned_codons = dict(pinned_codons)
        self.compile()

    def clear_pins(self) -> None:
        """
        Remove all codon pins from this graph view
        """
        self.pinned_codons.clear()
        self.compile()

    def set_banned_sequences(self, banned_sequences: Sequence[str]) -> None:
        """
        Set banned nucleotide sequences for this view.
        """
        self.banned_sequences = self._validate_banned_sequences(banned_sequences)
        self.compile()

    def contains(self, seq: str) -> bool:
        """
        Check whether a DNA sequence is contained in this view.

        Parameters
        ----------
        seq
            The sequence to check

        Returns
        -------
        True if and only if the sequence is contained in this coding space.
        """
        seq = self.graph.normalise_sequence(seq)

        if len(seq) != len(self.graph.aa_seq) * 3:
            return False

        node = self.graph.initial_node
        state = ''

        while node is not self.graph.final_node:
            if isinstance(node, CodonNode):
                pos = node.pos
                choice = seq[(pos - 1) * 3: pos * 3]

                if pos in self.pinned_codons:
                    if choice not in self.pinned_codons[pos]:
                        return False

                elif choice not in node.codons:
                    return False

            else:
                choice = next(iter(node.transitions.keys()))

            next_state = self._advance_banned_state(state, choice)

            if next_state is None:
                return False

            state = next_state

            if choice not in node.transitions:
                return False

            node = node.transitions[choice]

        return True

    def sequence_at(self, index: int) -> str:
        """
        Return the valid sequence at a given index.

        Parameters
        ----------
        index
            Zero-based sequence index.

        Returns
        -------
        str
            The indexed valid DNA sequence.
        """
        if index < 0 or index >= self.n_valid_sequences:
            raise IndexError(f'Sequence index {index} out of range for {self.n_valid_sequences} valid sequences.')

        node = self.graph.initial_node
        state = ''
        sequence = []

        while node is not self.graph.final_node:
            key = (node, state)
            choice_counts = self.valid_paths_by_choice[key]

            if isinstance(node, CodonNode):
                remaining = index

                for choice, count in choice_counts.items():
                    if remaining < count:
                        chosen = choice
                        break

                    remaining -= count
                else:
                    raise RuntimeError('Failed to resolve sequence index.')

                sequence.append(chosen)
                index = remaining

            else:
                chosen = next(iter(choice_counts.keys()))

            next_state = self._advance_banned_state(state, chosen)

            if next_state is None:
                raise RuntimeError('Resolved invalid banned-sequence transition.')

            state = next_state
            node = node.transitions[chosen]

        return ''.join(sequence)

    def sample(self) -> str:
        """
        Sample a DNA sequence from this graph view.
        """
        if self.n_valid_sequences == 0:
            raise ValueError('Cannot sample from an empty coding space.')

        node = self.graph.initial_node
        state = ''
        sequence = []

        while node is not self.graph.final_node:
            if not node.transitions:
                raise RuntimeError(f'Reached non-final node {node.id} with no outgoing transitions.')

            key = (node, state)

            if isinstance(node, CodonNode):
                choice = self.samplers[key].sample()
                sequence.append(choice)
            else:
                choice = next(iter(self.valid_paths_by_choice[key].keys()))

            next_state = self._advance_banned_state(state, choice)

            if next_state is None:
                raise RuntimeError('Sampled invalid banned-sequence transition.')

            state = next_state
            node = node.transitions[choice]

        return ''.join(sequence)

    def enumerate(self) -> Generator[str, None, None]:
        """
        Enumerate all valid sequences in this view.

        Yields
        ------
        str
            A valid DNA sequence.
        """
        for index in range(self.n_valid_sequences):
            yield self[index]

    def copy(self) -> 'CodonGraphView':
        """
        Copy this view and all its constraints.

        Returns
        -------
        A copy of the view.
        """
        view = self.graph.view()
        view.pinned_codons = self.pinned_codons.copy()
        view.banned_sequences = self.banned_sequences.copy()
        view.compile()
        return view

    def compile(self) -> None:
        """
        Calculate all graph properties that are derived from its structure plus constraints
        such as pins and banned sequences.

        Remember to do this after editing constraints!
        """
        # Calculate descendant counts!
        self._update_descendant_counts()

        # Update the samplers!
        self._update_samplers()

    @property
    def _max_banned_sequence_length(self) -> int:
        if not self.banned_sequences:
            return 0

        return max(len(sequence) for sequence in self.banned_sequences)

    def _advance_banned_state(self, state: str, emitted: str) -> Optional[str]:
        if not self.banned_sequences:
            return ''

        combined = state + emitted

        for banned_sequence in self.banned_sequences:
            if banned_sequence in combined:
                return None

        keep = self._max_banned_sequence_length - 1

        if keep <= 0:
            return ''

        return combined[-keep:]

    def _validate_banned_sequences(self, banned_sequences: Sequence[str]) -> List[str]:
        """
        Check the inputted banned sequences make sense.

        Parameters
        ----------
        banned_sequences
            The list of banned sequences.

        Returns
        -------
        A normalised, de-duplicated list of banned sequences.
        """
        banned_sequences = banned_sequences or []

        normalised = []
        for sequence in banned_sequences:
            sequence = self.graph.normalise_sequence(sequence)

            if len(sequence) == 0:
                raise ValueError('Banned sequences cannot be empty.')

            normalised.append(sequence)

        return sorted(set(normalised))

    def _update_samplers(self) -> None:
        """
        Make samplers for each codon node given the view's constraints. The probabilities
        are calculated by combining descendant counts and the base probabilities.
        """
        samplers = {}

        for key, choice_masses in self.weight_mass_by_choice.items():
            node, state = key

            if not isinstance(node, CodonNode):
                continue

            codons = list(choice_masses)
            weights = [choice_masses[codon] for codon in codons]

            if codons:
                samplers[key] = Sampler(codons, weights, rng=self._rng)

        self.samplers = samplers

    def _choices_for_node(self, node: CodonNode) -> List[str]:
        """
        Return codon choices available to this node in this view.
        """
        if node.pos in self.pinned_codons:
            return self.pinned_codons[node.pos]

        return node.codons

    def _update_descendant_counts(self) -> None:
        """
        Calculate valid path counts and weight masses for each outgoing transition.
        """
        valid_paths_by_choice = {}
        weight_mass_by_choice = {}

        # Final node: one completed path from any valid incoming state.
        next_counts = {('', self.graph.final_node): 1}
        next_masses = {('', self.graph.final_node): 1.0}

        # Right context layer.
        right_node = self.graph.right_context_node
        right_choice = right_node.sequence

        current_counts = {}
        current_masses = {}

        states = {state for state, node in next_counts if node is self.graph.final_node}

        for state in states:
            next_state = self._advance_banned_state(state, right_choice)

            key = (right_node, state)

            if next_state is None:
                valid_paths_by_choice[key] = {}
                weight_mass_by_choice[key] = {}
                continue

            child = right_node.transitions[right_choice]
            child_key = (next_state, child)

            count = next_counts.get(child_key, 0)
            mass = next_masses.get(child_key, 0.0)

            valid_paths_by_choice[key] = {right_choice: count} if count else {}
            weight_mass_by_choice[key] = {right_choice: mass} if mass else {}

            current_counts[(state, right_node)] = count
            current_masses[(state, right_node)] = mass

        next_counts = current_counts
        next_masses = current_masses

        # Codon layers, right to left.
        for pos in range(len(self.graph.aa_seq), 0, -1):
            current_counts = {}
            current_masses = {}

            states = {state for state, node in next_counts}

            for node in self.graph.codon_nodes_by_pos[pos]:
                for state in states:
                    key = (node, state)

                    choice_counts = {}
                    choice_masses = {}
                    total_count = 0
                    total_mass = 0.0

                    for codon in self._choices_for_node(node):
                        next_state = self._advance_banned_state(state, codon)

                        if next_state is None:
                            continue

                        child = node.transitions[codon]
                        child_key = (next_state, child)

                        child_count = next_counts.get(child_key, 0)
                        child_mass = next_masses.get(child_key, 0.0)

                        mass = self.graph.cw[codon] * child_mass

                        if child_count:
                            choice_counts[codon] = child_count
                            total_count += child_count

                        if mass:
                            choice_masses[codon] = mass
                            total_mass += mass

                    valid_paths_by_choice[key] = choice_counts
                    weight_mass_by_choice[key] = choice_masses

                    current_counts[(state, node)] = total_count
                    current_masses[(state, node)] = total_mass

            next_counts = current_counts
            next_masses = current_masses

        # Left context layer.
        left_node = self.graph.left_context_node
        left_choice = left_node.sequence
        initial_state = ''

        next_state = self._advance_banned_state(initial_state, left_choice)
        key = (left_node, initial_state)

        if next_state is None:
            valid_paths_by_choice[key] = {}
            weight_mass_by_choice[key] = {}
            total_count = 0
        else:
            child = left_node.transitions[left_choice]
            child_key = (next_state, child)

            total_count = next_counts.get(child_key, 0)
            total_mass = next_masses.get(child_key, 0.0)

            valid_paths_by_choice[key] = {left_choice: total_count} if total_count else {}
            weight_mass_by_choice[key] = {left_choice: total_mass} if total_mass else {}

        self.valid_paths_by_choice = valid_paths_by_choice
        self.weight_mass_by_choice = weight_mass_by_choice
        self.n_valid_sequences = total_count
