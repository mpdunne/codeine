import random

from typing import Dict, FrozenSet, Generator, List, Optional, Sequence, Tuple, Union

from codeine.utils.display import format_banned_sequences, format_count, format_restrictions
from codeine.graph.graph import CodonGraph, CodonNode
from codeine.utils.sampling import Sampler, Seedable


CodonRestriction = Union[str, Sequence[str]]
Watch = Tuple[int, int]
TrackerState = FrozenSet[Watch]
ViewKey = Tuple[object, TrackerState]
EMPTY_TRACKER_STATE: TrackerState = frozenset()


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

    def _view_key(self, node, state: TrackerState = EMPTY_TRACKER_STATE) -> ViewKey:
        """
        Return the cache key for a node plus banned-sequence tracker state.
        """
        return node, state

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
            if pos < 1 or pos > len(self.graph.codon_nodes):
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
        seq = seq.upper()

        if len(seq) != len(self.graph.aa_seq) * 3:
            return False

        current_node = self.graph.initial_node.transitions[self.graph.initial_node.sequence]

        while current_node is not self.graph.right_context_node:
            pos = current_node.pos
            codon = seq[(pos - 1) * 3: pos * 3]

            if pos in self.pinned_codons:
                if codon not in self.pinned_codons[pos]:
                    return False

            elif codon not in current_node.codons:
                return False

            current_node = current_node.transitions[codon]

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
        state = EMPTY_TRACKER_STATE
        sequence = []

        while node is not self.graph.final_node:
            choice_counts = self.valid_paths_by_choice[self._view_key(node, state)]

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
                # Context nodes only have one valid outgoing choice.
                chosen = next(iter(choice_counts.keys()))

            node = node.transitions[chosen]

        return ''.join(sequence)

    def sample(self) -> str:
        """
        Sample a DNA sequence from this graph view.
        """
        if self.n_valid_sequences == 0:
            raise ValueError('Cannot sample from an empty coding space.')

        node = self.graph.initial_node
        state = EMPTY_TRACKER_STATE
        sequence = []

        while node is not self.graph.final_node:
            if not node.transitions:
                raise RuntimeError(f'Reached non-final node {node.id} with no outgoing transitions.')

            if isinstance(node, CodonNode):
                choice = self.samplers[self._view_key(node, state)].sample()
                sequence.append(choice)

            else:
                choice = node.sequence

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
        # Placeholder!
        self._check_banned_sequence_support()

        # Calculate descendant counts!
        self._update_descendant_counts()

        # Update the samplers!
        self._update_samplers()

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

    def _check_banned_sequence_support(self) -> None:
        """
        For now, we don't support banned sequences. We mock this in tests.
        """
        if self.banned_sequences:
            raise NotImplementedError('Banned sequences are not yet supported.')

    def _update_samplers(self) -> None:
        """
        Make samplers for each codon node and banned-sequence tracker state.

        The probabilities are calculated by combining descendant weight masses
        and the base codon probabilities.
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

        right_choice = self.graph.right_context_node.sequence
        right_key = self._view_key(self.graph.right_context_node)

        valid_paths_by_choice[right_key] = {right_choice: 1}
        weight_mass_by_choice[right_key] = {right_choice: 1.0}

        next_counts = {right_key: 1}
        next_masses = {right_key: 1.0}

        for node in reversed(self.graph.codon_nodes):
            current_counts = {}
            current_masses = {}

            choice_counts = {}
            choice_masses = {}
            total_count = 0
            total_mass = 0.0

            for codon in self._choices_for_node(node):
                child = node.transitions[codon]
                child_key = self._view_key(child)

                count = next_counts.get(child_key, 0)
                mass = self.graph.cw[codon] * next_masses.get(child_key, 0.0)

                if count:
                    choice_counts[codon] = count
                    total_count += count

                if mass:
                    choice_masses[codon] = mass
                    total_mass += mass

            node_key = self._view_key(node)

            valid_paths_by_choice[node_key] = choice_counts
            weight_mass_by_choice[node_key] = choice_masses

            current_counts[node_key] = total_count
            current_masses[node_key] = total_mass

            next_counts = current_counts
            next_masses = current_masses

        left_choice = self.graph.left_context_node.sequence
        left_child = self.graph.left_context_node.transitions.get(left_choice)
        left_key = self._view_key(self.graph.left_context_node)

        if left_child is None:
            valid_paths_by_choice[left_key] = {}
            weight_mass_by_choice[left_key] = {}
            self.valid_paths_by_choice = valid_paths_by_choice
            self.weight_mass_by_choice = weight_mass_by_choice
            self.n_valid_sequences = 0
            return

        left_child_key = self._view_key(left_child)

        total_count = next_counts.get(left_child_key, 0)
        total_mass = next_masses.get(left_child_key, 0.0)

        valid_paths_by_choice[left_key] = {left_choice: total_count}
        weight_mass_by_choice[left_key] = {left_choice: total_mass}

        self.valid_paths_by_choice = valid_paths_by_choice
        self.weight_mass_by_choice = weight_mass_by_choice
        self.n_valid_sequences = total_count
