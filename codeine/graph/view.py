import random

from typing import Dict, FrozenSet, Generator, List, Optional, Sequence, Tuple, Union

from codeine.utils.display import format_banned_sequences, format_count, format_restrictions
from codeine.graph.graph import CodonGraph
from codeine.graph.nodes import CodonNode
from codeine.graph.tracking import AdvanceResult, BannedSequenceTracker
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
        self._banned_tracker = None

        self.valid_paths_by_choice = {}
        self.weight_mass_by_choice = {}
        self.next_node_by_choice = {}
        self.next_state_by_choice = {}
        self.next_key_by_choice = {}
        self.codon_pos_by_key = {}
        self.fixed_choice_by_key = {}

        self.n_valid_sequences = None
        self.samplers = {}
        self.sample_steps = {}
        self.initial_key = None

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

        key = self.initial_key
        next_key_by_choice = self.next_key_by_choice
        codon_pos_by_key = self.codon_pos_by_key
        fixed_choice_by_key = self.fixed_choice_by_key

        while key is not None:
            pos = codon_pos_by_key.get(key)

            if pos is None:
                choice = fixed_choice_by_key[key]
            else:
                start = (pos - 1) * 3
                choice = seq[start:start + 3]

            try:
                key = next_key_by_choice[key][choice]
            except KeyError:
                return False

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
        state = self._initial_banned_state()
        sequence = []

        final_node = self.graph.final_node
        valid_paths_by_choice = self.valid_paths_by_choice
        next_node_by_choice = self.next_node_by_choice
        next_state_by_choice = self.next_state_by_choice
        view_key = self._view_key

        while node is not final_node:
            key = view_key(node, state)
            choice_counts = valid_paths_by_choice[key]

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

            node = next_node_by_choice[key][chosen]
            state = next_state_by_choice[key][chosen]

        return ''.join(sequence)

    def sample(self) -> str:
        """
        Sample a DNA sequence from this graph view.
        """
        if self.n_valid_sequences == 0:
            raise ValueError('Cannot sample from an empty coding space.')

        key = self.initial_key
        sequence = []
        sample_steps = self.sample_steps

        while key is not None:
            emitted, key = sample_steps[key].sample()

            if emitted is not None:
                sequence.append(emitted)

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
        self._banned_tracker = BannedSequenceTracker(self.graph, self.banned_sequences)

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

    def _initial_banned_state(self) -> TrackerState:
        """
        Return the initial banned-sequence tracker state.
        """
        if not self.banned_sequences:
            return EMPTY_TRACKER_STATE

        return self._banned_tracker.initial_state

    def _advance_banned_state(self, state, node, choice) -> AdvanceResult:
        """
        Advance banned-sequence tracking after taking a graph step.
        """
        if self._banned_tracker.is_trivial:
            return AdvanceResult(banned=False, state=state)

        step = (node.pos, choice)
        return self._banned_tracker.advance(step, state)

    def _update_samplers(self) -> None:
        """
        Make samplers for each reachable graph state.

        ``samplers`` keeps the public/debug codon samplers keyed by
        ``(node, tracker_state)``. ``sample_steps`` is the compiled runtime
        sampler used by sample(); it returns ``(emitted_codon, next_key)`` so
        sampling does not need to recompute or look up graph transitions.
        """
        samplers = {}
        sample_steps = {}
        final_node = self.graph.final_node

        for key, choice_masses in self.weight_mass_by_choice.items():
            node, state = key

            if node is final_node:
                continue

            runtime_items = []
            runtime_weights = []

            for choice, mass in choice_masses.items():
                child = self.next_node_by_choice[key][choice]
                next_state = self.next_state_by_choice[key][choice]
                next_key = None if child is final_node else (child, next_state)
                emitted = choice if isinstance(node, CodonNode) else None

                runtime_items.append((emitted, next_key))
                runtime_weights.append(mass)

            if runtime_items:
                sample_steps[key] = Sampler(runtime_items, runtime_weights, rng=self._rng)

            if isinstance(node, CodonNode):
                codons = list(choice_masses)
                weights = [choice_masses[codon] for codon in codons]

                if codons:
                    samplers[key] = Sampler(codons, weights, rng=self._rng)

        self.samplers = samplers
        self.sample_steps = sample_steps

    def _choices_for_node(self, node: CodonNode) -> List[str]:
        """
        Return codon choices available to this node in this view.
        """
        if node.pos in self.pinned_codons:
            return self.pinned_codons[node.pos]

        return node.codons

    def _update_descendant_counts(self) -> None:
        """
        Calculate valid path counts and weight masses for each outgoing transition,
        tracking banned-sequence state as part of the cache key.

        Also compile legal transitions so contains(), sample(), and sequence_at()
        can follow precomputed state/node transitions instead of re-running the
        banned-sequence tracker at runtime.
        """

        valid_paths_by_choice = {}
        weight_mass_by_choice = {}
        next_node_by_choice = {}
        next_state_by_choice = {}
        next_key_by_choice = {}
        codon_pos_by_key = {}
        fixed_choice_by_key = {}
        total_cache = {}

        initial_state = self._banned_tracker.initial_state
        initial_key = self._view_key(self.graph.initial_node, initial_state)
        self.initial_key = initial_key

        stack = [(self.graph.initial_node, initial_state, False)]

        while stack:
            node, state, expanded = stack.pop()
            key = self._view_key(node, state)

            if key in total_cache:
                continue

            if node is self.graph.final_node:
                total_cache[key] = (1, 1.0)
                valid_paths_by_choice[key] = {}
                weight_mass_by_choice[key] = {}
                next_node_by_choice[key] = {}
                next_state_by_choice[key] = {}
                next_key_by_choice[key] = {}
                continue

            if isinstance(node, CodonNode):
                choices = self._choices_for_node(node)
                codon_pos_by_key[key] = node.pos
            else:
                choices = [node.sequence]
                fixed_choice_by_key[key] = node.sequence

            if not expanded:
                stack.append((node, state, True))

                for choice in choices:
                    child = node.transitions.get(choice)

                    if child is None:
                        continue

                    advance = self._advance_banned_state(state, node, choice)

                    if advance.banned:
                        continue

                    child_key = self._view_key(child, advance.state)

                    if child_key not in total_cache:
                        stack.append((child, advance.state, False))

                continue

            choice_counts = {}
            choice_masses = {}
            choice_next_nodes = {}
            choice_next_states = {}
            choice_next_keys = {}
            total_count = 0
            total_mass = 0.0

            for choice in choices:
                child = node.transitions.get(choice)

                if child is None:
                    continue

                advance = self._advance_banned_state(state, node, choice)

                if advance.banned:
                    continue

                child_key = self._view_key(child, advance.state)
                child_count, child_mass = total_cache[child_key]

                if isinstance(node, CodonNode):
                    choice_mass = self.graph.cw[choice] * child_mass
                else:
                    choice_mass = child_mass

                if child_count:
                    choice_counts[choice] = child_count
                    total_count += child_count

                    choice_next_nodes[choice] = child
                    choice_next_states[choice] = advance.state
                    choice_next_keys[choice] = None if child is self.graph.final_node else child_key

                if choice_mass:
                    choice_masses[choice] = choice_mass
                    total_mass += choice_mass

            valid_paths_by_choice[key] = choice_counts
            weight_mass_by_choice[key] = choice_masses
            next_node_by_choice[key] = choice_next_nodes
            next_state_by_choice[key] = choice_next_states
            next_key_by_choice[key] = choice_next_keys
            total_cache[key] = (total_count, total_mass)

        self.valid_paths_by_choice = valid_paths_by_choice
        self.weight_mass_by_choice = weight_mass_by_choice
        self.next_node_by_choice = next_node_by_choice
        self.next_state_by_choice = next_state_by_choice
        self.next_key_by_choice = next_key_by_choice
        self.codon_pos_by_key = codon_pos_by_key
        self.fixed_choice_by_key = fixed_choice_by_key
        self.n_valid_sequences = total_cache[initial_key][0]
