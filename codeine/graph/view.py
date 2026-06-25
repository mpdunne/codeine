import random

from itertools import islice
from typing import Dict, Generator, List, Optional, Sequence, Union

from codeine.constraints.base import PathConstraint
from codeine.constraints.banned import BannedSequenceTracker
from codeine.graph.base import CodonGraph, CodonRestriction
from codeine.translation.tables import TranslationTable
from codeine.translation.weights import CodonWeights
from codeine.utils.display import format_forbidden_motifs, format_count, format_restrictions
from codeine.utils.sampling import Seedable
from codeine.graph.compile import ViewCompiler


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
        self.path_constraint: Optional[PathConstraint] = None

        self.banned_tracker = BannedSequenceTracker(self.graph, self.banned_sequences)

        self._compiled = None
        self._requires_compile = True

        self.initial_state = None
        self.choices_by_state = {}
        self.choice_start_by_state = {}
        self.choice_results_by_state = {}
        self.codon_pos_by_state = {}
        self.fixed_choice_by_state = {}
        self.samplers = {}

    def __getitem__(self, index: Union[int, slice]) -> Union[str, List[str]]:
        """
        Return one valid sequence, or a list of valid sequences for a slice.

        Parameters
        ----------
        index
            Zero-based sequence index, or slice of sequence indices.

        Returns
        -------
        str or list of str
            The indexed valid coding sequence, or a list of valid coding sequences.
        """
        if isinstance(index, slice):
            return self.sequences_at(index)

        return self.sequence_at(index)

    def __iter__(self) -> Generator[str, None, None]:
        """
        Iterate over all valid sequences in this graph view.

        Yields
        ----------
        All valid sequences in the graph view, in order.
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
        if self._requires_compile:
            self.compile()

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
                *format_forbidden_motifs(
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

    def contains(self, seq: str) -> bool:
        """
        Check whether a coding sequence is contained in this view.

        Parameters
        ----------
        seq
            The sequence to check

        Returns
        -------
        True if and only if the sequence is contained in this coding space.
        """
        if self._requires_compile:
            self.compile()

        seq = seq.upper()

        if len(seq) != len(self.graph.aa_seq) * 3:
            return False

        state = self.initial_state
        choices_by_state = self.choices_by_state
        choice_start_by_state = self.choice_start_by_state
        fixed_choice_by_state = self.fixed_choice_by_state

        while state is not None:
            start = choice_start_by_state.get(state)

            if start is None:
                choice = fixed_choice_by_state[state]
            else:
                choice = seq[start:start + 3]

            result = choices_by_state[state].get(choice)

            if result is None:
                return False

            state = result.next_state

        return True

    def sample(self) -> str:
        """
        Sample a coding sequence from this graph view.

        Returns
        -------
        A random valid coding sequence that satisfies the provided constraints.
        """
        if self._requires_compile:
            self.compile()

        if self.n_valid_sequences == 0:
            raise ValueError('Cannot sample from an empty coding space.')

        state = self.initial_state
        sequence = []
        samplers = self.samplers

        while state is not None:
            choice, is_coding, state = samplers[state].sample()

            if is_coding:
                sequence.append(choice)

        return ''.join(sequence)

    def enumerate(self) -> Generator[str, None, None]:
        """
        Enumerate all valid sequences in this view.

        Yields
        ------
        str
            All valid coding sequences, one by one.
        """
        if self._requires_compile:
            self.compile()

        yield from self._iter_all_sequences()

    def enumerate_range(self, start: int = 0, stop: Optional[int] = None) -> Generator[str, None, None]:
        """
        Enumerate valid sequences from start up to, but not including, stop.

        Parameters
        ----------
        start
            The zero-based start from which to begin enumeration
        stop
            The zero-based enumeration stop.

        Yields
        -------
        str
            Sequences in the range, one by one.
        """
        if self._requires_compile:
            self.compile()

        n_sequences = self.n_valid_sequences

        if stop is None:
            stop = n_sequences

        if start < 0 or stop < start or stop > n_sequences:
            raise IndexError('Enumeration range is out of bounds.')

        if start == stop:
            return

        if start == 0 and stop == n_sequences:
            yield from self._iter_all_sequences()
            return

        if start == 0:
            yield from islice(self._iter_all_sequences(), stop)
            return

        yield from self._iter_sequence_range(start, stop)

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
            The indexed valid coding sequence.
        """
        if self._requires_compile:
            self.compile()

        if index < 0 or index >= self.n_valid_sequences:
            raise IndexError(
                f'Sequence index {index} out of range for '
                f'{self.n_valid_sequences} valid sequences.'
            )

        return next(self._iter_sequence_range(index, index + 1))

    def sequences_at(self, index_slice: slice) -> List[str]:
        """
        Return valid sequences from a slice.
        """
        if self._requires_compile:
            self.compile()

        n_sequences = self.n_valid_sequences
        start, stop, step = index_slice.indices(n_sequences)

        if start == stop:
            return []

        if step != 1:
            return [self.sequence_at(index) for index in range(start, stop, step)]

        if start == 0 and stop == n_sequences:
            return [*self._iter_all_sequences()]

        if start == 0:
            return [*islice(self._iter_all_sequences(), stop)]

        return [*self._iter_sequence_range(start, stop)]

    def copy(self) -> 'CodonGraphView':
        """
        Copy this view and all its constraints and attributes.

        Returns
        -------
        A copy of the view.
        """
        view = self.graph.view()
        view._rng.setstate(self._rng.getstate())

        view.pinned_codons = self.pinned_codons.copy()
        view.banned_sequences = self.banned_sequences.copy()
        view.path_constraint = self.path_constraint
        view.banned_tracker = self.banned_tracker

        view._compiled = self._compiled
        view._requires_compile = self._requires_compile

        view.initial_state = self.initial_state
        view.choices_by_state = self.choices_by_state
        view.choice_start_by_state = self.choice_start_by_state
        view.choice_results_by_state = self.choice_results_by_state
        view.codon_pos_by_state = self.codon_pos_by_state
        view.fixed_choice_by_state = self.fixed_choice_by_state
        view.samplers = self.samplers

        return view

    def compile(self) -> None:
        """
        Calculate all graph properties that are derived from its structure plus constraints
        such as pins and banned sequences.

        Remember to do this after editing any constraints!
        """
        compiler = ViewCompiler(self)
        compiled = compiler.compile()

        self._compiled = compiled
        self._requires_compile = False

        self.initial_state = compiled.initial_state
        self.choices_by_state = compiled.choices_by_state
        self.choice_results_by_state = compiled.choice_results_by_state
        self.choice_start_by_state = compiled.choice_start_by_state
        self.codon_pos_by_state = compiled.codon_pos_by_state
        self.fixed_choice_by_state = compiled.fixed_choice_by_state
        self.samplers = compiled.samplers

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
        self._requires_compile = True

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

        self._requires_compile = True

    def set_pinned_codons(self, pinned_codons: Dict[int, CodonRestriction]) -> None:
        """
        Pin (temporarily fix) a specified group codons, leaving all others unpinned.

        Parameters
        ----------
        pinned_codons:
            A dict specifying which codons to pin, by pos: codon
        """
        pinned_codons = self.graph.validate_codon_restrictions(pinned_codons)
        self.pinned_codons = dict(pinned_codons)
        self._requires_compile = True

    def clear_pins(self) -> None:
        """
        Remove all codon pins from this graph view
        """
        self.pinned_codons.clear()
        self._requires_compile = True

    def set_banned_sequences(self, banned_sequences: Sequence[str]) -> None:
        """
        Set banned nucleotide sequences for this view.

        Banned-sequence tracking depends only on the graph and banned sequences,
        not on temporary pins, so it is rebuilt only when the banned list changes.
        """
        self.banned_sequences = self._validate_banned_sequences(banned_sequences)
        self.banned_tracker = BannedSequenceTracker(self.graph, self.banned_sequences)
        self._requires_compile = True

    def clear_banned_sequences(self) -> None:
        """
        Remove all banned sequence restrictions from this view.
        """
        self.set_banned_sequences([])

    def set_path_constraint(self, path_constraint: Optional[PathConstraint]) -> None:
        """
        Set an additional generic path constraint for this view.

        Pass None to remove any path constraint.
        """
        self.path_constraint = path_constraint
        self._requires_compile = True

    def clear_path_constraint(self) -> None:
        """
        Remove the additional generic path constraint from this view.
        """
        self.set_path_constraint(None)

    @property
    def aa_seq(self) -> str:
        """
        The amino acid sequence.

        Returns
        -------
        The aa seq.
        """
        return self.graph.aa_seq

    @property
    def translation_table(self) -> TranslationTable:
        """
        The translation table used by the codon graph.
        """
        return self.graph.tt

    @property
    def codon_weights(self) -> CodonWeights:
        """
        The codon weights used by the codon graph.
        """
        return self.graph.cw

    @property
    def codon_restrictions(self) -> Dict[int, CodonRestriction]:
        """
        Any hard-fixed codon restrictions on the codon graph.
        """
        return self.graph.codon_restrictions

    @property
    def context_l(self) -> str:
        """
        The left context sequence.
        """
        return self.graph.context_l

    @property
    def context_r(self) -> str:
        """
        The right context sequence.
        """
        return self.graph.context_r

    @property
    def n_valid_sequences(self) -> int:
        """
        Number of valid coding sequences in this view given all constraints.
        """
        if self._requires_compile:
            self.compile()

        return self._compiled.n_valid_sequences

    def _validate_banned_sequences(self, banned_sequences: Optional[Sequence[str]]) -> List[str]:
        """
        Check the inputted banned sequences make sense, and return normalised versions of them.

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
            sequence = self.translation_table.normalise_sequence(sequence)

            if len(sequence) == 0:
                raise ValueError('Banned sequences cannot be empty.')

            normalised.append(sequence)

        return sorted(set(normalised))

    def _iter_all_sequences(self) -> Generator[str, None, None]:
        """
        Iterate over all valid sequences. Faster than _iter_sequence_range when
        we're starting at 0.

        Yields
        ------
        str
            All valid coding sequences, one by one
        """
        # Stack is:
        # (
        #       state,
        #       coding sequence constructed so far,
        # )
        choice_results_by_state = self.choice_results_by_state
        codon_pos_by_state = self.codon_pos_by_state

        stack = [(self.initial_state, '')]

        while stack:
            state, prefix = stack.pop()

            if state is None:
                yield prefix
                continue

            results = choice_results_by_state[state]

            if not results:
                continue

            if state not in codon_pos_by_state:
                stack.append((results[0].next_state, prefix))
                continue

            for result in reversed(results):
                stack.append((result.next_state, prefix + result.choice))

    def _iter_sequence_range(
        self,
        start: int,
        stop: int,
    ) -> Generator[str, None, None]:
        """
        Iterate over valid sequences in a given index range.

        Parameters
        ----------
        start
            0-based index of the first sequence.
        stop
            0-based index one past the final sequence.

        Yields
        ------
        str
            Valid coding sequences in the requested range.
        """
        # Stack is:
        # (
        #       state,
        #       sequence constructed so far,
        #       0-based index of the first sequence reachable from that state.
        # )
        choice_results_by_state = self.choice_results_by_state
        codon_pos_by_state = self.codon_pos_by_state

        stack = [(self.initial_state, '', 0)]

        while stack:
            state, prefix, offset = stack.pop()

            if state is None:
                if start <= offset < stop:
                    yield prefix
                continue

            results = choice_results_by_state[state]

            if not results:
                continue

            if state not in codon_pos_by_state:
                stack.append((results[0].next_state, prefix, offset))
                continue

            child_start = offset
            push = []

            for result in results:
                child_stop = child_start + result.descendant_count

                if child_stop > start and child_start < stop:
                    push.append((result, child_start))

                child_start = child_stop

            for result, child_start in reversed(push):
                stack.append((
                    result.next_state,
                    prefix + result.choice,
                    child_start,
                ))

