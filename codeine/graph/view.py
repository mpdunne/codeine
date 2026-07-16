import math
import random

from itertools import islice
from typing import Dict, Iterator, List, Optional, Sequence, Tuple, Union

from codeine.constraints.base import Constraint
from codeine.graph.base import CodonGraph, CodonRestriction
from codeine.graph.compile import ViewCompiler
from codeine.graph.nodes import CodonNode
from codeine.translation.tables import TranslationTable
from codeine.translation.weights import CodonWeights
from codeine.utils.display import format_count, format_restrictions
from codeine.utils.sampling import Seedable, Sampler, SingletonSampler, UniformSampler, WeightedSampler


# Sentinels for compile status, indicating what, if
# anything, needs compiling/recompiling.
COMPILED = 0
COMPILE_SHALLOW = 1
COMPILE_DEEP = 2


class CodonGraphView:
    """
    View of a codon graph with optional constraints and temporary codon pins.

    It is on this object that most operations (counting, sampling, enumeration....) take place.
    """

    def __init__(self,
                 graph: CodonGraph,
                 *,
                 constraints: Optional[Sequence[Constraint]] = None,
                 weights: Optional[CodonWeights] = None,
                 seed: Seedable = None,
                 ) -> None:
        """
        Constructor for the CodonGraphView

        Parameters
        ----------
        graph
            The underlying codon graph.
        constraints
            Any constraint trackers that we wish to use when traversing coding space.
        weights
            The codon weights to use when sampling.
        seed
            Seed used to initialise a random number generator.
        """

        self.graph = graph
        self.pinned_codons: Dict[int, CodonRestriction] = {}
        self.constraints: Tuple[Constraint, ...] = tuple(constraints or ())

        if weights is None:
            weights = CodonWeights.uniform(table=self.graph.tt)
        else:
            weights = weights.for_table(self.graph.tt)

        self._codon_weights = weights
        self._rng = random.Random(seed)

        self._compiled = None
        self._compile_status = COMPILE_DEEP

        self.initial_state_id = None
        self.choices_by_state_id = ()
        self.choice_results_by_state_id = ()
        self.samplers_by_state_id = []

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

    def __iter__(self) -> Iterator[str]:
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
        if self._compile_status:
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
        if self._compile_status:
            self.compile()

        seq = self.translation_table.normalise_sequence(seq)

        if len(seq) != len(self.graph.aa_seq) * 3:
            return False

        state_id = self.initial_state_id
        choices_by_state_id = self.choices_by_state_id
        states = self._compiled.states

        while state_id is not None:

            pos, _constraint_states = states[state_id]
            node = self.graph.codon_node_by_pos(pos)

            if isinstance(node, CodonNode):
                start = (pos - 1) * 3
                choice = seq[start:start + 3]
            else:
                choice = node.sequence

            result = choices_by_state_id[state_id].get(choice)

            if result is None:
                return False

            state_id = result.next_state_id

        return True

    def sample(self, n: Optional[int] = None) -> Union[str, List[str]]:
        """
        Sample one or more coding sequences from this graph view.

        Parameters
        ----------
        n
            Number of sequences to sample. If omitted, return a single sequence.

        Returns
        -------
        str or list of str
            One sampled coding sequence, or a list of sampled coding sequences.
        """
        if self._compile_status:
            self.compile()

        if self.n_valid_sequences == 0:
            raise ValueError('Cannot sample from an empty coding space.')

        if n is None:
            return self._sample()

        if n < 0:
            raise ValueError('n must be non-negative.')

        return [self._sample() for _ in range(n)]

    def enumerate(self) -> Iterator[str]:
        """
        Enumerate all valid sequences in this view.

        Yields
        ------
        str
            All valid coding sequences, one by one.
        """
        if self._compile_status:
            self.compile()

        yield from self._iter_all_sequences()

    def enumerate_range(self, start: int = 0, stop: Optional[int] = None) -> Iterator[str]:
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
        if self._compile_status:
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
        if self._compile_status:
            self.compile()

        if index < 0 or index >= self.n_valid_sequences:
            raise IndexError(
                f'Sequence index {index} out of range for '
                f'{self.n_valid_sequences} valid sequences.'
            )

        return self._sequence_at(index)

    def sequences_at(self, index_slice: slice) -> List[str]:
        """
        Return valid sequences from a slice.
        """
        if self._compile_status:
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
        view._codon_weights = self._codon_weights

        view.pinned_codons = self.pinned_codons.copy()
        view.constraints = self.constraints

        view._compiled = self._compiled
        view._compile_status = self._compile_status

        view.initial_state_id = self.initial_state_id
        view.choices_by_state_id = self.choices_by_state_id
        view.choice_results_by_state_id = self.choice_results_by_state_id
        view.samplers_by_state_id = [None] * len(self.samplers_by_state_id)

        return view

    def compile(self) -> None:
        """
        Compile or recompile the view as required.
        """
        if self._compile_status == COMPILED:
            return

        compiler = ViewCompiler(self)

        if self._compile_status == COMPILE_SHALLOW and self._compiled is not None:
            compiled = compiler.compile_shallow(self._compiled)
        else:
            compiled = compiler.compile()

        self._compiled = compiled
        self._compile_status = COMPILED

        self.initial_state_id = compiled.initial_state_id
        self.choices_by_state_id = compiled.choices_by_state_id
        self.choice_results_by_state_id = compiled.choice_results_by_state_id
        self.samplers_by_state_id = [None] * len(compiled.states)

    def _update_compile_status(self, status: int) -> None:
        """
        Mark this view as requiring at least the specified compile phase.
        """
        self._compile_status = max(self._compile_status, status)

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
        self._update_compile_status(COMPILE_SHALLOW)

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

        self._update_compile_status(COMPILE_SHALLOW)

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
        self._update_compile_status(COMPILE_SHALLOW)

    def clear_pins(self) -> None:
        """
        Remove all codon pins from this graph view
        """
        self.pinned_codons.clear()
        self._update_compile_status(COMPILE_SHALLOW)

    def set_constraints(self, constraints: Sequence[Constraint]) -> None:
        """
        Set the constraints for this view.

        Parameters
        ----------
        constraints
            Constraints to apply during graph traversal.
        """
        self.constraints = tuple(constraints)
        self._update_compile_status(COMPILE_DEEP)

    def clear_constraints(self) -> None:
        """
        Remove all constraints from this view.
        """
        self.set_constraints(())

    def set_weights(self, weights: Optional[CodonWeights] = None) -> None:
        """
        Set the codon weights used for sampling.

        Parameters
        ----------
        weights
            Codon weights to use. If omitted, use uniform weights.
        """
        if weights is None:
            weights = CodonWeights.uniform(table=self.graph.tt)
        else:
            weights = weights.for_table(self.graph.tt)

        self._codon_weights = weights
        self._update_compile_status(COMPILE_SHALLOW)

    def clear_weights(self) -> None:
        """
        Reset sampling to uniform codon weights.
        """
        self.set_weights()

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
        return self._codon_weights

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
        if self._compile_status:
            self.compile()

        return self._compiled.n_valid_sequences

    def _sampler_for_state_id(self, state_id: int) -> Optional[Sampler]:
        """
        Return the appropriate sampler for one compiled traversal state.

        Samplers are created lazily because many compiled states may never be
        visited during sampling.

        Parameters
        ----------
        state_id
            The state ID.

        Returns
        -------
        Sampler or None
            The cached sampler for this state, or None if the state has no
            choices with positive sampling mass.
        """
        sampler = self.samplers_by_state_id[state_id]

        if sampler is not None:
            return sampler

        choice_results = self.choice_results_by_state_id[state_id]

        if not choice_results:
            return None

        runtime_items = []
        runtime_log_masses = []

        for result in choice_results:
            if result.descendant_log_mass == -math.inf:
                continue

            runtime_items.append((result.choice, result.is_coding, result.next_state_id))
            runtime_log_masses.append(result.descendant_log_mass)

        if not runtime_items:
            return None

        if len(runtime_items) == 1:
            sampler = SingletonSampler(item=runtime_items[0])
        elif len(set(runtime_log_masses)) == 1:
            sampler = UniformSampler(items=runtime_items, rng=self._rng)
        else:
            runtime_weights = self._convert_log_masses_to_sampler_weights(runtime_log_masses)
            sampler = WeightedSampler(items=runtime_items, weights=runtime_weights, rng=self._rng,)

        self.samplers_by_state_id[state_id] = sampler

        return sampler

    @staticmethod
    def _convert_log_masses_to_sampler_weights(log_masses: Sequence[float]) -> List[float]:
        """
        Convert subtree log masses into relative weights for sampling.

        The returned weights are proportional to the true subtree probabilities but
        are rescaled to avoid numerical underflow. Only the relative values matter
        for weighted sampling.

        Parameters
        ----------
        log_masses
            Choice masses represented in log space.

        Returns
        -------
        list of float
            Relative non-log weights suitable for weighted sampling.
        """
        max_log_mass = max(log_masses)
        return [math.exp(log_mass - max_log_mass) for log_mass in log_masses]

    def _sample(self) -> str:
        """
        Sample one coding sequence from an already-compiled graph view.

        Returns
        -------
        A sampled sequence.
        """
        state_id = self.initial_state_id
        sequence = []

        while state_id is not None:
            sampler = self._sampler_for_state_id(state_id)

            if sampler is None:
                raise ValueError('No valid choices available for sampling.')

            choice, is_coding, state_id = sampler.sample()

            if is_coding:
                sequence.append(choice)

        return ''.join(sequence)

    def _sequence_at(self, index: int) -> str:
        """
        Return one valid sequence by directly descending through descendant counts.

        Parameters
        ----------
        index
            The index of the sequence in the graph.

        Returns
        -------
        The sequence at the desired index.
        """
        state_id = self.initial_state_id
        choice_results_by_state_id = self.choice_results_by_state_id
        sequence_parts = []

        while state_id is not None:
            results = choice_results_by_state_id[state_id]

            if not results:
                raise RuntimeError('Unexpected dead end during sequence index traversal.')

            if not results[0].is_coding:
                state_id = results[0].next_state_id
                continue

            for result in results:
                descendant_count = result.descendant_count

                if index < descendant_count:
                    sequence_parts.append(result.choice)
                    state_id = result.next_state_id
                    break

                index -= descendant_count
            else:
                raise RuntimeError('Invalid sequence index traversal state.')

        return ''.join(sequence_parts)

    def _iter_all_sequences(self) -> Iterator[str]:
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
        choice_results_by_state_id = self.choice_results_by_state_id
        sequence_parts = [''] * len(self.graph.aa_seq)

        stack = [(self.initial_state_id, 0, None)]

        while stack:
            state_id, codon_index, choice = stack.pop()

            if choice is not None:
                sequence_parts[codon_index - 1] = choice

            if state_id is None:
                yield ''.join(sequence_parts)
                continue

            results = choice_results_by_state_id[state_id]

            if not results:
                continue

            if not results[0].is_coding:
                stack.append((results[0].next_state_id, codon_index, None))
                continue

            next_codon_index = codon_index + 1

            for result in reversed(results):
                stack.append((result.next_state_id, next_codon_index, result.choice))

    def _iter_sequence_range(
        self,
        start: int,
        stop: int,
    ) -> Iterator[str]:
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
        choice_results_by_state_id = self.choice_results_by_state_id
        sequence_parts = [''] * len(self.graph.aa_seq)

        stack = [(self.initial_state_id, 0, None, 0)]

        while stack:
            state_id, codon_index, choice, offset = stack.pop()

            if choice is not None:
                sequence_parts[codon_index - 1] = choice

            if state_id is None:
                if start <= offset < stop:
                    yield ''.join(sequence_parts)
                continue

            results = choice_results_by_state_id[state_id]

            if not results:
                continue

            if not results[0].is_coding:
                stack.append((results[0].next_state_id, codon_index, None, offset))
                continue

            next_codon_index = codon_index + 1
            child_start = offset
            push = []

            for result in results:
                child_stop = child_start + result.descendant_count

                if child_stop > start and child_start < stop:
                    push.append((result, child_start))

                child_start = child_stop

            for result, child_start in reversed(push):
                stack.append((
                    result.next_state_id,
                    next_codon_index,
                    result.choice,
                    child_start,
                ))
