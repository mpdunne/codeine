import math
import random

from bisect import bisect_left
from dataclasses import dataclass, field
from itertools import islice
from typing import Dict, Iterator, List, NamedTuple, Optional, Protocol, Sequence, Tuple, Union

from codeine.constraints.base import ConstraintState
from codeine.graph.base import CodonGraph
from codeine.graph.nodes import CodonNode


class CompiledView(Protocol):
    """
    Sequence operations shared by compiled representations.

    Implementations own their traversal data and caches. Sampling uses the
    caller's random generator so a compiled result can be shared between views.
    """

    @property
    def n_valid_sequences(self) -> int:
        """
        The number of valid sequences, irrespective of sampling weights.
        """
        pass

    def contains(self, seq: str) -> bool:
        """
        Check membership, normalising the input to the graph's DNA/RNA alphabet.
        """
        pass

    def sample(self, rng: random.Random, n: Optional[int] = None) -> Union[str, List[str]]:
        """
        Sample one sequence, or n sequences, using the supplied generator.
        """
        pass

    def enumerate(self) -> Iterator[str]:
        """
        Yield valid sequences in graph order.
        """
        pass

    def enumerate_range(self, start: int = 0, stop: Optional[int] = None) -> Iterator[str]:
        """
        Enumerate valid sequences from start up to, but not including, stop.
        """
        pass

    def sequence_at(self, index: int) -> str:
        """
        Return a sequence by index, accepting negative indices.
        """
        pass

    def sequences_at(self, index_slice: slice) -> List[str]:
        """
        Return sequences selected by a Python slice.
        """
        pass


# The traversal state consists of the current graph position
# and the current state of each active constraint.
TraversalState = Tuple[int, Tuple[ConstraintState, ...]]


class ChoiceResult(NamedTuple):
    """
    Cached result of taking one graph choice from one compiled state. The
    "choice" is the graph edge label, i.e. a codon or a context sequence.

    Each ChoiceResult is specific to its location in the graph. The descendant
    counts and log mass are calculated iteratively by summing the values of
    downstream states.
    """
    choice: str
    descendant_count: int
    descendant_log_mass: float
    next_state_id: Optional[int]
    is_coding: bool


# A None cumulative distribution denotes equally weighted choices.
SamplingChoices = Tuple[Tuple[Tuple[str, bool, Optional[int]], ...], Optional[Tuple[float, ...]]]


@dataclass(frozen=True)
class FlatCompiledView:
    """
    A compiled snapshot supporting sequence queries and sampling.

    Topology and counts are fixed. Lazy sampling caches contain only choices and
    probabilities, so views can share this result while using independent random
    number generators.
    """
    graph: CodonGraph
    initial_state: TraversalState
    initial_state_id: int
    states: Tuple[TraversalState, ...]

    # Deep compiled transitions:
    # state ID -> ((choice, child state ID), ...)
    # These include every transition allowed by the graph and constraints,
    # before temporary view pins are applied.
    child_results_by_state_id: Tuple[Tuple[Tuple[str, int], ...], ...]

    # Compiled graph choices (lookup):
    # state ID -> choice -> ChoiceResult
    # Used for fast sequence validation and graph traversal.
    choices_by_state_id: Tuple[Dict[str, ChoiceResult], ...]

    # Compiled graph choices (iteration):
    # state ID -> ChoiceResults in graph order
    # Used for fast sampling and sequence enumeration.
    choice_results_by_state_id: Tuple[Tuple[ChoiceResult, ...], ...]

    n_valid_sequences: int

    _sampling_choices_by_state_id: List[Optional[SamplingChoices]] = field(
        default_factory=list,
        init=False,
        repr=False,
        compare=False,
    )

    def __post_init__(self) -> None:
        self._sampling_choices_by_state_id.extend([None] * len(self.states))

    def contains(self, seq: str) -> bool:
        """
        Check whether a DNA or RNA sequence belongs to this compiled space.
        """
        seq = self.graph.tt.normalise_sequence(seq)

        if len(seq) != len(self.graph.aa_seq) * 3:
            return False

        state_id = self.initial_state_id
        choices_by_state_id = self.choices_by_state_id
        states = self.states

        left_context_pos = self.graph.left_context_node.pos
        right_context_pos = self.graph.right_context_node.pos

        while state_id is not None:

            pos, _constraint_states = states[state_id]
            if pos == left_context_pos:
                node = self.graph.left_context_node
            elif pos == right_context_pos:
                node = self.graph.right_context_node
            else:
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

    def sample(self, rng: random.Random, n: Optional[int] = None) -> Union[str, List[str]]:
        """
        Sample one or more sequences using the caller's random number generator.
        """
        if self.n_valid_sequences == 0:
            raise ValueError('Cannot sample from an empty coding space.')

        if n is None:
            return self._sample(rng)

        if n < 0:
            raise ValueError('n must be non-negative.')

        return [self._sample(rng) for _ in range(n)]

    def enumerate(self) -> Iterator[str]:
        """
        Yield all valid coding sequences in graph order.
        """
        return self._iter_all_sequences()

    def enumerate_range(self, start: int = 0, stop: Optional[int] = None) -> Iterator[str]:
        """
        Enumerate valid sequences from start up to, but not including, stop.
        """
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
        Return the sequence at an index, accepting negative indices.
        """
        n_valid_sequences = self.n_valid_sequences

        if index < -n_valid_sequences or index >= n_valid_sequences:
            raise IndexError(f'Sequence index {index} out of range for {n_valid_sequences} valid sequences.')

        if index < 0:
            index += n_valid_sequences

        return self._sequence_at(index)

    def sequences_at(self, index_slice: slice) -> List[str]:
        """
        Return the sequences selected by a Python slice.
        """
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

    def _sampling_choices_for_state_id(self, state_id: int) -> Optional[SamplingChoices]:
        """
        Cache choices and cumulative probabilities without binding a random generator.
        """
        cached = self._sampling_choices_by_state_id[state_id]

        if cached is not None:
            return cached

        choices = tuple(
            result
            for result in self.choice_results_by_state_id[state_id]
            if result.descendant_log_mass != -math.inf
        )

        if not choices:
            return None

        log_masses = [result.descendant_log_mass for result in choices]
        cumulative = None

        if len(choices) > 1 and len(set(log_masses)) > 1:
            weights = self._convert_log_masses_to_sampler_weights(log_masses)
            total = sum(weights)
            running = 0.0
            probabilities = []

            for weight in weights:
                running += weight
                probabilities.append(running / total)

            cumulative = tuple(probabilities)

        items = tuple((result.choice, result.is_coding, result.next_state_id) for result in choices)
        cached = items, cumulative
        self._sampling_choices_by_state_id[state_id] = cached

        return cached

    def _sample(self, rng: random.Random) -> str:
        """
        Sample one sequence using cached choices and the caller's random generator.
        """
        state_id = self.initial_state_id
        sampling_choices_by_state_id = self._sampling_choices_by_state_id
        random_choice = rng.random
        sequence = []

        while state_id is not None:
            sampling_choices = sampling_choices_by_state_id[state_id]

            if sampling_choices is None:
                sampling_choices = self._sampling_choices_for_state_id(state_id)

                if sampling_choices is None:
                    raise ValueError('No valid choices available for sampling.')

            items, cumulative = sampling_choices

            if len(items) == 1:
                choice, is_coding, state_id = items[0]
            elif cumulative is None:
                choice, is_coding, state_id = items[int(random_choice() * len(items))]
            else:
                choice, is_coding, state_id = items[bisect_left(cumulative, random_choice())]

            if is_coding:
                sequence.append(choice)

        return ''.join(sequence)

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
