import math
import random

from bisect import bisect_left
from dataclasses import dataclass, field, replace
from itertools import islice
from typing import Dict, Iterator, List, NamedTuple, Optional, Protocol, Sequence, Tuple, TYPE_CHECKING, Union

from codeine.constraints.base import Constraint, ConstraintState, DEAD_STATE, SAFE_STATE
from codeine.graph.base import CodonGraph
from codeine.graph.nodes import CodonNode

if TYPE_CHECKING:
    from codeine.graph.view import CodonGraphView


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
TraversalStateKey = TraversalState


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


class ViewCompiler:
    """
    Compile a CodonGraphView into cached topology, choice, count, and sampling-mass data.
    """

    def __init__(self, view: 'CodonGraphView') -> None:
        self.view = view
        self.graph = view.graph

        self.constraints: Tuple[Constraint, ...] = ()
        self.constraint_advancers = ()

        self.state_ids: Dict[TraversalStateKey, int] = {}
        self.states: List[TraversalState] = []

        # Dynamic-programming totals:
        # state ID -> (descendant count, descendant log mass)
        # Avoids repeatedly recomputing subtree sizes and probability masses.
        self.totals_by_state_id: List[Optional[Tuple[int, float]]] = []

        # Deep compiled transitions:
        # state ID -> [(choice, child state ID), ...]
        # These account for graph restrictions and constraints, but not
        # temporary view pins.
        self.child_results_by_state_id: List[Optional[List[Tuple[str, int]]]] = []

        # Compiled graph choices (lookup):
        # state ID -> choice -> ChoiceResult
        # Used for fast sequence validation and graph traversal.
        self.choices_by_state_id: List[Optional[Dict[str, ChoiceResult]]] = []

        # The cached log-ified codon weights, to avoid repeated log calculations.
        self.log_codon_weights = {
            codon: math.log(weight) if weight > 0 else -math.inf
            for codon, weight in self.view.codon_weights.weights.items()
        }

        self.initial_pos = self.graph.initial_node.pos
        self.final_pos = self.graph.final_node.pos
        self.positions = tuple(node.pos for node in self.graph.nodes)
        self.seq_len = len(self.graph.aa_seq)

        # Graph transitions available at each position. Permanent graph-level codon
        # restrictions are already reflected in node.transitions. Temporary view
        # pins are applied later, during the shallow calculation pass.
        self.transitions_by_pos = {
            node.pos: tuple((choice, child.pos) for choice, child in node.transitions.items())
            for node in self.graph.nodes
            if node is not self.graph.final_node
        }

    def _set_constraints(self, constraints: Sequence[Constraint]) -> None:
        """
        Link and configure the constraints used by this compilation pass.

        Parameters
        ----------
        constraints
            Constraints whose states should be represented in the compiled
            topology.
        """
        constraints = tuple(constraints)

        for constraint in constraints:
            constraint.link(self.graph)

        self.constraints = tuple(
            constraint
            for constraint in constraints
            if not constraint.is_trivial
        )
        self.constraint_advancers = tuple(
            constraint.advance
            for constraint in self.constraints
        )

    def compile(self) -> FlatCompiledView:
        """
        Compile descendant counts, graph choices, and sampling masses.

        Returns
        -------
        FlatCompiledView
            A compiled view.
        """
        self._set_constraints(self.view.constraints)

        initial_state = self._initial_state()
        initial_pos, initial_constraint_states = initial_state
        initial_state_id, _ = self._get_or_register_state_id(
            initial_pos,
            initial_constraint_states,
        )

        self._compile_topology(initial_state_id)
        self._compile_choices()

        return self._compiled_view(initial_state_id)

    def compile_shallow(self, compiled: FlatCompiledView) -> FlatCompiledView:
        """
        Recompile choices, counts, and probability masses using an existing
        deep topology.

        Parameters
        ----------
        compiled
            The existing compiled view whose states and transitions should be reused.

        Returns
        -------
        FlatCompiledView
            The compiled view with updated shallow data.
        """
        self.states = list(compiled.states)
        self.child_results_by_state_id = list(compiled.child_results_by_state_id)

        self.totals_by_state_id = [None] * len(self.states)
        self.choices_by_state_id = [None] * len(self.states)

        self._compile_choices()

        choices_by_state_id = tuple(
            choices or {}
            for choices in self.choices_by_state_id
        )

        choice_results_by_state_id = tuple(
            tuple(choices.values()) if choices else ()
            for choices in self.choices_by_state_id
        )

        initial_total = self.totals_by_state_id[compiled.initial_state_id]
        assert initial_total is not None

        return replace(
            compiled,
            choices_by_state_id=choices_by_state_id,
            choice_results_by_state_id=choice_results_by_state_id,
            n_valid_sequences=initial_total[0],
        )

    def extend(
            self,
            compiled: FlatCompiledView,
            constraints: Sequence[Constraint],
    ) -> FlatCompiledView:
        """
        Extend an existing compiled topology with additional constraints.

        The existing compiled states and transitions already encode all previous
        constraints. Extension therefore traverses that topology directly and
        advances only the newly supplied constraints.

        Parameters
        ----------
        compiled
            The existing compiled view.
        constraints
            Additional constraints to compile.

        Returns
        -------
        FlatCompiledView
            An updated compiled view.
        """
        self._set_constraints(constraints)

        if not self.constraints:
            return self.compile_shallow(compiled)

        initial_new_states = tuple(constraint.initial_state for constraint in self.constraints)
        initial_pos, old_initial_states = compiled.initial_state
        initial_state_id, _ = self._get_or_register_state_id(initial_pos, old_initial_states + initial_new_states)

        self._compile_extended_topology(compiled, initial_state_id)
        self._compile_choices()

        return self._compiled_view(initial_state_id)

    def _get_or_register_state_id(
            self,
            pos: int,
            constraint_states: Tuple[ConstraintState, ...],
    ) -> Tuple[int, bool]:
        """
        Return the state ID and whether the state was newly registered.

        State lookup uses the graph position rather than the node object, keeping
        the hot dictionary key compact.
        """
        key = (pos, constraint_states)
        state_id = self.state_ids.get(key)

        if state_id is not None:
            return state_id, False

        state_id = len(self.states)
        self.state_ids[key] = state_id
        self.states.append(key)
        self.totals_by_state_id.append(None)
        self.choices_by_state_id.append(None)
        self.child_results_by_state_id.append(None)

        return state_id, True

    def _initial_state(self) -> TraversalState:
        """
        Return the starting traversal state.

        The initial state starts at the graph's initial position, with fresh
        constraint states.

        Returns
        -------
        TraversalState
            Starting state for graph compilation.
        """
        constraint_states = tuple(constraint.initial_state for constraint in self.constraints)

        return self.initial_pos, constraint_states

    def _compile_topology(self, initial_state_id: int) -> None:
        """
        Discover every reachable traversal state and transition.

        Temporary view pins are not applied during this pass.

        Parameters
        ----------
        initial_state_id
            ID of the state from which graph compilation should begin.
        """
        stack = [initial_state_id]

        while stack:
            state_id = stack.pop()

            if self.child_results_by_state_id[state_id] is not None:
                continue

            pos, _constraint_states = self.states[state_id]

            if pos == self.final_pos:
                self.child_results_by_state_id[state_id] = []
                continue

            stack.extend(self._uncompiled_children(state_id))

    def _compile_extended_topology(self, compiled: FlatCompiledView, initial_state_id: int) -> None:
        """
        Compile additional constraints over an existing compiled topology.

        Parameters
        ----------
        compiled
            The existing compiled view on which to build.
        initial_state_id
            The initial state ID.
        """
        n_new_constraints = len(self.constraints)

        old_state_ids = [compiled.initial_state_id]

        stack = [initial_state_id]

        while stack:
            state_id = stack.pop()

            if self.child_results_by_state_id[state_id] is not None:
                continue

            old_state_id = old_state_ids[state_id]
            pos, constraint_states = self.states[state_id]

            if pos == self.final_pos:
                self.child_results_by_state_id[state_id] = []
                continue

            new_constraint_states = constraint_states[-n_new_constraints:]
            child_results = []

            for choice, old_child_id in compiled.child_results_by_state_id[old_state_id]:
                next_new_states = self._advance_constraints(new_constraint_states, pos,  choice)

                if next_new_states is None:
                    continue

                child_pos, old_child_states = compiled.states[old_child_id]

                child_id, is_new = self._get_or_register_state_id(child_pos, old_child_states + next_new_states)

                child_results.append((choice, child_id))

                if is_new:
                    old_state_ids.append(old_child_id)
                    stack.append(child_id)

            self.child_results_by_state_id[state_id] = child_results

    def _compile_choices(self) -> None:
        """
        Compile active choices, descendant counts, and probability masses.

        Temporary view pins and codon weights are applied during this pass.

        States are processed from right to left through the graph, ensuring
        that every child has been compiled before its parent.
        """
        state_ids_by_pos = {pos: [] for pos in self.positions}

        for state_id, (pos, _constraint_states) in enumerate(self.states):
            state_ids_by_pos[pos].append(state_id)

        for pos in reversed(self.positions):
            for state_id in state_ids_by_pos[pos]:
                if pos == self.final_pos:
                    self._compile_final_state(state_id)
                else:
                    self._compile_state(state_id)

    def _compile_final_state(self, state_id: int) -> None:
        """
        Compile a terminal traversal state.

        By the time a terminal state is reached, all graph choices have already
        been processed, including the right context. Choices rejected by the
        constraints would not have reached this state.

        Parameters
        ----------
        state_id
            ID of the terminal traversal state being compiled.
        """
        self.totals_by_state_id[state_id] = (1, 0.0)
        self.choices_by_state_id[state_id] = {}

    def _compile_state(self, state_id: int) -> None:
        """
        Compile one non-final traversal state.

        For each outgoing graph choice allowed by the current pins, combine the
        previously compiled child state with the contribution from the current
        graph position to produce a ChoiceResult. The total descendant count and
        log mass are then cached for the current state.

        Parameters
        ----------
        state_id
            ID of the traversal state being compiled.
        """
        pos, _constraint_states = self.states[state_id]

        choice_results = {}
        descendant_count = 0

        max_log_mass = -math.inf
        relative_mass_sum = 0.0

        is_coding = 1 <= pos <= self.seq_len
        child_results = self.child_results_by_state_id[state_id] or ()

        pinned_codons = (self.view.pinned_codons.get(pos) if is_coding else None)

        for choice, child_id in child_results:
            if pinned_codons is not None and choice not in pinned_codons:
                continue

            child_pos, _child_constraint_states = self.states[child_id]
            child_total = self.totals_by_state_id[child_id]

            if child_total is None:
                continue

            child_count, subtree_log_mass = child_total

            if child_count == 0:
                continue

            if is_coding:
                codon_log_weight = self.log_codon_weights[choice]
                choice_log_mass = codon_log_weight + subtree_log_mass
            else:
                choice_log_mass = subtree_log_mass

            result = ChoiceResult(
                choice=choice,
                descendant_count=child_count,
                descendant_log_mass=choice_log_mass,
                next_state_id=None if child_pos == self.final_pos else child_id,
                is_coding=is_coding,
            )

            choice_results[choice] = result
            descendant_count += child_count

            if choice_log_mass == -math.inf:
                continue

            # Incremental log-sum-exp.
            if choice_log_mass <= max_log_mass:
                relative_mass_sum += math.exp(choice_log_mass - max_log_mass)
            else:
                if max_log_mass == -math.inf:
                    relative_mass_sum = 1.0
                else:
                    relative_mass_sum = relative_mass_sum * math.exp(max_log_mass - choice_log_mass) + 1.0

                max_log_mass = choice_log_mass

        if max_log_mass == -math.inf:
            descendant_log_mass = -math.inf
        else:
            descendant_log_mass = max_log_mass + math.log(relative_mass_sum)

        self.choices_by_state_id[state_id] = choice_results
        self.totals_by_state_id[state_id] = (descendant_count, descendant_log_mass)

    def _compiled_view(self, initial_state_id: int) -> FlatCompiledView:
        """
        Build an immutable FlatCompiledView from the compiler's current state.

        Parameters
        ----------
        initial_state_id
            The initial state ID.
        """
        initial_total = self.totals_by_state_id[initial_state_id]
        assert initial_total is not None

        choices_by_state_id = tuple(choices or {} for choices in self.choices_by_state_id)

        return FlatCompiledView(
            graph=self.graph,
            initial_state=self.states[initial_state_id],
            initial_state_id=initial_state_id,
            states=tuple(self.states),
            child_results_by_state_id=tuple(tuple(results or ()) for results in self.child_results_by_state_id),
            choices_by_state_id=choices_by_state_id,
            choice_results_by_state_id=tuple(tuple(choices.values()) for choices in choices_by_state_id),
            n_valid_sequences=initial_total[0],
        )

    def _uncompiled_children(self, state_id: int) -> List[int]:
        """
        Return child state IDs reached by taking each outgoing graph choice.

        Choices rejected by constraints are skipped. Temporary view pins are not
        applied during this pass. Only child states that have not yet been
        compiled are returned.

        Parameters
        ----------
        state_id
            ID of the traversal state whose children should be discovered.

        Returns
        -------
        list of int
            Child state IDs still needing compilation.
        """
        child_results = []
        uncompiled_children = []

        pos, constraint_states = self.states[state_id]

        for choice, child_pos in self.transitions_by_pos[pos]:
            next_constraint_states = self._advance_constraints(constraint_states, pos, choice)

            if next_constraint_states is None:
                continue

            child_id, is_new = self._get_or_register_state_id(child_pos, next_constraint_states)

            child_results.append((choice, child_id))

            if is_new:
                uncompiled_children.append(child_id)

        self.child_results_by_state_id[state_id] = child_results

        return uncompiled_children

    def _advance_constraints(
            self,
            constraint_states: Tuple[ConstraintState, ...],
            pos: int,
            choice: str,
    ) -> Optional[Tuple[ConstraintState, ...]]:
        """
        Advance all active constraints after taking one graph choice.

        Parameters
        ----------
        constraint_states
            Current state of each constraint.
        pos
            Current graph position.
        choice
            Graph choice taken from the current position.

        Returns
        -------
        tuple or None
            The updated constraint states, or None if any constraint rejects the
            graph choice.
        """
        next_states = []
        append = next_states.append

        for advance, state in zip(self.constraint_advancers, constraint_states):
            if state == SAFE_STATE:
                append(SAFE_STATE)
                continue

            next_state = advance(state, pos, choice)

            if next_state == DEAD_STATE:
                return None

            append(next_state)

        return tuple(next_states)
