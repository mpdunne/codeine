import math

from dataclasses import dataclass
from typing import Dict, Tuple, List, Optional, TYPE_CHECKING

from codeine.constraints.banned import TrackerState, AdvanceResult
from codeine.constraints.base import ConstraintState
from codeine.graph.nodes import CodonNode, Node
from codeine.utils.sampling import Sampler

if TYPE_CHECKING:
    from codeine.graph.view import CodonGraphView


NodeState = Tuple[Node, TrackerState, ConstraintState]


@dataclass(frozen=True)
class ChoiceResult:
    """
    Cached result of taking one graph choice from one compiled state. The
    "choice" is the graph edge label, i.e. a codon or a context sequence.

    Each ChoiceResult is specific to its location in the graph. The descendant
    counts and weight mass are calculated iteratively by summing the values of
    downstream nodes.
    """
    choice: str
    descendant_count: int
    descendant_weight_mass: float
    next_state: Optional[NodeState]
    is_coding: bool


@dataclass
class CompiledView:
    """
    Cached data for a compiled CodonGraphView, to speed up sampling and enumeration.
    """
    initial_state: NodeState
    n_valid_sequences: int
    choices_by_state: Dict[NodeState, Dict[str, ChoiceResult]]
    choice_start_by_state: Dict[NodeState, int]
    choice_results_by_state: Dict[NodeState, Tuple[ChoiceResult, ...]]
    codon_pos_by_state: Dict[NodeState, int]
    fixed_choice_by_state: Dict[NodeState, str]
    samplers: dict


class ViewCompiler:
    """
    Compile a CodonGraphView into cached choice, count, and sampling data.
    """

    def __init__(self, view: 'CodonGraphView') -> None:
        self.view = view
        self.graph = view.graph
        self.tracker = view._banned_tracker
        self.advance_cache = view._advance_cache
        self.path_constraint = view.path_constraint
        self.has_path_constraint = self.path_constraint is not None

        self.totals_by_state: Dict[NodeState, Tuple[int, float]] = {}
        self.choices_by_state: Dict[NodeState, Dict[str, ChoiceResult]] = {}
        self.choice_start_by_state: Dict[NodeState, int] = {}
        self.choice_results_by_state: Dict[NodeState, Tuple[ChoiceResult, ...]] = {}
        self.codon_pos_by_state: Dict[NodeState, int] = {}
        self.fixed_choice_by_state: Dict[NodeState, str] = {}

        self.choices_by_node = {
            node: tuple(self._get_choices_for_node(node))
            for node in self.graph.nodes
            if node is not self.graph.final_node
        }

    def compile(self) -> CompiledView:
        """
        Compile descendant counts, graph choices, and samplers.
        """
        initial_state = self._initial_state()
        self._compile_from(initial_state)

        self._compile_choice_result_tuples()
        self._compile_choice_starts()
        samplers = self._make_samplers()

        return CompiledView(
            initial_state=initial_state,
            n_valid_sequences=self.totals_by_state[initial_state][0],
            choices_by_state=self.choices_by_state,
            choice_results_by_state=self.choice_results_by_state,
            choice_start_by_state=self.choice_start_by_state,
            codon_pos_by_state=self.codon_pos_by_state,
            fixed_choice_by_state=self.fixed_choice_by_state,
            samplers=samplers,
        )

    def _compile_from(self, initial_state: NodeState) -> None:
        """
        Walk the reachable graph states and compile each one after its children.
        """
        initial_node, initial_tracker_state, initial_constraint_state = initial_state
        stack = [(initial_node, initial_tracker_state, initial_constraint_state, False)]

        while stack:
            node, tracker_state, constraint_state, expanded = stack.pop()
            state = self._state(node, tracker_state, constraint_state)

            if state in self.totals_by_state:
                continue

            if node is self.graph.final_node:
                self._compile_final_state(state, constraint_state)
                continue

            if not expanded:
                stack.append((node, tracker_state, constraint_state, True))
                stack.extend(self._uncompiled_children(node, tracker_state, constraint_state))
                continue

            self._compile_state(node, tracker_state, constraint_state)

    def _compile_final_state(self, state: NodeState, constraint_state: ConstraintState) -> None:
        """
        Compile the final graph state.
        """
        if (
            not self.has_path_constraint
            or self.path_constraint.is_satisfied(constraint_state)
        ):
            self.totals_by_state[state] = (1, 0.0)
        else:
            self.totals_by_state[state] = (0, -math.inf)

        self.choices_by_state[state] = {}

    def _compile_state(
        self,
        node,
        tracker_state: TrackerState,
        constraint_state: ConstraintState,
    ) -> None:
        """
        Compile one non-final graph state after all valid children have been compiled.
        """
        state = self._state(node, tracker_state, constraint_state)
        raw_results = []

        self._record_state_kind(state, node)

        for choice in self.choices_by_node[node]:
            result = self._choice_result(node, tracker_state, constraint_state, choice)

            if result is None:
                continue

            raw_results.append(result)

        results = raw_results

        choice_results = {}
        descendant_count = 0
        descendant_weight_masses = []

        for result in results:
            choice_results[result.choice] = result
            descendant_count += result.descendant_count
            descendant_weight_masses.append(result.descendant_weight_mass)

        descendant_weight_mass = self._sum_weight_masses(descendant_weight_masses)

        self.choices_by_state[state] = choice_results
        self.totals_by_state[state] = (descendant_count, descendant_weight_mass)

    def _compile_choice_result_tuples(self) -> None:
        """
        Store each state's choice results as tuples for fast indexed traversal.
        """
        self.choice_results_by_state = {
            state: tuple(choice_results.values())
            for state, choice_results in self.choices_by_state.items()
        }

    def _get_choices_for_node(self, node) -> List[str]:
        """
        Return choices available to this node in this view.
        """
        if isinstance(node, CodonNode):
            if node.pos in self.view.pinned_codons:
                return self.view.pinned_codons[node.pos]

            return node.codons

        return [node.sequence]

    def _choice_result(
        self,
        node,
        tracker_state: TrackerState,
        constraint_state: ConstraintState,
        choice: str,
    ) -> Optional[ChoiceResult]:
        """
        Compile the result of taking one outgoing choice from a graph node.
        """
        child = node.transitions.get(choice)

        if child is None:
            return None

        advance = self._advance_tracker(tracker_state, node, choice)

        if advance.banned:
            return None

        if self.has_path_constraint:
            next_constraint_state = self.path_constraint.advance(
                constraint_state,
                node.pos,
                choice,
            )

            if next_constraint_state is None:
                return None
        else:
            next_constraint_state = ()

        child_state = self._state(child, advance.state, next_constraint_state)
        child_count, child_weight_mass = self.totals_by_state[child_state]

        if child_count == 0:
            return None

        descendant_weight_mass = self._choice_weight_mass(node, choice, child_weight_mass)

        if descendant_weight_mass == -math.inf:
            return None

        return ChoiceResult(
            choice=choice,
            descendant_count=child_count,
            descendant_weight_mass=descendant_weight_mass,
            next_state=None if child is self.graph.final_node else child_state,
            is_coding=isinstance(node, CodonNode),
        )

    def _uncompiled_children(
        self,
        node,
        tracker_state: TrackerState,
        constraint_state: ConstraintState,
    ) -> List[Tuple[object, TrackerState, ConstraintState, bool]]:
        """
        Return uncompiled child states reachable from a graph state.
        """
        children = []

        for choice in self.choices_by_node[node]:
            child = node.transitions.get(choice)

            if child is None:
                continue

            advance = self._advance_tracker(tracker_state, node, choice)

            if advance.banned:
                continue

            if self.has_path_constraint:
                next_constraint_state = self.path_constraint.advance(
                    constraint_state,
                    node.pos,
                    choice,
                )

                if next_constraint_state is None:
                    continue
            else:
                next_constraint_state = ()

            child_state = self._state(child, advance.state, next_constraint_state)

            if child_state not in self.totals_by_state:
                children.append((child, advance.state, next_constraint_state, False))

        return children

    def _normalise_choice_results(self, results: List[ChoiceResult]) -> List[ChoiceResult]:
        """
        Rescale descendant weight masses within one state.

        Only relative weights matter for sampling, so this prevents long paths
        from underflowing toward zero.
        """
        max_mass = max((result.descendant_weight_mass for result in results), default=0.0)

        if max_mass <= 0:
            return results

        return [
            ChoiceResult(
                choice=result.choice,
                descendant_count=result.descendant_count,
                descendant_weight_mass=result.descendant_weight_mass / max_mass,
                next_state=result.next_state,
                is_coding=result.is_coding,
            )
            for result in results
        ]

    def _normalise_weights(self, weights: List[float]) -> List[float]:
        """
        Rescale sampler weights defensively.
        """
        if not weights:
            return weights

        max_weight = max(weights)

        if max_weight == -math.inf:
            return [1.0] * len(weights)

        return [math.exp(weight - max_weight) for weight in weights]

    def _make_samplers(self) -> dict:
        """
        Make samplers for each reachable graph state.
        """
        samplers = {}

        for state, choice_results in self.choice_results_by_state.items():
            node, _, _ = state

            if node is self.graph.final_node:
                continue

            runtime_items = []
            runtime_weights = []

            for result in choice_results:
                runtime_items.append((result.choice, result.is_coding, result.next_state))
                runtime_weights.append(result.descendant_weight_mass)

            if runtime_items:
                runtime_weights = self._normalise_weights(runtime_weights)
                samplers[state] = Sampler(runtime_items, runtime_weights, rng=self.view._rng)

        return samplers

    def _record_state_kind(self, state: NodeState, node) -> None:
        """
        Record whether a graph state consumes a codon from the user sequence
        or follows a fixed context sequence.
        """
        if isinstance(node, CodonNode):
            self.codon_pos_by_state[state] = node.pos
        else:
            self.fixed_choice_by_state[state] = node.sequence

    def _choice_weight_mass(self, node, choice: str, child_weight_mass: float) -> float:
        """
        Return the weighted mass for a choice.
        """
        if isinstance(node, CodonNode):
            weight = self.graph.cw[choice]

            if weight <= 0:
                return -math.inf

            return math.log(weight) + child_weight_mass

        return child_weight_mass

    def _sum_weight_masses(self, weights: List[float]) -> float:
        """
        Sum weight masses defensively.
        """
        weights = [weight for weight in weights if weight != -math.inf]

        if not weights:
            return -math.inf

        max_weight = max(weights)

        return max_weight + math.log(
            sum(math.exp(weight - max_weight) for weight in weights)
        )

    def _initial_state(self) -> NodeState:
        """
        Return the initial compiled graph state.
        """
        if self.has_path_constraint:
            constraint_state = self.path_constraint.initial_state
        else:
            constraint_state = ()

        return self._state(
            self.graph.initial_node,
            self._initial_tracker_state(),
            constraint_state,
        )

    def _initial_tracker_state(self) -> TrackerState:
        """
        Return the initial banned-sequence tracker state.
        """
        if not self.view.banned_sequences:
            return frozenset()

        return self.tracker.initial_state

    def _compile_choice_starts(self) -> None:
        """
        Store sequence slice starts for coding states.
        """
        self.choice_start_by_state = {
            state: (pos - 1) * 3
            for state, pos in self.codon_pos_by_state.items()
        }

    def _advance_tracker(self, tracker_state: TrackerState, node: Node, choice: str) -> AdvanceResult:
        """
        Advance banned-sequence tracking after taking a graph step. Results are cached.
        """
        key = (node, tracker_state, choice)

        if key in self.advance_cache:
            return self.advance_cache[key]

        if self.tracker.is_trivial:
            result = AdvanceResult(banned=False, state=tracker_state)
        else:
            step = (node.pos, choice)
            result = self.tracker.advance(step, tracker_state)

        self.advance_cache[key] = result
        return result

    def _state(
            self,
            node,
            tracker_state: TrackerState = frozenset(),
            constraint_state: ConstraintState = (),
    ) -> NodeState:
        """
        Return the compiled state for a graph node plus tracker states.
        """
        return node, tracker_state, constraint_state
