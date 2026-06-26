import math

from dataclasses import dataclass
from typing import Dict, NamedTuple, Tuple, List, Optional, TYPE_CHECKING

from codeine.constraints.banned import BannedTrackerState, AdvanceResult
from codeine.constraints.base import ConstraintState
from codeine.graph.nodes import CodonNode, Node, ContextNode
from codeine.utils.sampling import Sampler

if TYPE_CHECKING:
    from codeine.graph.view import CodonGraphView


class TraversalState(NamedTuple):
    """
    The traversal state consists of the current node, plus a summary of the relevant
    parts of how we got there. For example, it tracks whether we have seen parts of
    banned sequences, and can track nucleotide or codon properties.

    Different graph traversal histories that produce the same traversal state are
    collapsed and are equivalent under this framework.
    """
    node: Node
    banned_tracker_state: BannedTrackerState
    constraint_state: ConstraintState


@dataclass(frozen=True)
class ChoiceResult:
    """
    Cached result of taking one graph choice from one compiled state. The
    "choice" is the graph edge label, i.e. a codon or a context sequence.

    Each ChoiceResult is specific to its location in the graph. The descendant
    counts and log mass are calculated iteratively by summing the values of
    downstream nodes.
    """
    choice: str
    descendant_count: int
    descendant_log_mass: float
    next_state: Optional[TraversalState]
    next_state_id: Optional[int]
    is_coding: bool


@dataclass(frozen=True)
class CompiledView:
    """
    Cached data for a compiled CodonGraphView, to speed up sampling and enumeration.
    """
    initial_state: TraversalState
    initial_state_id: int
    states: Tuple[TraversalState, ...]

    # Compiled graph choices (lookup):
    # state/ID -> choice -> ChoiceResult
    # Used for fast sequence validation and graph traversal in the graph view.
    choices_by_state: Dict[TraversalState, Dict[str, ChoiceResult]]
    choices_by_state_id: Tuple[Dict[str, ChoiceResult], ...]

    # Compiled graph choices (iteration):
    # state/ID -> ChoiceResults in graph order
    # Used for fast sampling and sequence enumeration in the graph view.
    choice_results_by_state: Dict[TraversalState, Tuple[ChoiceResult, ...]]
    choice_results_by_state_id: Tuple[Tuple[ChoiceResult, ...], ...]

    n_valid_sequences: int
    samplers: dict
    samplers_by_state_id: tuple


class ViewCompiler:
    """
    Compile a CodonGraphView into cached choice, count, and sampling data.
    """

    def __init__(self, view: 'CodonGraphView') -> None:
        self.view = view
        self.graph = view.graph
        self.banned_tracker = view.banned_tracker
        self.path_constraint = view.path_constraint
        self.has_path_constraint = self.path_constraint is not None

        self.state_ids: Dict[TraversalState, int] = {}
        self.states: List[TraversalState] = []

        # Dynamic-programming totals:
        # state -> (descendant count, descendant log mass)
        # Avoids repeatedly recomputing subtree sizes and probability masses.
        self.totals_by_state: Dict[TraversalState, Tuple[int, float]] = {}
        self.totals_by_state_id: List[Optional[Tuple[int, float]]] = []

        # Banned-sequence tracker transitions:
        # (node, tracker state, choice) -> tracker result
        # Avoids recomputing tracker advances during compilation.
        self.banned_advance_cache: Dict[Tuple[Node, BannedTrackerState, str], AdvanceResult] = {}

        # Traversal transition table:
        # (state, choice) -> child state
        # Avoids rediscovering successor states during compilation.
        self.child_state_by_state_choice: Dict[Tuple[TraversalState, str], TraversalState] = {}

        # Compiled graph choices (lookup):
        # state -> choice -> ChoiceResult
        # Used for fast sequence validation and graph traversal in the graph view.
        self.choices_by_state: Dict[TraversalState, Dict[str, ChoiceResult]] = {}

        # Compiled graph choices (iteration):
        # state -> ChoiceResults in graph order
        # Used for fast sampling and sequence enumeration in the graph view.
        self.choice_results_by_state: Dict[TraversalState, Tuple[ChoiceResult, ...]] = {}

        # The cached log-ified codon weights, to avoid repeated log calculations
        self.log_codon_weights = {
            codon: math.log(weight)
            for codon, weight in self.graph.cw.weights.items()
            if weight > 0
        }

        # Cached version of the choices available at each node, taking into account fixed codons & pins
        self.choices_by_node = {
            node: tuple(self._get_choices_for_node(node))
            for node in self.graph.nodes
            if node is not self.graph.final_node
        }

    def compile(self) -> CompiledView:
        """
        Compile descendant counts, graph choices, and samplers.

        Returns
        -------
        CompiledView
            A compiled view.
        """
        initial_state = self._initial_state()
        initial_state_id = self._get_or_register_state_id(initial_state)

        self._compile_from(initial_state)

        self.choice_results_by_state = {
            state: tuple(choice_results.values())
            for state, choice_results in self.choices_by_state.items()
        }

        choices_by_state_id = tuple(
            self.choices_by_state.get(state, {})
            for state in self.states
        )

        choice_results_by_state_id = tuple(
            self.choice_results_by_state.get(state, ())
            for state in self.states
        )

        samplers = self._make_samplers()

        samplers_by_state_id = tuple(
            samplers.get(state)
            for state in self.states
        )

        return CompiledView(
            initial_state=initial_state,
            initial_state_id=initial_state_id,
            states=tuple(self.states),
            n_valid_sequences=self.totals_by_state[initial_state][0],
            choices_by_state=self.choices_by_state,
            choice_results_by_state=self.choice_results_by_state,
            choices_by_state_id=choices_by_state_id,
            choice_results_by_state_id=choice_results_by_state_id,
            samplers=samplers,
            samplers_by_state_id=samplers_by_state_id,
        )

    def _get_or_register_state_id(self, state: TraversalState) -> int:
        """
        Return the stable integer ID for a traversal state, creating one if needed.
        """
        state_id = self.state_ids.get(state)

        if state_id is None:
            state_id = len(self.states)
            self.state_ids[state] = state_id
            self.states.append(state)
            self.totals_by_state_id.append(None)

        return state_id

    def _initial_state(self) -> TraversalState:
        """
        Return the starting traversal state.

        The initial state starts at the graph's initial node, with fresh banned-
        sequence and path-constraint states if those systems are active.

        Returns
        -------
        TraversalState
            Starting state for graph compilation.
        """
        if self.has_path_constraint:
            constraint_state = self.path_constraint.initial_state
        else:
            constraint_state = ()

        if self.view.banned_sequences:
            banned_tracker_state = self.banned_tracker.initial_state
        else:
            banned_tracker_state = frozenset()

        return TraversalState(self.graph.initial_node, banned_tracker_state, constraint_state)

    def _compile_from(self, initial_state: TraversalState) -> None:
        """
        Compile every reachable traversal state starting from an initial state.

        Uses an explicit depth-first stack so that each non-final state is compiled
        only after its child states have been compiled.

        Parameters
        ----------
        initial_state
            State from which graph compilation should begin.
        """
        stack = [(initial_state, False)]

        while stack:
            state, expanded = stack.pop()
            node = state.node

            _ = self._get_or_register_state_id(state)

            state_id = self.state_ids[state]

            if self.totals_by_state_id[state_id] is not None:
                continue

            if node is self.graph.final_node:
                self._compile_final_state(state)
                continue

            if not expanded:
                stack.append((state, True))
                stack.extend(self._uncompiled_children(state))
                continue

            self._compile_state(state)

    def _compile_final_state(self, state: TraversalState) -> None:
        """
        Compile a terminal traversal state.

        By the time a terminal state is reached, all graph choices have already
        been processed, including the right context. Choices rejected by the
        banned-sequence tracker or path constraint would not have reached this
        state.

        The only remaining decision is whether the final path-constraint state is
        acceptable. If so, this terminal state contributes one complete sequence
        with log mass 0.0. Otherwise, it contributes no sequences.

        Parameters
        ----------
        state
            Terminal traversal state being compiled.
        """
        state_id = self.state_ids[state]

        if not self.has_path_constraint or self.path_constraint.is_satisfied(state.constraint_state):
            total = (1, 0.0)
        else:
            total = (0, -math.inf)

        self.totals_by_state[state] = total
        self.totals_by_state_id[state_id] = total
        self.choices_by_state[state] = {}

    def _compile_state(self, state: TraversalState) -> None:
        """
        Compile one non-final traversal state.

        For each outgoing graph choice, combine the previously compiled child state
        with the contribution from the current node to produce a ChoiceResult.
        The total descendant count and log mass are then cached for the current state.

        Parameters
        ----------
        state
            Traversal state being compiled.
        """
        node = state.node
        choice_results = {}
        descendant_count = 0
        descendant_log_masses = []
        is_coding = isinstance(node, CodonNode)

        for choice in self.choices_by_node[node]:
            child_state = self.child_state_by_state_choice.get((state, choice))

            if child_state is None:
                continue

            child = child_state.node
            child_id = self.state_ids[child_state]
            child_total = self.totals_by_state_id[child_id]

            if child_total is None:
                continue

            child_count, subtree_log_mass = child_total

            if child_count == 0:
                continue

            choice_log_mass = self._accumulate_log_mass(node, choice, subtree_log_mass)

            if choice_log_mass == -math.inf:
                continue

            result = ChoiceResult(
                choice=choice,
                descendant_count=child_count,
                descendant_log_mass=choice_log_mass,
                next_state=None if child is self.graph.final_node else child_state,
                next_state_id=None if child is self.graph.final_node else self.state_ids[child_state],
                is_coding=is_coding,
            )

            choice_results[choice] = result
            descendant_count += child_count
            descendant_log_masses.append(choice_log_mass)

        descendant_log_mass = self._sum_log_masses(descendant_log_masses)

        state_id = self.state_ids[state]
        total = (descendant_count, descendant_log_mass)

        self.choices_by_state[state] = choice_results
        self.totals_by_state[state] = total
        self.totals_by_state_id[state_id] = total

    def _uncompiled_children(self, state: TraversalState) -> List[Tuple[TraversalState, bool]]:
        """
        Return child states reached by taking each outgoing graph choice.

        Choices rejected by the banned-sequence tracker or path constraint are skipped.
        Only child states that have not yet been compiled are returned.

        Parameters
        ----------
        state
            Traversal state whose children should be discovered.

        Returns
        -------
        list of tuple
            Stack entries for child states still needing compilation.
        """
        children = []
        node = state.node

        for choice in self.choices_by_node[node]:
            child = node.transitions.get(choice)

            if child is None:
                continue

            advance = self._advance_banned_tracker(state.banned_tracker_state, node, choice)

            if advance.banned:
                continue

            if self.has_path_constraint:
                next_constraint_state = self.path_constraint.advance(state.constraint_state, node.pos, choice)

                if next_constraint_state is None:
                    continue
            else:
                next_constraint_state = ()

            child_state = TraversalState(child, advance.state, next_constraint_state)
            child_id = self._get_or_register_state_id(child_state)
            self.child_state_by_state_choice[(state, choice)] = child_state

            if self.totals_by_state_id[child_id] is None:
                children.append((child_state, False))

        return children

    def _make_samplers(self) -> dict:
        """
        Build weighted samplers for every compiled traversal state.

        Each sampler chooses between the state's valid outgoing choices with
        probabilities proportional to their descendant probability masses.

        Returns
        -------
        dict
            Mapping from traversal state to a sampler over its outgoing choices.
        """
        samplers = {}

        for state, choice_results in self.choice_results_by_state.items():
            node = state.node

            if node is self.graph.final_node:
                continue

            runtime_items = []
            runtime_log_masses = []

            for result in choice_results:
                runtime_items.append((result.choice, result.is_coding, result.next_state_id))
                runtime_log_masses.append(result.descendant_log_mass)

            if runtime_items:
                runtime_weights = self._convert_log_masses_to_sampler_weights(runtime_log_masses)
                samplers[state] = Sampler(runtime_items, runtime_weights, rng=self.view._rng)

        return samplers

    def _get_choices_for_node(self, node: Node) -> List[str]:
        """
        Return the graph choices available from a node in this view.

        Codon nodes respect any pinned codons defined by the view. Context nodes
        always have a single fixed sequence.

        Parameters
        ----------
        node
            Graph node whose available choices are required.

        Returns
        -------
        list of str
            The choices that may be taken from this node in the current view.
        """
        if isinstance(node, CodonNode):
            if node.pos in self.view.pinned_codons:
                return self.view.pinned_codons[node.pos]
            else:
                return node.codons

        elif isinstance(node, ContextNode):
            return [node.sequence]

    def _advance_banned_tracker(
            self,
            banned_tracker_state: BannedTrackerState,
            node: Node,
            choice: str,
    ) -> AdvanceResult:
        """
        Advance the banned-sequence tracker after taking one graph choice.

        Results are cached because the same tracker transition may be encountered
        from many traversal states during compilation.

        Parameters
        ----------
        banned_tracker_state
            Current banned-sequence tracker state.
        node
            Current graph node.
        choice
            Graph choice taken from the current node.

        Returns
        -------
        AdvanceResult
            Whether the choice enters a banned state and the resulting tracker
            state.
        """
        key = (node, banned_tracker_state, choice)

        if key in self.banned_advance_cache:
            return self.banned_advance_cache[key]

        if self.banned_tracker.is_trivial:
            result = AdvanceResult(banned=False, state=banned_tracker_state)
        else:
            step = (node.pos, choice)
            result = self.banned_tracker.advance(step, banned_tracker_state)

        self.banned_advance_cache[key] = result
        return result

    def _accumulate_log_mass(
            self,
            node: Node,
            choice: str,
            subtree_log_mass: float,
    ) -> float:
        """
        Accumulate the log probability mass contributed by one graph choice.

        The subtree log mass has already been computed for the child state. Codon
        nodes contribute the log of their codon weight, whereas context nodes
        contribute no additional mass.

        Parameters
        ----------
        node
            The graph node from which the choice is taken.
        choice
            The outgoing graph choice.
        subtree_log_mass
            The total log mass reachable from the child state.

        Returns
        -------
        float
            The total log mass reachable after taking this choice.
        """
        if isinstance(node, CodonNode):
            codon_log_weight = self.log_codon_weights.get(choice)

            if codon_log_weight is None:
                return -math.inf

            return codon_log_weight + subtree_log_mass

        return subtree_log_mass

    def _sum_log_masses(self, log_masses: List[float]) -> float:
        """
        Combine several subtree log masses into a single log mass.

        The calculation is performed using the log-sum-exp trick to avoid numerical
        underflow when the subtree probabilities are extremely small.

        Parameters
        ----------
        log_masses
            Log-space masses to sum.

        Returns
        -------
        float
            The log of the summed masses, or -inf if no finite masses exist.
        """
        log_masses = [m for m in log_masses if m != -math.inf]

        if not log_masses:
            return -math.inf

        max_log_mass = max(log_masses)

        total_relative_mass = sum(math.exp(log_mass - max_log_mass) for log_mass in log_masses)

        return max_log_mass + math.log(total_relative_mass)

    def _convert_log_masses_to_sampler_weights(self, log_masses: List[float]) -> List[float]:
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
        if not log_masses:
            return log_masses

        max_log_mass = max(log_masses)

        if max_log_mass == -math.inf:
            return [1.0] * len(log_masses)

        return [math.exp(log_mass - max_log_mass) for log_mass in log_masses]
