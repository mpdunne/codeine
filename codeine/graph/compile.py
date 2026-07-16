import math

from typing import Dict, List, NamedTuple, Optional, Tuple, TYPE_CHECKING

from codeine.constraints.base import ConstraintState, DEAD_STATE, SAFE_STATE
from codeine.graph.nodes import CodonNode, Node

if TYPE_CHECKING:
    from codeine.graph.view import CodonGraphView


# The traversal state consists of the current graph node and the current state
# of each active constraint. Plain tuples are immutable and hashable, while
# avoiding NamedTuple construction in the compiler's hottest path.
TraversalState = Tuple[Node, Tuple[ConstraintState, ...]]
TraversalStateKey = Tuple[int, Tuple[ConstraintState, ...]]


class ChoiceResult(NamedTuple):
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
    next_state_id: Optional[int]
    is_coding: bool


class CompiledView(NamedTuple):
    """
    Cached data for a compiled CodonGraphView, to speed up sampling and enumeration.
    """
    initial_state: TraversalState
    initial_state_id: int
    states: Tuple[TraversalState, ...]

    # Deep compiled transitions:
    # state ID -> (choice, child state ID)
    # These include every transition allowed by the graph and constraints,
    # before temporary view pins are applied.
    child_results_by_state_id: Tuple[Tuple[Tuple[str, int], ...], ...]

    # Compiled graph choices (lookup):
    # state ID -> choice -> ChoiceResult
    # Used for fast sequence validation and graph traversal in the graph view.
    choices_by_state_id: Tuple[Dict[str, ChoiceResult], ...]

    # Compiled graph choices (iteration):
    # state ID -> ChoiceResults in graph order
    # Used for fast sampling and sequence enumeration in the graph view.
    choice_results_by_state_id: Tuple[Tuple[ChoiceResult, ...], ...]

    n_valid_sequences: int


class ViewCompiler:
    """
    Compile a CodonGraphView into cached choice, count, and sampling data.
    """

    def __init__(self, view: 'CodonGraphView') -> None:
        self.view = view
        self.graph = view.graph

        constraints = view.constraints

        for constraint in constraints:
            constraint.link(self.graph)

        self.constraints = tuple(constraint for constraint in constraints if not constraint.is_trivial)
        self.constraint_advancers = tuple(constraint.advance for constraint in self.constraints)

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
        # Used for fast sequence validation and graph traversal in the graph view.
        self.choices_by_state_id: List[Optional[Dict[str, ChoiceResult]]] = []

        # The cached log-ified codon weights, to avoid repeated log calculations.
        self.log_codon_weights = {
            codon: math.log(weight) if weight > 0 else -math.inf
            for codon, weight in self.view.codon_weights.weights.items()
        }

        # Graph transitions available at each node. Permanent graph-level codon
        # restrictions are already reflected in node.transitions. Temporary view
        # pins are applied later, during the shallow calculation pass.
        self.transitions_by_node = {
            node: tuple(node.transitions.items())
            for node in self.graph.nodes
            if node is not self.graph.final_node
        }

    def compile(self) -> CompiledView:
        """
        Compile descendant counts, graph choices, and sampling masses.

        Returns
        -------
        CompiledView
            A compiled view.
        """
        initial_state = self._initial_state()
        initial_node, initial_constraint_states = initial_state
        initial_state_id, _ = self._get_or_register_state_id(
            initial_node,
            initial_constraint_states,
        )

        self._compile_topology(initial_state_id)
        self._calculate_results()

        child_results_by_state_id = tuple(
            tuple(results or ())
            for results in self.child_results_by_state_id
        )

        choices_by_state_id = tuple(
            choices or {}
            for choices in self.choices_by_state_id
        )

        choice_results_by_state_id = tuple(
            tuple(choices.values()) if choices else ()
            for choices in self.choices_by_state_id
        )

        initial_total = self.totals_by_state_id[initial_state_id]
        assert initial_total is not None

        return CompiledView(
            initial_state=initial_state,
            initial_state_id=initial_state_id,
            states=tuple(self.states),
            child_results_by_state_id=child_results_by_state_id,
            choices_by_state_id=choices_by_state_id,
            choice_results_by_state_id=choice_results_by_state_id,
            n_valid_sequences=initial_total[0],
        )

    def _get_or_register_state_id(
            self,
            node: Node,
            constraint_states: Tuple[ConstraintState, ...],
    ) -> Tuple[int, bool]:
        """
        Return the state ID and whether the state was newly registered.

        State lookup uses the node position rather than the node object, keeping
        the hot dictionary key compact. The full node is stored only when a state
        is genuinely new.
        """
        key = (node.pos, constraint_states)
        state_id = self.state_ids.get(key)

        if state_id is not None:
            return state_id, False

        state_id = len(self.states)
        self.state_ids[key] = state_id
        self.states.append((node, constraint_states))
        self.totals_by_state_id.append(None)
        self.choices_by_state_id.append(None)
        self.child_results_by_state_id.append(None)

        return state_id, True

    def _initial_state(self) -> TraversalState:
        """
        Return the starting traversal state.

        The initial state starts at the graph's initial node, with fresh
        constraint states.

        Returns
        -------
        TraversalState
            Starting state for graph compilation.
        """
        constraint_states = tuple(constraint.initial_state for constraint in self.constraints)

        return self.graph.initial_node, constraint_states

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

            node, _constraint_states = self.states[state_id]

            if node is self.graph.final_node:
                self.child_results_by_state_id[state_id] = []
                continue

            stack.extend(child_id for child_id, _expanded in self._uncompiled_children(state_id))

    def _calculate_results(self) -> None:
        """
        Resolve descendant counts, probability masses, and choice results.

        Temporary view pins and codon weights are applied during this pass.

        States are processed from right to left through the graph, ensuring
        that every child has been resolved before its parent.
        """
        state_ids_by_node = {node: [] for node in self.graph.nodes}

        for state_id, (node, _constraint_states) in enumerate(self.states):
            state_ids_by_node[node].append(state_id)

        for node in reversed(self.graph.nodes):
            for state_id in state_ids_by_node[node]:
                if node is self.graph.final_node:
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
        node to produce a ChoiceResult. The total descendant count and log mass
        are then cached for the current state.

        Parameters
        ----------
        state_id
            ID of the traversal state being compiled.
        """
        state = self.states[state_id]
        node, _constraint_states = state

        choice_results = {}
        descendant_count = 0

        max_log_mass = -math.inf
        relative_mass_sum = 0.0

        is_coding = isinstance(node, CodonNode)
        child_results = self.child_results_by_state_id[state_id] or ()

        pinned_codons = (self.view.pinned_codons.get(node.pos) if is_coding else None)

        for choice, child_id in child_results:

            if pinned_codons is not None and choice not in pinned_codons:
                continue

            child, _child_constraint_states = self.states[child_id]
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
                next_state_id= None if child is self.graph.final_node else child_id,
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

    def _uncompiled_children(self, state_id: int) -> List[Tuple[int, bool]]:
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
        list of tuple
            Stack entries for child state IDs still needing compilation.
        """
        child_results = []
        uncompiled_children = []

        node, constraint_states = self.states[state_id]
        pos = node.pos

        for choice, child in self.transitions_by_node[node]:
            next_constraint_states = self._advance_constraints(constraint_states, pos, choice)

            if next_constraint_states is None:
                continue

            child_id, is_new = self._get_or_register_state_id(child, next_constraint_states)

            child_results.append((choice, child_id))

            if is_new:
                uncompiled_children.append((child_id, False))

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
            Position of the current graph node.
        choice
            Graph choice taken from the current node.

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