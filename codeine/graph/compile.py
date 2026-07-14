import math

from typing import Dict, NamedTuple, Tuple, List, Optional, TYPE_CHECKING

from codeine.constraints.base import ConstraintState, DEAD_STATE
from codeine.graph.nodes import CodonNode, Node, ContextNode

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

        self.state_ids: Dict[TraversalStateKey, int] = {}
        self.states: List[TraversalState] = []

        # Dynamic-programming totals:
        # state ID -> (descendant count, descendant log mass)
        # Avoids repeatedly recomputing subtree sizes and probability masses.
        self.totals_by_state_id: List[Optional[Tuple[int, float]]] = []

        # Traversal transition table:
        # state ID -> [(choice, child state ID), ...]
        # Avoids rediscovering successor states during compilation.
        self.child_results_by_state_id: List[
            Optional[List[Tuple[str, int]]]
        ] = []

        # Compiled graph choices (lookup):
        # state ID -> choice -> ChoiceResult
        # Used for fast sequence validation and graph traversal in the graph view.
        self.choices_by_state_id: List[
            Optional[Dict[str, ChoiceResult]]
        ] = []

        # The cached log-ified codon weights, to avoid repeated log calculations
        self.log_codon_weights = {
            codon: math.log(weight)
            for codon, weight in self.graph.cw.weights.items()
            if weight > 0
        }

        # Cached version of the choices available at each node, taking into
        # account fixed codons & pins.
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
        initial_node, initial_constraint_states = initial_state
        initial_state_id, _ = self._get_or_register_state_id(
            initial_node,
            initial_constraint_states,
        )

        self._compile_from(initial_state_id)

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
            n_valid_sequences=initial_total[0],
            choices_by_state_id=choices_by_state_id,
            choice_results_by_state_id=choice_results_by_state_id,
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

        The initial state starts at the graph's initial node, with fresh constraint states.

        Returns
        -------
        TraversalState
            Starting state for graph compilation.
        """
        constraint_states = tuple(constraint.initial_state for constraint in self.constraints)

        return (
            self.graph.initial_node,
            constraint_states,
        )

    def _compile_from(self, initial_state_id: int) -> None:
        """
        Compile every reachable traversal state starting from an initial state ID.

        Uses an explicit depth-first stack so that each non-final state is compiled
        only after its child states have been compiled.

        Parameters
        ----------
        initial_state_id
            ID of the state from which graph compilation should begin.
        """
        stack = [(initial_state_id, False)]

        while stack:
            state_id, expanded = stack.pop()
            state = self.states[state_id]
            node, _constraint_states = state

            if self.totals_by_state_id[state_id] is not None:
                continue

            if node is self.graph.final_node:
                self._compile_final_state(state_id)
                continue

            if not expanded:
                stack.append((state_id, True))
                stack.extend(self._uncompiled_children(state_id))
                continue

            self._compile_state(state_id)

    def _compile_final_state(self, state_id: int) -> None:
        """
        Compile a terminal traversal state.

        By the time a terminal state is reached, all graph choices have already
        been processed, including the right context. Choices rejected by the constraint
        would not have reached this state.

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

        For each outgoing graph choice, combine the previously compiled child state
        with the contribution from the current node to produce a ChoiceResult.
        The total descendant count and log mass are then cached for the current state.

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

        for choice, child_id in child_results:
            child, _child_constraint_states = self.states[child_id]
            child_total = self.totals_by_state_id[child_id]

            if child_total is None:
                continue

            child_count, subtree_log_mass = child_total

            if child_count == 0:
                continue

            choice_log_mass = self._accumulate_log_mass(
                node,
                choice,
                subtree_log_mass,
            )

            if choice_log_mass == -math.inf:
                continue

            result = ChoiceResult(
                choice=choice,
                descendant_count=child_count,
                descendant_log_mass=choice_log_mass,
                next_state_id=None if child is self.graph.final_node else child_id,
                is_coding=is_coding,
            )

            choice_results[choice] = result
            descendant_count += child_count

            # Incremental log-sum-exp.
            if choice_log_mass <= max_log_mass:
                relative_mass_sum += math.exp(choice_log_mass - max_log_mass)
            else:
                if max_log_mass == -math.inf:
                    relative_mass_sum = 1.0
                else:
                    relative_mass_sum = (
                            relative_mass_sum
                            * math.exp(max_log_mass - choice_log_mass)
                            + 1.0
                    )

                max_log_mass = choice_log_mass

        if max_log_mass == -math.inf:
            descendant_log_mass = -math.inf
        else:
            descendant_log_mass = max_log_mass + math.log(relative_mass_sum)

        self.choices_by_state_id[state_id] = choice_results
        self.totals_by_state_id[state_id] = (
            descendant_count,
            descendant_log_mass,
        )

    def _uncompiled_children(self, state_id: int) -> List[Tuple[int, bool]]:
        """
        Return child state IDs reached by taking each outgoing graph choice.

        Choices rejected by constraints are skipped. Only child states that have not yet been compiled are returned.

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

        for choice in self.choices_by_node[node]:
            child = node.transitions.get(choice)

            if child is None:
                continue

            next_constraint_states = self._advance_constraints(
                constraint_states,
                pos,
                choice,
            )

            if next_constraint_states is None:
                continue

            child_id, is_new = self._get_or_register_state_id(
                child,
                next_constraint_states,
            )

            child_results.append((choice, child_id))

            if is_new:
                uncompiled_children.append((child_id, False))

        self.child_results_by_state_id[state_id] = child_results

        return uncompiled_children

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

        for constraint, state in zip(self.constraints, constraint_states):
            next_state = constraint.advance(state, pos, choice)
            if next_state == DEAD_STATE:
                return None

            next_states.append(next_state)

        return tuple(next_states)

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
