import math

from typing import Dict, List, NamedTuple, Optional, Sequence, Tuple, TYPE_CHECKING

from codeine.constraints.base import Constraint, ConstraintState, DEAD_STATE, SAFE_STATE

if TYPE_CHECKING:
    from codeine.graph.view import CodonGraphView


# The traversal state consists of the current graph position
# and the current state of each active constraint.
ConstraintStates = Tuple[ConstraintState, ...]
ConstraintStateId = int
TraversalState = Tuple[int, ConstraintStateId]
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


class CompiledView(NamedTuple):
    """
    Cached data for a compiled CodonGraphView, to speed up sampling and enumeration.
    """
    initial_state: TraversalState
    initial_state_id: int
    states: Tuple[TraversalState, ...]
    constraint_states: Tuple[ConstraintStates, ...]

    # Deep compiled transitions:
    # state ID -> ((choice, child state ID), ...)
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
    Compile a CodonGraphView into cached topology, choice, count, and sampling-mass data.
    """

    def __init__(self, view: 'CodonGraphView') -> None:
        self.view = view
        self.graph = view.graph

        self.constraints: Tuple[Constraint, ...] = ()
        self.constraint_advancers = ()

        self.constraint_state_ids: Dict[ConstraintStates, ConstraintStateId] = {}
        self.constraint_states: List[ConstraintStates] = []

        self.state_ids: Dict[TraversalStateKey, int] = {}
        self.states: List[TraversalState] = []

        # Dynamic-programming totals:
        # state ID -> descendant count
        # state ID -> descendant log mass
        # Avoids repeatedly recomputing subtree sizes and probability masses.
        self.descendant_counts: List[Optional[int]] = []
        self.descendant_log_masses: List[Optional[float]] = []

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

    def compile(self) -> CompiledView:
        """
        Compile descendant counts, graph choices, and sampling masses.

        Returns
        -------
        CompiledView
            A compiled view.
        """
        self._set_constraints(self.view.constraints)

        initial_pos, initial_constraint_states = self._initial_state()
        initial_state_id, _ = self._get_or_register_state_id(
            initial_pos,
            initial_constraint_states,
        )

        self._compile_topology(initial_state_id)
        self._compile_choices()

        return self._compiled_view(initial_state_id)

    def compile_shallow(self, compiled: CompiledView) -> CompiledView:
        """
        Recompile choices, counts, and probability masses using an existing
        deep topology.

        Parameters
        ----------
        compiled
            The existing compiled view whose states and transitions should be reused.

        Returns
        -------
        CompiledView
            The compiled view with updated shallow data.
        """
        self.states = list(compiled.states)
        self.constraint_states = list(compiled.constraint_states)
        self.child_results_by_state_id = list(compiled.child_results_by_state_id)

        self.descendant_counts = [None] * len(self.states)
        self.descendant_log_masses = [None] * len(self.states)
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

        initial_count = self.descendant_counts[compiled.initial_state_id]
        assert initial_count is not None

        return compiled._replace(
            choices_by_state_id=choices_by_state_id,
            choice_results_by_state_id=choice_results_by_state_id,
            n_valid_sequences=initial_count,
        )

    def extend(
            self,
            compiled: CompiledView,
            constraints: Sequence[Constraint],
    ) -> CompiledView:
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
        CompiledView
            An updated compiled view.
        """
        self._set_constraints(constraints)

        if not self.constraints:
            return self.compile_shallow(compiled)

        initial_new_states = tuple(constraint.initial_state for constraint in self.constraints)
        initial_pos, old_initial_constraint_state_id = compiled.initial_state
        old_initial_states = compiled.constraint_states[old_initial_constraint_state_id]
        initial_state_id, _ = self._get_or_register_state_id(
            initial_pos,
            old_initial_states + initial_new_states,
        )

        self._compile_extended_topology(compiled, initial_state_id)
        self._compile_choices()

        return self._compiled_view(initial_state_id)

    def _get_or_register_constraint_state_id(
            self,
            constraint_states: ConstraintStates,
    ) -> ConstraintStateId:
        """
        Return the dense ID for a tuple of constraint states, registering it if
        necessary.
        """
        constraint_state_id = self.constraint_state_ids.get(constraint_states)

        if constraint_state_id is not None:
            return constraint_state_id

        constraint_state_id = len(self.constraint_states)
        self.constraint_state_ids[constraint_states] = constraint_state_id
        self.constraint_states.append(constraint_states)

        return constraint_state_id

    def _get_or_register_state_id(
            self,
            pos: int,
            constraint_states: ConstraintStates,
    ) -> Tuple[int, bool]:
        """
        Return the state ID and whether the state was newly registered.

        Traversal-state lookup uses only the graph position and a dense integer
        ID for the full constraint-state tuple.
        """
        constraint_state_id = self._get_or_register_constraint_state_id(
            constraint_states,
        )
        key = (pos, constraint_state_id)
        state_id = self.state_ids.get(key)

        if state_id is not None:
            return state_id, False

        state_id = len(self.states)
        self.state_ids[key] = state_id
        self.states.append(key)
        self.descendant_counts.append(None)
        self.descendant_log_masses.append(None)
        self.choices_by_state_id.append(None)
        self.child_results_by_state_id.append(None)

        return state_id, True

    def _initial_state(self) -> Tuple[int, ConstraintStates]:
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

            pos, _constraint_state_id = self.states[state_id]

            if pos == self.final_pos:
                self.child_results_by_state_id[state_id] = []
                continue

            stack.extend(child_id for child_id, _expanded in self._uncompiled_children(state_id))

    def _compile_extended_topology(self, compiled: CompiledView, initial_state_id: int) -> None:
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
            pos, constraint_state_id = self.states[state_id]

            if pos == self.final_pos:
                self.child_results_by_state_id[state_id] = []
                continue

            constraint_states = self.constraint_states[constraint_state_id]
            new_constraint_states = constraint_states[-n_new_constraints:]
            child_results = []

            for choice, old_child_id in compiled.child_results_by_state_id[old_state_id]:
                next_new_states = self._advance_constraints(new_constraint_states, pos,  choice)

                if next_new_states is None:
                    continue

                child_pos, old_child_constraint_state_id = compiled.states[old_child_id]
                old_child_states = compiled.constraint_states[old_child_constraint_state_id]

                child_id, is_new = self._get_or_register_state_id(
                    child_pos,
                    old_child_states + next_new_states,
                )

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

        for state_id, (pos, _constraint_state_id) in enumerate(self.states):
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
        self.descendant_counts[state_id] = 1
        self.descendant_log_masses[state_id] = 0.0
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
        pos, _constraint_state_id = self.states[state_id]

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

            child_pos, _child_constraint_state_id = self.states[child_id]
            child_count = self.descendant_counts[child_id]
            subtree_log_mass = self.descendant_log_masses[child_id]

            if child_count is None or subtree_log_mass is None:
                continue

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
        self.descendant_counts[state_id] = descendant_count
        self.descendant_log_masses[state_id] = descendant_log_mass

    def _compiled_view(self, initial_state_id: int) -> CompiledView:
        """
        Build an immutable CompiledView from the compiler's current state.

        Parameters
        ----------
        initial_state_id
            The initial state ID.
        """
        initial_count = self.descendant_counts[initial_state_id]
        assert initial_count is not None

        choices_by_state_id = tuple(choices or {} for choices in self.choices_by_state_id)

        return CompiledView(
            initial_state=self.states[initial_state_id],
            initial_state_id=initial_state_id,
            states=tuple(self.states),
            constraint_states=tuple(self.constraint_states),
            child_results_by_state_id=tuple(tuple(results or ()) for results in self.child_results_by_state_id),
            choices_by_state_id=choices_by_state_id,
            choice_results_by_state_id=tuple(tuple(choices.values()) for choices in choices_by_state_id),
            n_valid_sequences=initial_count,
        )

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

        pos, constraint_state_id = self.states[state_id]
        constraint_states = self.constraint_states[constraint_state_id]

        for choice, child_pos in self.transitions_by_pos[pos]:
            next_constraint_states = self._advance_constraints(constraint_states, pos, choice)

            if next_constraint_states is None:
                continue

            child_id, is_new = self._get_or_register_state_id(child_pos, next_constraint_states)

            child_results.append((choice, child_id))

            if is_new:
                uncompiled_children.append((child_id, False))

        self.child_results_by_state_id[state_id] = child_results

        return uncompiled_children

    def _advance_constraints(
            self,
            constraint_states: ConstraintStates,
            pos: int,
            choice: str,
    ) -> Optional[ConstraintStates]:
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
