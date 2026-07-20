from codeine.constraints.base import Constraint, DEAD_STATE, SAFE_STATE
from codeine.graph.base import CodonGraph
from codeine.graph.compile import ViewCompiler


class RejectChoicesConstraint(Constraint):
    """
    Test constraint that rejects specified graph choices.
    """

    def __init__(self, rejected_choices):
        self.rejected_choices = frozenset(rejected_choices)

    @property
    def initial_state(self):
        return ()

    @property
    def is_trivial(self):
        return not self.rejected_choices

    def link(self, graph):
        pass

    def advance(self, state, pos, choice):
        if state == SAFE_STATE:
            return SAFE_STATE

        if choice in self.rejected_choices:
            return DEAD_STATE

        return state


def test_existing_state_does_not_interrupt_dense_state_ids():
    compiler = ViewCompiler(CodonGraph('MIKEY').view())

    assert compiler._get_or_register_state_id(compiler.initial_pos, ()) == (0, True)
    assert compiler._get_or_register_state_id(1, ()) == (1, True)
    assert compiler._get_or_register_state_id(compiler.initial_pos, ()) == (0, False)
    assert compiler._get_or_register_state_id(2, ()) == (2, True)

    assert compiler.state_positions == [
        compiler.initial_pos,
        1,
        2,
    ]
    assert compiler.state_constraint_state_ids == [0, 0, 0]


def test_extend_matches_full_compile():
    graph = CodonGraph('MIKEY')
    constraint = RejectChoicesConstraint({'ATA'})

    base_compiled = ViewCompiler(graph.view()).compile()
    extended_compiled = ViewCompiler(graph.view()).extend(base_compiled, [constraint])
    full_compiled = ViewCompiler(graph.view(constraints=(constraint,))).compile()

    assert extended_compiled == full_compiled


def test_extend_with_multiple_constraints_matches_full_compile():
    graph = CodonGraph('MIKEY')
    base_constraint = RejectChoicesConstraint({'ATG'})
    new_constraints = (
        RejectChoicesConstraint({'AAG'}),
        RejectChoicesConstraint({'TAT'}),
    )

    base_compiled = ViewCompiler(
        graph.view(constraints=(base_constraint,))
    ).compile()

    extended_view = graph.view(
        constraints=(base_constraint, *new_constraints),
    )
    extended_compiled = ViewCompiler(extended_view).extend(
        base_compiled,
        new_constraints,
    )
    full_compiled = ViewCompiler(extended_view).compile()

    assert extended_compiled == full_compiled


def test_extend_from_constrained_compile_matches_full_compile():
    graph = CodonGraph('MIKEY')
    base_constraint = RejectChoicesConstraint({'ATA'})
    new_constraint = RejectChoicesConstraint({'AAG'})

    base_compiled = ViewCompiler(graph.view(constraints=[base_constraint])).compile()

    extended_compiled = ViewCompiler(graph.view()).extend(base_compiled, [new_constraint])
    full_compiled = ViewCompiler(graph.view(constraints=[base_constraint, new_constraint])).compile()

    assert extended_compiled == full_compiled


def test_extend_with_no_constraints_matches_shallow_compile():
    graph = CodonGraph('MIKEY')

    constraint = RejectChoicesConstraint({'ATA'})
    constrained_view = graph.view(constraints=[constraint])

    base_compiled = ViewCompiler(constrained_view).compile()

    extended_compiled = ViewCompiler(graph.view()).extend(base_compiled, [])
    shallow_compiled = ViewCompiler(graph.view()).compile_shallow(base_compiled)

    assert extended_compiled == shallow_compiled


def test_extend_ignores_trivial_constraints():
    view = CodonGraph('MIKEY').view()

    base_compiled = ViewCompiler(view).compile()
    extended_compiled = ViewCompiler(view).extend(base_compiled, [RejectChoicesConstraint(set())])
    shallow_compiled = ViewCompiler(view).compile_shallow(base_compiled)

    assert extended_compiled == shallow_compiled


def test_extend_does_not_add_sequences():
    graph = CodonGraph('MIKEY')

    base_compiled = ViewCompiler(graph.view()).compile()
    extended_compiled = ViewCompiler(graph.view()).extend(base_compiled, [RejectChoicesConstraint({'AAG'})])

    assert extended_compiled.n_valid_sequences <= base_compiled.n_valid_sequences


def test_extend_can_remove_all_sequences():
    graph = CodonGraph('MIKEY')
    unconstrained_view = graph.view()
    constraint = RejectChoicesConstraint({'ATG'})

    base_compiled = ViewCompiler(unconstrained_view).compile()
    extended_compiled = ViewCompiler(unconstrained_view).extend(base_compiled, [constraint])

    constrained_view = graph.view(constraints=[constraint])
    full_compiled = ViewCompiler(constrained_view).compile()

    assert extended_compiled == full_compiled
    assert extended_compiled.n_valid_sequences == 0


def test_compiled_state_ids_are_dense():
    view = CodonGraph('MIKEY').view(constraints=[RejectChoicesConstraint({'ATA'})])

    compiled = ViewCompiler(view).compile()

    referenced_state_ids = set(compiled.child_transition_state_ids)

    assert compiled.initial_state_id == 0
    assert referenced_state_ids == set(range(1, len(compiled.state_positions)))


def test_extended_state_ids_are_dense():
    view = CodonGraph('MIKEY').view()

    base_compiled = ViewCompiler(view).compile()
    extended_compiled = ViewCompiler(view).extend(base_compiled, [RejectChoicesConstraint({'ATA'})])

    referenced_state_ids = set(extended_compiled.child_transition_state_ids)

    assert extended_compiled.initial_state_id == 0
    assert referenced_state_ids == set(range(1, len(extended_compiled.state_positions)))
