from codeine.constraints.base import PathConstraint


def test_base_path_constraint_does_nothing():
    constraint = PathConstraint()

    state = constraint.initial_state
    assert state == ()

    assert constraint.advance(state, 1, 'ATA') == state
    assert constraint.is_satisfied(state)
