import pytest

from codeine.constraints.base import Constraint
from codeine.graph.factors import ChoiceFactor


@pytest.mark.parametrize('assignments,rejected', [
    ({1: 'AAA', 2: 'AAA', 3: 'CCC'}, True),
    ({1: 'AAA', 2: 'AAA', 3: 'CCA'}, False),
    ({1: 'AAG', 2: 'AAA', 3: 'CCC'}, False),
])
def test_choice_factor_requires_all_relations_to_match(assignments, rejected):
    factor = ChoiceFactor([
        (1, 2, {'AAA': frozenset({'AAA'})}),
        (2, 3, {'AAA': frozenset({'CCC'})}),
    ])

    assert factor.scope == frozenset({1, 2, 3})
    assert factor.rejects(assignments) is rejected


def test_choice_factor_can_compare_a_position_with_itself():
    factor = ChoiceFactor([(1, 1, {'AAA': frozenset({'AAA'})})])

    assert factor.rejects({1: 'AAA'})
    assert not factor.rejects({1: 'AAG'})


def test_choice_factor_without_relations_always_rejects():
    factor = ChoiceFactor([])

    assert factor.scope == frozenset()
    assert factor.rejects({})


def test_constraints_without_factor_support_return_none():
    class UnsupportedConstraint(Constraint):
        initial_state = 0

        def link(self, graph):
            pass

        def advance(self, state, pos, choice):
            return state

    constraint = UnsupportedConstraint()

    assert constraint.factors() is None
