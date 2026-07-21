import pytest

from itertools import product

from codeine.constraints.base import DEAD_STATE
from codeine.constraints._count import _CountConstraint, INITIAL_STATE, SAFE_STATE


class BinaryCountConstraint(_CountConstraint):
    """
    Count choices equal to 'B'.

    Each choice contributes one countable position.
    """

    _positions_per_choice = 1
    _choice_counts = {'A': 0, 'B': 1}


@pytest.mark.parametrize(
    ('name', 'value'),
    [
        ('min_frac', -0.01),
        ('max_frac', 1.01),
        ('min_perc', -1),
        ('max_perc', 101),
        ('min_frac', float('inf')),
        ('max_perc', float('nan')),
    ],
)
def test_invalid_relative_bounds(name, value):
    with pytest.raises(ValueError):
        BinaryCountConstraint(**{name: value})


@pytest.mark.parametrize(
    ('name', 'value'),
    [
        ('min_frac', '0.5'),
        ('max_perc', object()),
        ('min_frac', True),
    ],
)
def test_non_numeric_relative_bounds(name, value):
    with pytest.raises(TypeError):
        BinaryCountConstraint(**{name: value})


@pytest.mark.parametrize(
    ('name', 'value'),
    [
        ('min_count', -1),
        ('max_count', -1),
    ],
)
def test_negative_count_bounds(name, value):
    with pytest.raises(ValueError):
        BinaryCountConstraint(**{name: value})


@pytest.mark.parametrize(
    ('name', 'value'),
    [
        ('min_count', 1.5),
        ('max_count', '2'),
        ('min_count', True),
    ],
)
def test_non_integer_count_bounds(name, value):
    with pytest.raises(TypeError):
        BinaryCountConstraint(**{name: value})


class MockCodonNode:
    def __init__(self, codons):
        self.codons = codons


class MockCodonGraph:
    def __init__(self, codons_by_position):
        self.codon_nodes = tuple(MockCodonNode(codons) for codons in codons_by_position)

    def codon_node_by_pos(self, pos):
        return self.codon_nodes[pos - 1]


def test_no_bounds_is_safe():
    constraint = BinaryCountConstraint()
    constraint.link(MockCodonGraph([['A', 'B']] * 4))

    assert constraint.initial_state == SAFE_STATE


def test_fraction_bounds_resolve_to_counts():
    constraint = BinaryCountConstraint(min_frac=0.26, max_frac=0.74)
    constraint.link(MockCodonGraph([['A', 'B']] * 4))

    assert constraint._resolved_min_count == 2
    assert constraint._resolved_max_count == 2


def test_percentage_bounds_resolve_to_counts():
    constraint = BinaryCountConstraint(min_perc=25, max_perc=75,)
    constraint.link(MockCodonGraph([['A', 'B']] * 4))

    assert constraint._resolved_min_count == 1
    assert constraint._resolved_max_count == 3


def test_multiple_minimum_bounds_use_strictest():
    constraint = BinaryCountConstraint(
        min_frac=0.25,
        min_perc=50,
        min_count=3,
    )
    constraint.link(MockCodonGraph([['A', 'B']] * 4))

    assert constraint._resolved_min_count == 3


def test_multiple_maximum_bounds_use_strictest():
    constraint = BinaryCountConstraint(
        max_frac=0.75,
        max_perc=50,
        max_count=3,
    )
    constraint.link(MockCodonGraph([['A', 'B']] * 4))

    assert constraint._resolved_max_count == 2


def test_incompatible_bounds_kill_constraint():
    constraint = BinaryCountConstraint(
        min_count=3,
        max_count=2,
    )
    constraint.link(MockCodonGraph([['A', 'B']] * 4))

    assert constraint.initial_state == DEAD_STATE


def test_initial_state_dead_when_minimum_is_unreachable():
    graph = MockCodonGraph([
        ['A'],
        ['A', 'B'],
        ['A'],
    ])

    constraint = BinaryCountConstraint(min_count=2)
    constraint.link(graph)

    assert constraint.initial_state == DEAD_STATE


def test_initial_state_dead_when_maximum_is_unavoidable():
    graph = MockCodonGraph([
        ['B'],
        ['A', 'B'],
        ['B'],
    ])

    constraint = BinaryCountConstraint(max_count=1)
    constraint.link(graph)

    assert constraint.initial_state == DEAD_STATE


def test_initial_state_safe_when_every_sequence_is_valid():
    graph = MockCodonGraph([
        ['B'],
        ['A', 'B'],
        ['A'],
    ])

    constraint = BinaryCountConstraint(min_count=1, max_count=2)
    constraint.link(graph)

    assert constraint.initial_state == SAFE_STATE


def test_initial_state_tracks_when_some_sequences_are_valid():
    graph = MockCodonGraph([
        ['A', 'B'],
        ['A', 'B'],
    ])

    constraint = BinaryCountConstraint(min_count=1, max_count=1)
    constraint.link(graph)

    assert constraint.initial_state == INITIAL_STATE


def test_dead_state_is_absorbing():
    constraint = BinaryCountConstraint(min_count=1)
    constraint.link(MockCodonGraph([['A', 'B']]))

    assert constraint.advance(DEAD_STATE, 1, 'B') == DEAD_STATE


def test_safe_state_is_absorbing():
    constraint = BinaryCountConstraint()
    constraint.link(MockCodonGraph([['A', 'B']]))

    assert constraint.advance(SAFE_STATE, 1, 'B') == SAFE_STATE


def test_advance_returns_dead_when_minimum_becomes_unreachable():
    graph = MockCodonGraph([
        ['A', 'B'],
        ['A'],
        ['A', 'B'],
    ])

    constraint = BinaryCountConstraint(min_count=2)
    constraint.link(graph)

    # Choosing A at position zero leaves only one possible B later.
    assert constraint.advance(INITIAL_STATE, 1, 'A') == DEAD_STATE


def test_advance_returns_dead_when_maximum_becomes_unavoidable():
    graph = MockCodonGraph([
        ['A', 'B'],
        ['B'],
        ['B'],
    ])

    constraint = BinaryCountConstraint(max_count=1)
    constraint.link(graph)

    # Even choosing A now leaves two unavoidable Bs.
    assert constraint.advance(INITIAL_STATE, 1, 'A') == DEAD_STATE


def test_advance_returns_safe_when_all_completions_are_valid():
    graph = MockCodonGraph([
        ['A', 'B'],
        ['A', 'B'],
        ['A'],
    ])

    constraint = BinaryCountConstraint(min_count=1,max_count=2)
    constraint.link(graph)

    state = constraint.advance(INITIAL_STATE, 1, 'B')

    assert state == SAFE_STATE


def test_advance_keeps_count_when_future_outcome_still_matters():
    graph = MockCodonGraph([
        ['A', 'B'],
        ['A', 'B'],
    ])

    constraint = BinaryCountConstraint(min_count=1, max_count=1)
    constraint.link(graph)

    assert constraint.advance(INITIAL_STATE, 1, 'A') == 0
    assert constraint.advance(INITIAL_STATE, 1, 'B') == 1


def helper_walk_constraint(constraint, sequence):
    state = constraint.initial_state

    for pos, choice in enumerate(sequence, 1):
        state = constraint.advance(state, pos, choice)

    return state


@pytest.mark.parametrize(
    'kwargs',
    (
        {},
        {'min_count': 2},
        {'max_count': 1},
        {'min_count': 1, 'max_count': 2},
        {'min_frac': 0.5},
        {'max_perc': 50},
        {'min_frac': 0.25, 'min_count': 2},
        {'max_perc': 75, 'max_count': 2},
        {'min_count': 3, 'max_count': 1},
    ),
)
@pytest.mark.parametrize(
    'graph',
    (
        MockCodonGraph([]),
        MockCodonGraph([
            ['A'],
            ['A'],
            ['A'],
            ['A'],
        ]),
        MockCodonGraph([
            ['B'],
            ['B'],
            ['B'],
            ['B'],
        ]),
        MockCodonGraph([
            ['A'],
            ['A'],
            ['A', 'B'],
            ['A', 'B'],
        ]),
        MockCodonGraph([
            ['A', 'B'],
            ['A', 'B'],
            ['A', 'B'],
            ['A', 'B'],
            ['A', 'B'],
            ['A', 'B'],
        ]),
        
    )
)
def test_constraint_matches_enumerate_reject(graph, kwargs):

    constraint = BinaryCountConstraint(**kwargs)
    constraint.link(graph)

    codons = (node.codons for node in graph.codon_nodes)
    sequences = product(*codons)

    for sequence in sequences:
        count = sum(choice == 'B' for choice in sequence)

        minimum = constraint._resolved_min_count
        maximum = constraint._resolved_max_count
        should_reject = not (minimum <= count <= maximum)

        state = helper_walk_constraint(constraint, sequence)
        did_reject = (state == DEAD_STATE)

        assert did_reject == should_reject, sequence
