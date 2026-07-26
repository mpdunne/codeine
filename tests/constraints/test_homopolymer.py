import pytest

from codeine.graph.base import CodonGraph
from codeine.constraints import HomopolymerConstraint


def test_homopolymer_constraint_builds_expected_banned_sequences():
    constraint = HomopolymerConstraint(max_length=6)

    assert constraint.max_length == 6

    assert set(constraint.banned_sequences) == {
        'AAAAAAA',
        'CCCCCCC',
        'GGGGGGG',
        'TTTTTTT',
    }


@pytest.mark.parametrize('max_length', [1, 2, 6, 20, 1000])
def test_homopolymer_constraint_accepts_positive_integer_max_length(max_length):
    constraint = HomopolymerConstraint(max_length=max_length)

    assert constraint.max_length == max_length


@pytest.mark.parametrize('max_length', [0, -1, -10])
def test_homopolymer_constraint_rejects_non_positive_max_length(max_length):
    with pytest.raises(ValueError, match='max_length must be at least 1'):
        HomopolymerConstraint(max_length=max_length)


@pytest.mark.parametrize('max_length', [1.5, '6', None])
def test_homopolymer_constraint_rejects_non_integer_max_length(max_length):
    with pytest.raises(TypeError, match='max_length must be an integer'):
        HomopolymerConstraint(max_length=max_length)


def test_homopolymer_constraint_filters_homopolymers():
    view = CodonGraph('KK').view()
    view.add_constraints([HomopolymerConstraint(5)])

    assert 'AAAAAG' in view
    assert 'AAAAAA' not in view
