import pytest

from codeine.graph.base import CodonGraph
from codeine.constraints import MaxHomopolymer


def test_homopolymer_constraint_builds_expected_banned_sequences():
    constraint = MaxHomopolymer(
        max_length=6,
        allow_single_interruption=False,
    )

    assert constraint.max_length == 6
    assert constraint.allow_single_interruption is False

    assert set(constraint.forbidden_sequences) == {
        'AAAAAAA',
        'CCCCCCC',
        'GGGGGGG',
        'TTTTTTT',
    }


def test_homopolymer_constraint_allows_single_interruption_by_default():
    constraint = MaxHomopolymer(max_length=5)

    assert constraint.allow_single_interruption is True
    assert 'CCCACCC' in constraint.forbidden_sequences


@pytest.mark.parametrize('max_length', [1, 2, 6, 20, 1000])
def test_homopolymer_constraint_accepts_positive_integer_max_length(max_length):
    constraint = MaxHomopolymer(
        max_length=max_length,
        allow_single_interruption=False,
    )

    assert constraint.max_length == max_length


@pytest.mark.parametrize('max_length', [0, -1, -10])
def test_homopolymer_constraint_rejects_non_positive_max_length(max_length):
    with pytest.raises(ValueError, match='max_length must be at least 1'):
        MaxHomopolymer(max_length=max_length)


@pytest.mark.parametrize('max_length', [1.5, '6', None])
def test_homopolymer_constraint_rejects_non_integer_max_length(max_length):
    with pytest.raises(TypeError, match='max_length must be an integer'):
        MaxHomopolymer(max_length=max_length)


@pytest.mark.parametrize('allow_single_interruption', [True, False])
def test_homopolymer_constraint_accepts_boolean_allow_single_interruption(
    allow_single_interruption,
):
    constraint = MaxHomopolymer(
        max_length=5,
        allow_single_interruption=allow_single_interruption,
    )

    assert constraint.allow_single_interruption is allow_single_interruption


@pytest.mark.parametrize('allow_single_interruption', [0, 1, 'true', None])
def test_homopolymer_constraint_rejects_non_boolean_allow_single_interruption(
    allow_single_interruption,
):
    with pytest.raises(
        TypeError,
        match='allow_single_interruption must be a boolean',
    ):
        MaxHomopolymer(
            max_length=5,
            allow_single_interruption=allow_single_interruption,
        )


def test_homopolymer_constraint_filters_homopolymers():
    view = CodonGraph('KK').view()
    view.add_constraints([MaxHomopolymer(5)])

    assert 'AAAAAG' in view
    assert 'AAAAAA' not in view


def test_homopolymer_constraint_filters_single_interruption():
    view = CodonGraph('PTQ').view()
    view.add_constraints([MaxHomopolymer(5)])

    assert 'CCCACCCAA' not in view


def test_homopolymer_constraint_allows_max_length_with_single_interruption():
    view = CodonGraph('PT').view()
    view.add_constraints([MaxHomopolymer(5)])

    assert 'CCCACC' in view


def test_homopolymer_constraint_does_not_filter_single_interruption_when_disabled():
    view = CodonGraph('PTQ').view()
    view.add_constraints([
        MaxHomopolymer(5, allow_single_interruption=False),
    ])

    assert 'CCCACCCAA' in view


def test_homopolymer_constraint_does_not_bridge_two_base_interruption():
    view = CodonGraph('PNP').view()
    view.add_constraints([MaxHomopolymer(5)])

    assert 'CCCAACCCA' in view
