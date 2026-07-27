import pytest

from codeine.constraints.tandem import TandemRepeatConstraint


def test_tandem_repeat_constraint_stores_parameters():
    constraint = TandemRepeatConstraint(repeat_length=4, min_copies=3)

    assert constraint.repeat_length == 4
    assert constraint.min_copies == 3


def test_tandem_repeat_constraint_defaults_to_two_copies():
    constraint = TandemRepeatConstraint(repeat_length=4)

    assert constraint.min_copies == 2


@pytest.mark.parametrize('repeat_length', [0, -1])
def test_repeat_length_must_be_positive(repeat_length):
    with pytest.raises(ValueError, match='repeat_length must be at least 1'):
        TandemRepeatConstraint(repeat_length=repeat_length)


@pytest.mark.parametrize('repeat_length', [None, 1.5, '4', True,],)
def test_repeat_length_must_be_an_integer(repeat_length):
    with pytest.raises(TypeError, match='repeat_length must be an integer'):
        TandemRepeatConstraint(repeat_length=repeat_length)


@pytest.mark.parametrize('min_copies', [0, 1, -1])
def test_min_copies_must_be_at_least_two(min_copies):
    with pytest.raises(ValueError, match='min_copies must be at least 2'):
        TandemRepeatConstraint(repeat_length=4, min_copies=min_copies)


@pytest.mark.parametrize('min_copies', [None, 2.5, '2', True])
def test_min_copies_must_be_an_integer(min_copies):
    with pytest.raises(TypeError, match='min_copies must be an integer'):
        TandemRepeatConstraint(repeat_length=4, min_copies=min_copies)