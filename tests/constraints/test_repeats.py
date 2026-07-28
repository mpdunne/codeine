import pytest

from codeine.constraints.repeats import DirectRepeatConstraint, InvertedRepeatConstraint, RepeatConstraint


class ExampleRepeatConstraint(RepeatConstraint):
    pass


def test_repeat_constraint_stores_args_correctly():
    constraint = ExampleRepeatConstraint(repeat_length=10, min_distance=3, max_distance=20, inverted=True)

    assert constraint.repeat_length == 10
    assert constraint.min_distance == 3
    assert constraint.max_distance == 20
    assert constraint.inverted is True


def test_repeat_constraint_uses_defaults():
    constraint = ExampleRepeatConstraint(repeat_length=10)

    assert constraint.repeat_length == 10
    assert constraint.min_distance == 0
    assert constraint.max_distance is None
    assert constraint.inverted is False


@pytest.mark.parametrize('repeat_length', [None, 1.0, '1', [], True, False])
def test_repeat_length_must_be_integer(repeat_length):
    with pytest.raises(TypeError, match='repeat_length must be an integer'):
        ExampleRepeatConstraint(repeat_length)


@pytest.mark.parametrize('repeat_length', [0, -1, -100])
def test_repeat_length_must_be_positive(repeat_length):
    with pytest.raises(ValueError, match='repeat_length must be at least 1'):
        ExampleRepeatConstraint(repeat_length)


@pytest.mark.parametrize('min_distance', [None, 1.0, '1', [], True, False])
def test_min_distance_must_be_integer(min_distance):
    with pytest.raises(TypeError, match='min_distance must be an integer'):
        ExampleRepeatConstraint(repeat_length=10, min_distance=min_distance)


@pytest.mark.parametrize('min_distance', [-1, -100])
def test_min_distance_must_be_non_negative(min_distance):
    with pytest.raises(ValueError, match='min_distance must be at least 0'):
        ExampleRepeatConstraint(repeat_length=10, min_distance=min_distance)


@pytest.mark.parametrize('max_distance', [1.0, '10', [], True, False])
def test_max_distance_must_be_integer_or_none(max_distance):
    with pytest.raises(TypeError, match='max_distance must be an integer or None'):
        ExampleRepeatConstraint(repeat_length=10, max_distance=max_distance)


def test_max_distance_must_not_be_less_than_minimum():
    with pytest.raises(ValueError, match='max_distance must be at least min_distance'):
        ExampleRepeatConstraint(repeat_length=10, min_distance=5, max_distance=4)


@pytest.mark.parametrize('inverted', [None, 0, 1, 'yes'])
def test_inverted_must_be_boolean(inverted):
    with pytest.raises(TypeError, match='inverted must be a boolean'):
        ExampleRepeatConstraint(repeat_length=10, inverted=inverted)


def test_direct_repeat_stores_values_correctly():
    constraint = DirectRepeatConstraint(repeat_length=10, min_distance=3, max_distance=20)

    assert constraint.repeat_length == 10
    assert constraint.min_distance == 3
    assert constraint.max_distance == 20
    assert constraint.inverted is False


def test_inverted_repeat_stores_values_correctly():
    constraint = InvertedRepeatConstraint(repeat_length=10, min_distance=3, max_distance=20)

    assert constraint.repeat_length == 10
    assert constraint.min_distance == 3
    assert constraint.max_distance == 20
    assert constraint.inverted is True
