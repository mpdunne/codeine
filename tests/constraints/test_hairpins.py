import pytest

from codeine.constraints.hairpins import Hairpins


def test_hairpin_constraint_stores_values_correctly():
    constraint = Hairpins(stem_length=6, min_loop_length=4, max_loop_length=20)

    assert constraint.stem_length == 6
    assert constraint.min_loop_length == 4
    assert constraint.max_loop_length == 20

    assert constraint.repeat_length == 6
    assert constraint.min_distance == 4
    assert constraint.max_distance == 20
    assert constraint.inverted is True


def test_hairpin_constraint_uses_defaults():
    constraint = Hairpins(stem_length=6)

    assert constraint.stem_length == 6
    assert constraint.min_loop_length == 3
    assert constraint.max_loop_length is None

    assert constraint.repeat_length == 6
    assert constraint.min_distance == 3
    assert constraint.max_distance is None
    assert constraint.inverted is True


@pytest.mark.parametrize('stem_length', [None, 1.0, '1', [], True, False])
def test_stem_length_must_be_integer(stem_length):
    with pytest.raises(TypeError, match='repeat_length must be an integer'):
        Hairpins(stem_length)


@pytest.mark.parametrize('stem_length', [0, -1, -100])
def test_stem_length_must_be_positive(stem_length):
    with pytest.raises(ValueError, match='repeat_length must be at least 1'):
        Hairpins(stem_length)


@pytest.mark.parametrize('min_loop_length', [None, 1.0, '1', [], True, False])
def test_min_loop_length_must_be_integer(min_loop_length):
    with pytest.raises(TypeError, match='min_distance must be an integer'):
        Hairpins(stem_length=6, min_loop_length=min_loop_length)


@pytest.mark.parametrize('min_loop_length', [-1, -100])
def test_min_loop_length_must_be_non_negative(min_loop_length):
    with pytest.raises(ValueError, match='min_distance must be at least 0'):
        Hairpins(stem_length=6, min_loop_length=min_loop_length)


@pytest.mark.parametrize('max_loop_length', [1.0, '10', [], True, False])
def test_max_loop_length_must_be_integer_or_none(max_loop_length):
    with pytest.raises(TypeError, match='max_distance must be an integer or None'):
        Hairpins(stem_length=6, max_loop_length=max_loop_length)


def test_max_loop_length_must_not_be_less_than_minimum():
    with pytest.raises(ValueError, match='max_distance must be at least min_distance'):
        Hairpins(stem_length=6, min_loop_length=5, max_loop_length=4)