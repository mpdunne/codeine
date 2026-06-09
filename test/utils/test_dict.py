import pytest

from codeine.utils.dict import FrozenDict


def test_getitem():
    fd = FrozenDict({'A': 420, 'B': 69})
    assert fd['A'] == 420
    assert fd['B'] == 69


def test_len():
    fd = FrozenDict({'A': 420, 'B': 69})
    assert len(fd) == 2


def test_iter():
    fd = FrozenDict({'A': 420, 'B': 69})
    assert list(fd) == ['A', 'B']


def test_contains():
    fd = FrozenDict({'A': 420, 'B': 69})
    assert 'A' in fd
    assert 'C' not in fd


def test_repr_matches_dict():
    data = {'A': 420, 'B': 69}
    fd = FrozenDict(data)
    assert repr(fd) == repr(data)


def test_items():
    fd = FrozenDict({'A': 420, 'B': 69})
    assert list(fd.items()) == [('A', 420), ('B', 69)]


def test_keys():
    fd = FrozenDict({'A': 420, 'B': 69})
    assert list(fd.keys()) == ['A', 'B']


def test_values():
    fd = FrozenDict({'A': 420, 'B': 69})
    assert list(fd.values()) == [420, 69]


def test_missing_key_raises():
    fd = FrozenDict({'A': 420})
    with pytest.raises(KeyError):
        _ = fd['B']


def test_data_doesnt_change_if_input_dict_does():
    data = {'A': 420}
    fd = FrozenDict(data)
    data['A'] = 69
    assert fd['A'] == 420


def test_item_assignment_not_allowed():
    fd = FrozenDict({'A': 420})
    with pytest.raises(TypeError):
        fd['A'] = 2


def test_item_deletion_not_allowed():
    fd = FrozenDict({'A': 420})
    with pytest.raises(TypeError):
        del fd['A']
