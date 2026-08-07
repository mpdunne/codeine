from codeine.utils.tuples import tuplify


def test_tuplify_single_item():
    assert tuplify(3, int) == (3,)


def test_tuplify_iterable():
    assert tuplify([1, 2, 3], int) == (1, 2, 3)


def test_tuplify_tuple():
    assert tuplify((1, 2, 3), int) == (1, 2, 3)


def test_tuplify_none():
    assert tuplify(None, int) == ()


def test_tuplify_multiple_item_types():
    assert tuplify(3, (int, str)) == (3,)
    assert tuplify('hello', (int, str)) == ('hello',)


def test_tuplify_iterable_with_multiple_item_types():
    assert tuplify([1, 'hello'], (int, str)) == (1, 'hello')