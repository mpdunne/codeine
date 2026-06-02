import pytest

from collections import Counter

from codeine.utils.sampling import Sampler


def test_sampler_single_value_returns_item():
    s = Sampler([1])
    for _ in range(100):
        assert s.sample() == 1


def test_sampler_weights_work():
    s = Sampler([1, 2], weights=[5, 10], seed=8675309)

    sampled_items = []
    for _ in range(100000):
        sampled_items.append(s.sample())

    counts = Counter(sampled_items)
    assert counts[2] > counts[1] * 1.5


def test_sampler_multiple_values_return_all_items():
    items = [1, 2, 3, 4, 5]
    s = Sampler(items, seed=8675309)

    sampled_items = []
    for _ in range(100000):
        sampled_items.append(s.sample())

    assert set(sampled_items) == set(items)


@pytest.mark.parametrize("seed", [0, 5318008, 'hello', 0.401])
def test_sample_seed_consistent(seed):
    sampled_item_lists = []
    for rep in range(10):
        s = Sampler([1, 2, 3, 4, 5], weights=[5, 10, 1, 2, 13], seed=seed)
        sampled_items = []
        for _ in range(100):
            sampled_items.append(s.sample())
        sampled_item_lists.append(sampled_items)

    assert all(l == sampled_item_lists[0] for l in sampled_item_lists)


def test_sample_no_seed_is_not_consistent():
    sampled_item_lists = []
    for rep in range(10):
        s = Sampler([1, 2, 3, 4, 5], weights=[5, 10, 1, 2, 13])
        sampled_items = []
        for _ in range(1000):
            sampled_items.append(s.sample())
        sampled_item_lists.append(sampled_items)

    assert not all(l == sampled_item_lists[0] for l in sampled_item_lists)


def test_invalid_inputs_raise_errors():
    with pytest.raises(ValueError):
        Sampler([])

    with pytest.raises(ValueError):
        Sampler([1], weights=[-1])

    with pytest.raises(ValueError):
        Sampler([1], weights=[0])

    with pytest.raises(ValueError):
        Sampler([1, 2], weights=[1, 2, 4])

    with pytest.raises(ValueError):
        Sampler(['a', 'b'], weights=[1])

    with pytest.raises(ValueError):
        Sampler(['a', 'b'], weights=[0, 0])

    with pytest.raises(ValueError):
        Sampler(['a', 'b'], weights=[1, -2])
