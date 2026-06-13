import pickle
import pytest
import random

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


@pytest.mark.parametrize('seed', [0, 5318008, 'hello', 0.401])
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


def test_sampler_seed_and_rng_cannot_both_be_provided():
    import random
    with pytest.raises(ValueError):
        Sampler([1, 2, 3], seed=123, rng=random.Random(123))


def test_sampler_seed_strategies_exhibit_similar_behaviour():
    items = [1, 2, 3, 4, 5]

    s1 = Sampler(items)
    s2 = Sampler(items, seed=5318008)
    s3 = Sampler(items, rng=random.Random(8675309))

    sampled1 = [s1.sample() for _ in range(10000)]
    sampled2 = [s2.sample() for _ in range(10000)]
    sampled3 = [s3.sample() for _ in range(10000)]

    assert set(sampled1) == set(items)
    assert set(sampled2) == set(items)
    assert set(sampled3) == set(items)

    counts1 = Counter(sampled1)
    counts2 = Counter(sampled1)
    counts3 = Counter(sampled1)

    assert [500 < counts1[i] < 5000 for i in items]
    assert [500 < counts2[i] < 5000 for i in items]
    assert [500 < counts3[i] < 5000 for i in items]


def test_sampler_rng_consistent():
    rng1 = random.Random(8675309)
    s1 = Sampler([1, 2, 3], rng=rng1)

    rng2 = random.Random(8675309)
    s2 = Sampler([1, 2, 3], rng=rng2)

    for _ in range(100):
        assert s1.sample() == s2.sample()


def test_weights_can_be_ints_or_floats():
    s = Sampler(['a', 'b', 'c'], weights=[1, 0.5, 2.5], seed=123)
    sampled_items = [s.sample() for _ in range(1000)]
    assert set(sampled_items) == {'a', 'b', 'c'}


def test_zero_weight_items_are_never_sampled():
    s = Sampler(['a', 'b', 'c'], weights=[1, 0, 1], seed=123)
    sampled_items = [s.sample() for _ in range(1000)]
    assert 'b' not in sampled_items
    assert set(sampled_items) == {'a', 'c'}


def test_default_weights_are_uniform():
    s = Sampler(['a', 'b'], seed=8675309)
    sampled_items = [s.sample() for _ in range(10000)]
    counts = Counter(sampled_items)
    assert 0.8 < counts['a'] / counts['b'] < 1.25


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


def test_sampler_pickle_preserves_random_state():
    sampler = Sampler(['a', 'b', 'c'], weights=[1, 2, 3], seed=8675309)
    _ = [sampler.sample() for _ in range(100)]

    dumped = pickle.dumps(sampler)
    samples_orig = [sampler.sample() for _ in range(100)]

    loaded = pickle.loads(dumped)
    samples_loaded = [loaded.sample() for _ in range(100)]

    assert samples_orig == samples_loaded


def test_single_value_sampler_pickle():
    sampler = Sampler(['a'])
    loaded = pickle.loads(pickle.dumps(sampler))
    assert loaded.sample() == 'a'
    assert loaded.items == ('a',)