import pickle
import pytest
import random

from collections import Counter

from codeine.utils.sampling import SingletonSampler, WeightedSampler, UniformSampler


###################################
# Singleton sampler
###################################

def test_singleton_sampler_single_value_returns_item():
    s = SingletonSampler(1)
    for _ in range(100):
        assert s.sample() == 1


def test_singleton_sampler_pickle():
    sampler = SingletonSampler('a')
    loaded = pickle.loads(pickle.dumps(sampler))
    assert loaded.sample() == 'a'


###################################
# Random samplers
###################################

@pytest.mark.parametrize('sampler', (WeightedSampler, UniformSampler))
def test_random_sampler_empty_items_raise_error(sampler):
    with pytest.raises(ValueError):
        sampler([])


@pytest.mark.parametrize('sampler', (WeightedSampler, UniformSampler))
def test_random_sampler_seed_and_rng_cannot_both_be_provided(sampler):
    with pytest.raises(ValueError):
        sampler([1, 2, 3], seed=123, rng=random.Random(123))

    with pytest.raises(ValueError):
        sampler([55], seed=123, rng=random.Random(123))


@pytest.mark.parametrize('sampler', (WeightedSampler, UniformSampler))
def test_random_sampler_seed_strategies_exhibit_similar_behaviour(sampler):
    items = [1, 2, 3, 4, 5]

    s1 = sampler(items)
    s2 = sampler(items, seed=5318008)
    s3 = sampler(items, rng=random.Random(8675309))

    sampled1 = [s1.sample() for _ in range(10000)]
    sampled2 = [s2.sample() for _ in range(10000)]
    sampled3 = [s3.sample() for _ in range(10000)]

    assert set(sampled1) == set(items)
    assert set(sampled2) == set(items)
    assert set(sampled3) == set(items)

    counts1 = Counter(sampled1)
    counts2 = Counter(sampled2)
    counts3 = Counter(sampled3)

    assert all(500 < counts1[i] < 5000 for i in items)
    assert all(500 < counts2[i] < 5000 for i in items)
    assert all(500 < counts3[i] < 5000 for i in items)


@pytest.mark.parametrize('sampler', (WeightedSampler, UniformSampler))
def test_random_sampler_rng_consistent(sampler):
    rng1 = random.Random(8675309)
    s1 = sampler([1, 2, 3], rng=rng1)

    rng2 = random.Random(8675309)
    s2 = sampler([1, 2, 3], rng=rng2)

    for _ in range(100):
        assert s1.sample() == s2.sample()


@pytest.mark.parametrize('sampler', (WeightedSampler, UniformSampler))
def test_random_sampler_rng_object_is_advanced(sampler):
    rng = random.Random(8675309)
    s = sampler(['x', 'y'], rng=rng)

    before = rng.getstate()
    s.sample()
    after = rng.getstate()

    assert before != after


@pytest.mark.parametrize('sampler', (WeightedSampler, UniformSampler))
def test_weighted_sampler_items_are_copied(sampler):
    items = ['a', 'b']
    s = sampler(items, seed=8675309)
    items += ['c']

    assert s._items == ('a', 'b')


@pytest.mark.parametrize('sampler', (WeightedSampler, UniformSampler))
def test_weighted_sampler_multiple_values_return_all_items(sampler):
    items = [1, 2, 3, 4, 5]
    s = sampler(items, seed=8675309)

    sampled_items = []
    for _ in range(100000):
        sampled_items.append(s.sample())

    assert set(sampled_items) == set(items)


@pytest.mark.parametrize('sampler', (WeightedSampler, UniformSampler))
def test_weighted_sampler_pickle_preserves_random_state(sampler):
    s = sampler(['a', 'b', 'c'], seed=8675309)
    _ = [s.sample() for _ in range(100)]

    dumped = pickle.dumps(s)
    samples_orig = [s.sample() for _ in range(100)]

    loaded = pickle.loads(dumped)
    samples_loaded = [loaded.sample() for _ in range(100)]

    assert samples_orig == samples_loaded


###################################
# Uniform samplers
###################################


def test_uniform_sampler_samples_items_uniformly():
    sampler = UniformSampler(['a', 'b', 'c'], seed=8675309)
    counts = Counter(sampler.sample() for _ in range(30000))
    assert all(9000 < counts[item] < 11000 for item in ('a', 'b', 'c'))


@pytest.mark.parametrize('seed', [0, 5318008, 'hello', 0.401])
def test_uniform_sampler_sample_seed_consistent(seed):
    sampled_item_lists = []
    for rep in range(10):
        s = UniformSampler([1, 2, 3, 4, 5], seed=seed)
        sampled_items = []
        for _ in range(100):
            sampled_items.append(s.sample())
        sampled_item_lists.append(sampled_items)

    assert all(l == sampled_item_lists[0] for l in sampled_item_lists)


def test_uniform_sampler_no_seed_is_not_consistent():
    sampled_item_lists = []
    for rep in range(10):
        s = UniformSampler([1, 2, 3, 4, 5])
        sampled_items = []
        for _ in range(1000):
            sampled_items.append(s.sample())
        sampled_item_lists.append(sampled_items)

    assert not all(l == sampled_item_lists[0] for l in sampled_item_lists)


###################################
# Weighted samplers
###################################

def test_weighted_sampler_weights_work():
    s = WeightedSampler([1, 2], weights=[5, 10], seed=8675309)

    sampled_items = []
    for _ in range(100000):
        sampled_items.append(s.sample())

    counts = Counter(sampled_items)
    assert counts[2] > counts[1] * 1.5


def test_weighted_sampler_weights_are_copied():
    weights = [1, 100]
    sampler = WeightedSampler(['a', 'b'], weights=weights, seed=8675309)
    weights[0] = 5318008

    sampled = [sampler.sample() for _ in range(1000)]
    assert sampled.count('b') > sampled.count('a')


@pytest.mark.parametrize('seed', [0, 5318008, 'hello', 0.401])
def test_weighted_sampler_sample_seed_consistent(seed):
    sampled_item_lists = []
    for rep in range(10):
        s = WeightedSampler([1, 2, 3, 4, 5], weights=[5, 10, 1, 2, 13], seed=seed)
        sampled_items = []
        for _ in range(100):
            sampled_items.append(s.sample())
        sampled_item_lists.append(sampled_items)

    assert all(l == sampled_item_lists[0] for l in sampled_item_lists)


def test_weighted_sampler_no_seed_is_not_consistent():
    sampled_item_lists = []
    for rep in range(10):
        s = WeightedSampler([1, 2, 3, 4, 5], weights=[5, 10, 1, 2, 13])
        sampled_items = []
        for _ in range(1000):
            sampled_items.append(s.sample())
        sampled_item_lists.append(sampled_items)

    assert not all(l == sampled_item_lists[0] for l in sampled_item_lists)


def test_weighted_sampler_weights_can_be_ints_or_floats():
    s = WeightedSampler(['a', 'b', 'c'], weights=[1, 0.5, 2.5], seed=123)
    sampled_items = [s.sample() for _ in range(1000)]
    assert set(sampled_items) == {'a', 'b', 'c'}


def test_weighted_sampler_zero_weight_items_are_never_sampled():
    s = WeightedSampler(['a', 'b', 'c'], weights=[1, 0, 1], seed=123)
    sampled_items = [s.sample() for _ in range(1000)]
    assert 'b' not in sampled_items
    assert set(sampled_items) == {'a', 'c'}


def test_weighted_sampler_default_weights_are_uniform():
    s = WeightedSampler(['a', 'b'], seed=8675309)
    sampled_items = [s.sample() for _ in range(10000)]
    counts = Counter(sampled_items)
    assert 0.8 < counts['a'] / counts['b'] < 1.25


def test_weighted_sampler_invalid_inputs_raise_errors():

    with pytest.raises(ValueError):
        WeightedSampler([1], weights=[-1])

    with pytest.raises(ValueError):
        WeightedSampler([1], weights=[0])

    with pytest.raises(ValueError):
        WeightedSampler([1, 2], weights=[1, 2, 4])

    with pytest.raises(ValueError):
        WeightedSampler(['a', 'b'], weights=[1])

    with pytest.raises(ValueError):
        WeightedSampler(['a', 'b'], weights=[0, 0])

    with pytest.raises(ValueError):
        WeightedSampler(['a', 'b'], weights=[1, -2])


def test_weighted_sampler_pickle_preserves_config():
    sampler = WeightedSampler(['a', 'b', 'c'], weights=[1, 2, 3], seed=123)
    loaded = pickle.loads(pickle.dumps(sampler))

    assert loaded._items == ('a', 'b', 'c')
    assert loaded._cumulative == sampler._cumulative
