import math
import random

from collections import Counter
from itertools import product
from unittest.mock import Mock, patch

import pytest

from codeine.constraints.repeats import DirectRepeats, InvertedRepeats
from codeine.graph.base import CodonGraph
from codeine.graph.factor_count import ComponentModelCounter
from codeine.graph.factor_sampling import FactorSampler
from codeine.graph.factors import ChoiceFactor
from codeine.tools.repeats import contains_direct_repeat


def make_counter():
    return ComponentModelCounter(
        {pos: ('AAA', 'AAG') for pos in range(5)},
        [
            ChoiceFactor([
                (0, 0, {'AAA': frozenset({'AAA'})}),
                (1, 2, {'AAA': frozenset({'AAA'})}),
            ]),
            ChoiceFactor([
                (0, 0, {'AAA': frozenset({'AAA'})}),
                (3, 4, {'AAA': frozenset({'AAA'})}),
            ]),
        ],
    )


@pytest.mark.parametrize('weighted', [False, True])
def test_sampling_matches_exhaustive_probability_distribution(weighted):
    counter = make_counter()
    weights = {pos: {'AAA': pos + 1, 'AAG': 2} for pos in counter.domains} if weighted else {}
    sampler = FactorSampler(counter, weights)
    masses = {}

    for choices in product(*counter.domains.values()):
        if counter.contains(dict(zip(counter.domains, choices))):
            masses[choices] = math.prod(weights.get(pos, {}).get(choice, 1) for pos, choice in enumerate(choices))

    total = sum(masses.values())
    plan = sampler.prepare_sampling()
    log_mass = sampler._mass_component(counter._initial_component_id)

    assert math.exp(log_mass) == pytest.approx(total)

    samples = sampler.sample(random.Random(42), 30000)
    observed = Counter(tuple(assignment[pos] for pos in counter.domains) for assignment in samples)

    assert set(observed) == set(masses)
    assert sampler.prepare_sampling() is plan

    for choices, mass in masses.items():
        probability = mass / total
        tolerance = 6 * math.sqrt(probability * (1 - probability) / len(samples)) + 1 / len(samples)

        assert abs(observed[choices] / len(samples) - probability) < tolerance


def test_sampling_reuses_plan_without_recomputing_components():
    counter = make_counter()
    sampler = FactorSampler(counter)
    sampler.prepare_sampling()

    with patch.object(counter, '_split_components', side_effect=AssertionError('Recomputed components')):
        samples = sampler.sample(random.Random(1), 20)

    assert all(counter.contains(assignment) for assignment in samples)


def test_shared_plan_uses_callers_random_generator():
    sampler = FactorSampler(make_counter())
    first_rng = random.Random(1)
    second_rng = random.Random(1)

    first = sampler.sample(first_rng, 20)
    sampler.sample(random.Random(9), 10)
    second = [sampler.sample(second_rng) for _ in range(20)]

    assert first == second
    assert sampler.sample(first_rng, 0) == []

    with pytest.raises(ValueError, match='non-negative'):
        sampler.sample(first_rng, -1)


def test_different_weights_share_counts_without_changing_each_other():
    counter = ComponentModelCounter({1: ('AAA', 'AAG')}, ())
    weights = {1: {'AAA': 1, 'AAG': 0}}
    first = FactorSampler(counter, weights)
    second = FactorSampler(counter, {1: {'AAA': 0, 'AAG': 1}})
    weights[1]['AAA'] = 0
    rng = Mock()

    assert first.sample(rng) == {1: 'AAA'}
    assert second.sample(rng) == {1: 'AAG'}
    assert counter.count() == 2
    assert counter.contains({1: 'AAG'})
    rng.random.assert_not_called()


@pytest.mark.parametrize('weight', [1e-300, 1e300])
def test_sampling_large_spaces_avoids_overflow_and_underflow(weight):
    counter = ComponentModelCounter({pos: ('AAA', 'AAG') for pos in range(1100)}, ())
    weights = {pos: {'AAA': weight, 'AAG': weight} for pos in counter.domains}
    sampler = FactorSampler(counter, weights)

    assignment = sampler.sample(random.Random(42))

    assert counter.contains(assignment)
    assert set(assignment.values()) == {'AAA', 'AAG'}


@pytest.mark.parametrize('weight', [-1, math.inf, -math.inf, math.nan])
def test_sampling_rejects_invalid_weights(weight):
    counter = ComponentModelCounter({1: ('AAA',)}, ())

    with pytest.raises(ValueError, match='finite and non-negative'):
        FactorSampler(counter, {1: {'AAA': weight}})


@pytest.mark.parametrize('factors,weights,error', [
    ([ChoiceFactor([])], {}, 'empty'),
    ([], {1: {'AAA': 0, 'AAG': 0}}, 'zero total weight'),
    ([ChoiceFactor([(1, 1, {'AAA': frozenset({'AAA'})})])],
     {1: {'AAA': 1, 'AAG': 0}}, 'zero total weight'),
])
def test_sampling_rejects_models_without_positive_mass(factors, weights, error):
    sampler = FactorSampler(ComponentModelCounter({1: ('AAA', 'AAG')}, factors), weights)

    with pytest.raises(ValueError, match=error):
        sampler.sample(random.Random(1))


def test_sampling_zero_random_draw_never_selects_zero_weight_branch():
    counter = make_counter()
    weights = {0: {'AAA': 0, 'AAG': 1}}
    sampler = FactorSampler(counter, weights)
    rng = Mock()
    rng.random.return_value = 0.0

    assignment = sampler.sample(rng)

    assert assignment[0] == 'AAG'
    assert counter.contains(assignment)


@pytest.mark.parametrize('constraint_type', [DirectRepeats, InvertedRepeats])
def test_sampling_repeat_factors_matches_flat_membership(constraint_type):
    graph = CodonGraph('LRFK', context_l='AT', context_r='GC')
    constraint = constraint_type(3)
    constraint.link(graph)
    domains = {node.pos: tuple(node.transitions) for node in graph.nodes if node is not graph.end_node}
    counter = ComponentModelCounter(domains, constraint.factors())
    view = graph.view(constraints=[constraint_type(3)])
    samples = FactorSampler(counter).sample(random.Random(1), 100)

    for assignment in samples:
        sequence = ''.join(assignment[pos] for pos in range(1, len(graph.aa_seq) + 1))

        assert counter.contains(assignment)
        assert sequence in view
        assert assignment[0] == 'AT'
        assert assignment[len(graph.aa_seq) + 1] == 'GC'


def test_sampling_hard_repeat_case():
    graph = CodonGraph('MIKEYSASSAFRASMIKEYSASSAFRAS')
    constraint = DirectRepeats(15)
    constraint.link(graph)
    domains = {node.pos: tuple(node.transitions) for node in graph.nodes if node is not graph.end_node}
    counter = ComponentModelCounter(domains, constraint.factors())
    samples = FactorSampler(counter).sample(random.Random(42), 100)

    for assignment in samples:
        sequence = ''.join(assignment[pos] for pos in sorted(assignment))

        assert counter.contains(assignment)
        assert not contains_direct_repeat(sequence, 15)


@pytest.mark.parametrize('draws,expected', [
    ([0.54], {0: 'AAA', 1: 'AAG'}),
    ([0.55, 0.59], {0: 'AAG', 1: 'AAA'}),
    ([0.55, 0.61], {0: 'AAG', 1: 'AAG'}),
])
def test_sampling_uses_completion_mass_for_branch_probabilities(draws, expected):
    counter = ComponentModelCounter(
        {0: ('AAA', 'AAG'), 1: ('AAA', 'AAG')},
        [ChoiceFactor([(0, 1, {'AAA': frozenset({'AAA'})})])],
    )
    sampler = FactorSampler(counter, {
        0: {'AAA': 3, 'AAG': 1},
        1: {'AAA': 3, 'AAG': 2},
    })
    rng = Mock()
    rng.random.side_effect = draws

    # The three valid assignments have masses 6, 3 and 2. The first branch
    # therefore has probability 6/11, rather than the unconditioned 3/4.
    assert sampler.sample(rng) == expected
