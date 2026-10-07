import random

from itertools import product

import pytest

from codeine.constraints.repeats import DirectRepeats, InvertedRepeats
from codeine.graph.base import CodonGraph
from codeine.graph.factor_count import ComponentModelCounter
from codeine.graph.factors import ChoiceFactor
from codeine.translation.tables import TranslationTable


@pytest.mark.parametrize('domains,expected', [
    ({}, 1),
    ({0: ('',), 1: ('AAA', 'AAG')}, 2),
    ({1: ()}, 0),
    ({pos: ('AAA', 'AAG') for pos in range(1500)}, 2 ** 1500),
])
def test_count_without_factors(domains, expected):
    counter = ComponentModelCounter(domains, ())

    assert counter.count() == expected


def test_count_impossible_factor():
    counter = ComponentModelCounter({1: ('AAA', 'AAG')}, [ChoiceFactor([])])

    assert counter.count() == 0


def test_count_multiplies_independent_groups_and_free_choices():
    factors = [
        ChoiceFactor([(1, 2, {'AAA': frozenset({'AAA'})})]),
        ChoiceFactor([(3, 4, {'AAA': frozenset({'AAA'})})]),
    ]
    domains = {pos: ('AAA', 'AAG') for pos in range(1, 6)}
    counter = ComponentModelCounter(domains, factors)

    assert counter.count() == 3 * 3 * 2
    assert counter.n_splits > 0

    calls = counter.n_calls

    assert counter.count() == 18
    assert counter.n_calls == calls


def test_count_splits_after_assigning_a_shared_position():
    factors = [
        ChoiceFactor([
            (0, 0, {'AAA': frozenset({'AAA'})}),
            (1, 2, {'AAA': frozenset({'AAA'})}),
        ]),
        ChoiceFactor([
            (0, 0, {'AAA': frozenset({'AAA'})}),
            (3, 4, {'AAA': frozenset({'AAA'})}),
        ]),
    ]
    counter = ComponentModelCounter({pos: ('AAA', 'AAG') for pos in range(5)}, factors)

    # AAA at position 0 leaves two independent rules; AAG makes both harmless.
    assert counter.count() == 3 * 3 + 2 ** 4
    assert counter.n_splits > 0


@pytest.mark.parametrize('seed', range(30))
def test_count_matches_exhaustive_assignments(seed):
    rng = random.Random(seed)
    domains = {pos: ('AAA', 'AAG', 'AAC')[:rng.randint(1, 3)] for pos in range(5)}
    factors = []

    for _ in range(rng.randint(1, 6)):
        relations = []

        for _ in range(rng.randint(1, 3)):
            reference_pos, compare_pos = rng.randrange(5), rng.randrange(5)
            matching = {
                choice: frozenset(other for other in domains[compare_pos] if rng.random() < 0.5)
                for choice in domains[reference_pos]
            }
            relations.append((reference_pos, compare_pos, matching))

        factors.append(ChoiceFactor(relations))

    # Duplicate factors must not change the number of solutions.
    factors.append(factors[0])
    expected = sum(
        not any(factor.rejects(dict(zip(domains, choices))) for factor in factors)
        for choices in product(*domains.values())
    )

    assert ComponentModelCounter(domains, factors).count() == expected


@pytest.mark.parametrize('rna', [False, True])
@pytest.mark.parametrize('repeat_specs', [
    [(DirectRepeats, 3, 0, None)],
    [(InvertedRepeats, 3, 1, 4)],
    [(DirectRepeats, 2, 0, 3), (InvertedRepeats, 2, 0, None)],
])
@pytest.mark.parametrize('aa_seq,context_l,context_r,restricted', [
    ('KK', '', '', False),
    ('LRF', 'AT', 'GC', False),
    ('KK', 'AAA', 'AAA', False),
    ('LRF', 'AT', 'GC', True),
])
def test_repeat_counts_match_flat_compiler(rna, repeat_specs, aa_seq, context_l, context_r, restricted):
    graph = CodonGraph(
        aa_seq,
        context_l=context_l,
        context_r=context_r,
        translation_table=TranslationTable(table_id=1, rna=rna),
    )
    constraints = [kind(length, minimum, maximum) for kind, length, minimum, maximum in repeat_specs]
    factors = []

    for constraint in constraints:
        constraint.link(graph)
        factors.extend(constraint.factors())

    domains = {node.pos: tuple(node.transitions) for node in graph.nodes if node is not graph.end_node}
    fixed_codons = None

    if restricted:
        domains[1] = domains[1][:1]
        fixed_codons = {1: domains[1]}

    # Build the reference separately, so restrictions are applied after factor extraction.
    reference = CodonGraph(
        aa_seq,
        fixed_codons=fixed_codons,
        context_l=context_l,
        context_r=context_r,
        translation_table=graph.tt,
    ).view(constraints=[kind(length, minimum, maximum) for kind, length, minimum, maximum in repeat_specs])

    assert ComponentModelCounter(domains, factors).count() == reference.n_valid_sequences


def test_count_overlapping_unbounded_repeats():
    graph = CodonGraph('MIKEYSASSAFRASMIKEYSASSAFRAS')
    constraint = DirectRepeats(15)
    constraint.link(graph)
    domains = {node.pos: tuple(node.transitions) for node in graph.nodes if node is not graph.end_node}
    counter = ComponentModelCounter(domains, constraint.factors())

    assert counter.count() == 542376004976640
    assert counter.n_calls < 3000
    assert counter.n_splits > 0
