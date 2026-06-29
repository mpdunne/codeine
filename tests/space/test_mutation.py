import itertools
import pytest
import random

from collections import Counter
from scipy.stats import chi2_contingency

from codeine.space.coding import CodingSpace
from codeine.space.mutation import MutationSpace
from codeine.motifs.restriction import RestrictionSite
from codeine.translation.tables import TranslationTable
from codeine.translation.weights import CodonWeights


from tests.data import NORMAL_PROTEINS, DIFFICULT_PROTEINS, ANTIBODIES, LARGE_PROTEINS


def test_mutation_space_requires_cds_in_space():
    space = CodingSpace('MIKEY')

    with pytest.raises(ValueError):
        MutationSpace(space=space, cds='AAAAAAAAAAAA')

    with pytest.raises(ValueError):
        MutationSpace(space=space, cds='')

    with pytest.raises(ValueError):
        MutationSpace(space=space, cds='MIKEY')

    for seq in space:
        _ = MutationSpace(space=space, cds=seq)
        _ = MutationSpace(space=space, cds=seq.lower())
        _ = MutationSpace(space=space, cds=seq.upper())


def test_mutation_space_raises_if_positions_are_invalid():
    space = CodingSpace('MIKEY')

    with pytest.raises(ValueError):
        _ = space.mutants('ATGATTAAAGAATATATG', [0])

    with pytest.raises(ValueError):
        _ = space.mutants('ATGATTAAAGAATATATG', [-1])

    with pytest.raises(ValueError):
        _ = space.mutants('ATGATTAAAGAATATATG', [1, 6])


def test_mutation_space_rejects_invalid_distance_constraints():
    space = CodingSpace('MIKEY')
    cds = 'ATGATTAAAGAATAT'

    with pytest.raises(TypeError, match='min_nts must be an integer'):
        MutationSpace(space, cds, min_nts=1.5)

    with pytest.raises(TypeError, match='max_nts must be an integer'):
        MutationSpace(space, cds, max_nts='2')

    with pytest.raises(ValueError, match='min_nts must be non-negative'):
        MutationSpace(space, cds, min_nts=-1)

    with pytest.raises(ValueError, match='max_codons must be non-negative'):
        MutationSpace(space, cds, max_codons=-1)


def test_mutation_space_rejects_inverted_distance_constraints():
    space = CodingSpace('MIKEY')
    cds = 'ATGATTAAAGAATAT'

    with pytest.raises(ValueError, match='min_nts cannot be greater than max_nts'):
        MutationSpace(space, cds, min_nts=3, max_nts=2)

    with pytest.raises(ValueError, match='min_codons cannot be greater than max_codons'):
        MutationSpace(space, cds, min_codons=2, max_codons=1)


@pytest.mark.parametrize('aa_seq,positions',
                         (
                                 ('MIKEY', [2, 3]),
                                 ('MILDRED', [2, 3, 4]),
                                 ('STEVEN', [1, 2]),
                                 ('WILLIAM', [2, 5, 7]),
                         ))
def test_mutation_space_mutates_only_specified_positions(aa_seq, positions):
    space = CodingSpace(aa_seq)
    ref_cds = space.sample()

    mut = space.mutants(ref_cds, positions)
    sampled_seqs = [mut.sample() for _ in range(1000)]
    assert all(space.contains(s) for s in sampled_seqs)
    assert all(mut.contains(s) for s in sampled_seqs)

    fixed_positions = [pos for pos in range(1, len(aa_seq) + 1) if pos not in positions]

    values_at_fixed_positions = [tuple(seq[(pos - 1) * 3: pos * 3] for pos in fixed_positions) for seq in sampled_seqs]
    values_at_unfixed_positions = [tuple(seq[(pos - 1) * 3: pos * 3] for pos in positions) for seq in sampled_seqs]

    assert len(set(values_at_fixed_positions)) == 1
    assert len(set(values_at_unfixed_positions)) != 1


def test_mutation_space_n_valid_sequences():
    aa_seq = 'MIKEY'
    ref_cds = 'ATGATTAAAGAATAT'

    space = CodingSpace(aa_seq)
    assert space.n_valid_sequences == 24

    mut = space.mutants(ref_cds, [1])
    assert mut.n_valid_variants == 1

    mut = space.mutants(ref_cds, [2])
    assert mut.n_valid_variants == 3

    mut = space.mutants(ref_cds, [1, 2, 3])
    assert mut.n_valid_variants == 6

    assert space.n_valid_sequences == 24


def test_mutation_space_inherits_existing_pins():
    space = CodingSpace('MIKEY')
    cds = space[0]
    space.pin_codons({3: cds[6:9]})
    mutation_space = space.mutants(cds)
    assert mutation_space.view.pinned_codons[3] == [cds[6:9]]
    assert mutation_space.contains(cds)

    space = CodingSpace('MIKEY')
    cds = space[0]
    space.pin_codons({3: cds[6:9]})
    mutation_space = MutationSpace(space, cds)
    assert mutation_space.view.pinned_codons[3] == [cds[6:9]]
    assert mutation_space.contains(cds)


def test_mutation_space_remains_unchanged_if_space_changes():
    space = CodingSpace('MIKEY')
    cds = space[0]
    space.pin_codons({3: cds[6:9]})
    mutation_space = space.mutants(cds)
    space.pin_codons({4: cds[9:12]})
    assert 4 not in mutation_space.view.pinned_codons

    space = CodingSpace('MIKEY')
    cds = space[0]
    space.pin_codons({3: cds[6:9]})
    mutation_space = MutationSpace(space, cds)
    space.pin_codons({4: cds[9:12]})
    assert 4 not in mutation_space.view.pinned_codons


def test_mutation_space_combines_original_pins_with_frozen_positions():
    space = CodingSpace('MIKEY')
    cds = space[0]
    space.pin_codons({3: cds[6:9]})
    mutation_space = space.mutants(cds, free_positions=[2, 3])
    expected_pins = {1: [cds[0:3]], 3: [cds[6:9]], 4: [cds[9:12]], 5: [cds[12:15]]}
    assert mutation_space.view.pinned_codons == expected_pins

    space = CodingSpace('MIKEY')
    cds = space[0]
    space.pin_codons({3: cds[6:9]})
    mutation_space = MutationSpace(space, cds, free_positions=[2, 3])
    expected_pins = {1: [cds[0:3]], 3: [cds[6:9]], 4: [cds[9:12]], 5: [cds[12:15]]}
    assert mutation_space.view.pinned_codons == expected_pins


def test_mutation_space_original_space_untouched_by_freezing():
    space = CodingSpace('MIKEY')
    cds = space[0]
    space.pin_codons({3: cds[6:9]})
    original_pins = dict(space.view.pinned_codons)
    mutation_space = space.mutants(cds, free_positions=[1, 3, 4])
    assert space.view.pinned_codons == original_pins
    assert mutation_space.view.pinned_codons != original_pins

    space = CodingSpace('MIKEY')
    cds = space[0]
    space.pin_codons({3: cds[6:9]})
    original_pins = dict(space.view.pinned_codons)
    mutation_space = MutationSpace(space, cds, free_positions=[1, 3, 4])
    assert space.view.pinned_codons == original_pins
    assert mutation_space.view.pinned_codons != original_pins


def test_mutation_space_defaults_to_all_positions_free():
    space = CodingSpace('MIKEY')
    cds = space[0]
    muts = MutationSpace(space, cds)

    assert muts.free_positions == frozenset(range(1, len(space.view.aa_seq) + 1))
    assert muts.frozen_positions == frozenset()


def test_mutation_space_can_set_free_positions():
    space = CodingSpace('MIKEY')
    cds = space[0]
    muts = MutationSpace(space, cds)

    muts.set_free_positions([1, 3])

    assert muts.free_positions == frozenset({1, 3})
    assert muts.frozen_positions == frozenset({2, 4, 5})


def test_mutation_space_can_freeze_positions():
    space = CodingSpace('MIKEY')
    cds = space[0]
    muts = MutationSpace(space, cds)

    muts.freeze_positions([1])

    assert 1 not in muts.free_positions
    assert 1 in muts.frozen_positions


def test_mutation_space_can_unfreeze_positions():
    space = CodingSpace('MIKEY')
    cds = space[0]
    muts = MutationSpace(space, cds, free_positions=[3])

    muts.unfreeze_positions([1])

    assert muts.free_positions == frozenset({1, 3})


def test_mutation_space_freezes_non_free_positions():
    space = CodingSpace('FF')
    muts = space.mutants('TTTTTT', free_positions=[2])
    seqs = [*muts]
    assert len(seqs) == 2
    assert set(seqs) == {'TTTTTT', 'TTTTTC'}


def test_mutation_space_can_freeze_all():
    space = CodingSpace('MIKEY')
    cds = space[0]
    muts = MutationSpace(space, cds)

    muts.freeze_all()

    assert muts.free_positions == frozenset()
    assert muts.frozen_positions == frozenset(range(1, len(space.view.aa_seq) + 1))


def test_mutation_space_can_unfreeze_all():
    space = CodingSpace('MIKEY')
    cds = space[0]
    muts = MutationSpace(space, cds, free_positions=[])

    muts.unfreeze_all()

    assert muts.free_positions == frozenset(range(1, len(space.view.aa_seq) + 1))
    assert muts.frozen_positions == frozenset()


def test_mutation_space_rejects_invalid_positions():
    space = CodingSpace('MIKEY')
    cds = space[0]
    muts = MutationSpace(space, cds)

    with pytest.raises(ValueError):
        muts.set_free_positions([-1])

    with pytest.raises(ValueError):
        muts.set_free_positions([len(space.view.aa_seq) + 1])


def test_mutation_space_updates_pins_when_free_positions_change():
    space = CodingSpace('MIKEY')
    cds = space[0]

    muts = MutationSpace(space, cds)
    muts.set_free_positions([2])

    expected_pins = {1: [cds[0:3]], 3: [cds[6:9]], 4: [cds[9:12]], 5: [cds[12:15]]}
    assert muts.view.pinned_codons == expected_pins


def test_mutation_space_freeze_all_pins_every_codon():
    space = CodingSpace('MIKEY')
    cds = space.sample()
    muts = MutationSpace(space, cds)

    muts.freeze_all()

    expected_pins = {
        pos: [cds[3 * (pos - 1): 3 * pos]]
        for pos in range(1, len(space.view.aa_seq) + 1)
    }

    assert muts.view.pinned_codons == expected_pins


def test_mutation_space_unfreeze_all_clears_pins():
    space = CodingSpace('MIKEY')
    cds = space.sample()
    muts = MutationSpace(space, cds, free_positions=[])

    muts.unfreeze_all()
    assert space.view.pinned_codons == {}
    assert muts.view.pinned_codons == {}


def test_mutation_space_exposes_graph_properties():
    cw = CodonWeights.ecoli()
    space = CodingSpace(
        'MIKEY',
        codon_restrictions={2: 'ATC'},
        codon_weights=cw,
        context_l='AAA',
        context_r='CCC',
    )
    muts = space.mutants('ATGATCAAAGAGTAT')

    assert muts.aa_seq == space.aa_seq
    assert muts.translation_table is space.translation_table
    assert muts.codon_weights is space.codon_weights
    assert muts.codon_restrictions == space.codon_restrictions
    assert muts.context_l == space.context_l
    assert muts.context_r == space.context_r


def test_mutation_space_exposes_pins_including_frozen_positions():
    space = CodingSpace('MIKEY')
    muts = space.mutants('ATGATCAAAGAGTAT', free_positions=[1, 3, 4, 5])

    assert muts.pinned_codons == {2: ['ATC']}

    muts.unfreeze_positions([2])

    assert muts.pinned_codons == {}


def test_mutation_space_samples_correctly():
    space = CodingSpace('MIKEY')
    cds = space.sample()

    muts = MutationSpace(space, cds, free_positions=[2])
    assert space.view.pinned_codons == {}
    assert set(muts.view.pinned_codons.keys()) == {1, 3, 4, 5}
    variants = [muts.sample() for _ in range(1000)]
    assert len(set(variants)) == muts.n_valid_variants == 3
    for variant in variants:
        assert variant[0:3] == cds[0:3]
        assert variant[6:9] == cds[6:9]
        assert variant[9:12] == cds[9:12]
        assert variant[12:15] == cds[12:15]

    muts = MutationSpace(space, cds, free_positions=[2, 3])
    assert space.view.pinned_codons == {}
    assert set(muts.view.pinned_codons.keys()) == {1, 4, 5}
    variants = [muts.sample() for _ in range(1000)]
    assert len(set(variants)) == muts.n_valid_variants == 6
    for variant in variants:
        assert variant[0:3] == cds[0:3]
        assert variant[9:12] == cds[9:12]
        assert variant[12:15] == cds[12:15]

    muts = MutationSpace(space, cds, free_positions=[5])
    assert space.view.pinned_codons == {}
    assert set(muts.view.pinned_codons.keys()) == {1, 2, 3, 4}
    variants = [muts.sample() for _ in range(1000)]
    assert len(set(variants)) == muts.n_valid_variants == 2
    for variant in variants:
        assert variant[0:3] == cds[0:3]
        assert variant[3:6] == cds[3:6]
        assert variant[6:9] == cds[6:9]
        assert variant[9:12] == cds[9:12]


def test_mutation_space_sample_many():
    space = CodingSpace('MIKEY', seed=88)
    reference = space[0]
    muts = space.mutants(reference)
    seqs = muts.sample(n=10)
    assert len(seqs) == 10
    assert all(seq in space for seq in seqs)


def test_mutation_space_does_not_modify_original_space():
    space = CodingSpace('MIKEY')
    cds = space[0]

    muts = space.mutants(cds, free_positions=[2])
    muts.freeze_all()

    assert space.view.pinned_codons == {}
    assert muts.view.pinned_codons != {}


def test_mutation_space_has_no_distance_constraints_by_default():
    space = CodingSpace('MIKEY')
    cds = space.sample()
    muts = MutationSpace(space, cds)

    assert muts.min_nts is None
    assert muts.max_nts is None
    assert muts.min_codons is None
    assert muts.max_codons is None
    assert not muts.has_distance_constraints


def test_mutation_space_accepts_distance_constraints_in_constructor():
    space = CodingSpace('MIKEY')
    cds = space.sample()
    muts = MutationSpace(space, cds, min_nts=3, max_nts=4, min_codons=1, max_codons=2)

    assert muts.min_nts == 3
    assert muts.max_nts == 4
    assert muts.min_codons == 1
    assert muts.max_codons == 2
    assert muts.has_distance_constraints


def test_mutation_space_accepts_exact_nt_distance_in_constructor():
    space = CodingSpace('MIKEY')
    cds = space.sample()
    muts = MutationSpace(space, cds, min_nts=3, max_nts=5)

    assert muts.min_nts == 3
    assert muts.max_nts == 5


def test_mutation_space_accepts_exact_codon_distance_in_constructor():
    space = CodingSpace('MIKEY')
    cds = space.sample()
    muts = MutationSpace(space, cds, min_codons=3, max_codons=5)

    assert muts.min_codons == 3
    assert muts.max_codons == 5


def test_set_distance_sets_nt_range():
    space = CodingSpace('MIKEY')
    cds = space.sample()
    muts = MutationSpace(space, cds)

    muts.set_distance_constraints(min_nts=3, max_nts=4)

    assert muts.min_nts == 3
    assert muts.max_nts == 4
    assert muts.min_codons is None
    assert muts.max_codons is None


def test_set_distance_sets_codon_range():
    space = CodingSpace('MIKEY')
    cds = space.sample()
    muts = MutationSpace(space, cds)

    muts.set_distance_constraints(min_codons=3, max_codons=10)

    assert muts.min_codons == 3
    assert muts.max_codons == 10
    assert muts.min_nts is None
    assert muts.max_nts is None


def test_set_distance_can_set_nt_and_codon_constraints_together():
    space = CodingSpace('MIKEY')
    cds = space.sample()
    muts = MutationSpace(space, cds)

    muts.set_distance_constraints(min_nts=3, max_nts=6, min_codons=2, max_codons=3)

    assert muts.min_nts == 3
    assert muts.max_nts == 6
    assert muts.min_codons == 2
    assert muts.max_codons == 3


@pytest.mark.parametrize(
    'kwargs',
    [
        {'min_nts': -1},
        {'max_nts': -1},
        {'min_codons': -1},
        {'max_codons': -1},
    ],
)
def test_set_distance_rejects_negative_distances(kwargs):
    space = CodingSpace('MIKEY')
    cds = space.sample()
    muts = MutationSpace(space, cds)
    with pytest.raises(ValueError):
        muts.set_distance_constraints(**kwargs)


@pytest.mark.parametrize(
    'kwargs',
    [
        {'min_nts': 1.5},
        {'max_nts': 1.5},
        {'min_codons': 1.5},
        {'max_codons': 1.5},
    ],
)
def test_set_distance_rejects_noninteger_distances(kwargs):
    space = CodingSpace('MIKEY')
    cds = space.sample()
    muts = MutationSpace(space, cds)
    with pytest.raises(TypeError):
        muts.set_distance_constraints(**kwargs)


def test_set_distance_rejects_min_nts_greater_than_max_nts():
    space = CodingSpace('MIKEY')
    cds = space.sample()
    muts = MutationSpace(space, cds)
    with pytest.raises(ValueError):
        muts.set_distance_constraints(min_nts=5, max_nts=3)


def test_set_distance_rejects_min_codons_greater_than_max_codons():
    space = CodingSpace('MIKEY')
    cds = space.sample()
    muts = MutationSpace(space, cds)
    with pytest.raises(ValueError):
        muts.set_distance_constraints(min_codons=5, max_codons=3)


def test_set_distance_with_no_args_clears_distance_constraints():
    space = CodingSpace('MIKEY')
    cds = space.sample()
    muts = MutationSpace(space, cds, min_nts=1, max_nts=10, min_codons=1, max_codons=5)

    muts.set_distance_constraints()
    assert muts.min_nts is None
    assert muts.max_nts is None
    assert muts.min_codons is None
    assert muts.max_codons is None
    assert not muts.has_distance_constraints


def test_clear_distance_constraints():
    space = CodingSpace('MIKEY')
    cds = space.sample()
    muts = MutationSpace(space, cds, min_nts=1, max_nts=10, min_codons=1, max_codons=5)

    muts.clear_distance_constraints()
    assert muts.min_nts is None
    assert muts.max_nts is None
    assert muts.min_codons is None
    assert muts.max_codons is None
    assert not muts.has_distance_constraints


standard_table = TranslationTable().aa_to_codons


def sample_naive(aa_seq, n=1):
    seqs = []

    for _ in range(n):
        codons = []
        for aa in aa_seq:
            codon = random.choice(standard_table[aa])
            codons.append(codon)
        seq = ''.join(codons)
        seqs.append(seq)

    if n == 1:
        return seqs[0]
    else:
        return seqs


def enumerate_naive(aa_seq):
    all_codons = [standard_table[aa] for aa in aa_seq]
    all_seqs = [''.join(entry) for entry in itertools.product(*all_codons)]
    return all_seqs


def get_nt_diffs(cds1, cds2):
    if len(cds1) != len(cds2):
        raise ValueError('Sequences must be same length')

    diffs = 0
    for nt1, nt2 in zip(cds1.upper(), cds2.upper()):
        if nt1 != nt2:
            diffs += 1

    return diffs


def get_codon_diffs(cds1, cds2):
    if len(cds1) != len(cds2):
        raise ValueError('Sequences must be same length')

    if len(cds1) % 3 != 0:
        raise ValueError('Sequence lengths must both be multiples of 3')

    diffs = 0
    for pos in range(1, 1 + (len(cds1) // 3)):
        codon1 = cds1[3 * (pos - 1): 3 * pos]
        codon2 = cds2[3 * (pos - 1): 3 * pos]

        if codon1 != codon2:
            diffs += 1

    return diffs


def test_max_nts_zero_only_allows_reference_sequence():
    space = CodingSpace('MIKEY')
    muts = MutationSpace(space, 'ATGATTAAAGAATAT')

    muts.set_distance_constraints(max_nts=0)

    assert [*muts.enumerate()] == ['ATGATTAAAGAATAT']


def test_exact_nt_distance_one():
    space = CodingSpace('MIKEY')
    muts = MutationSpace(space, 'ATGATTAAAGAATAT')

    muts.set_distance_constraints(min_nts=1, max_nts=1)

    expected = [
        'ATGATTAAAGAATAC',
        'ATGATTAAAGAGTAT',
        'ATGATTAAGGAATAT',
        'ATGATCAAAGAATAT',
        'ATGATAAAAGAATAT'
    ]

    observed = [*muts]

    assert len(observed) == len(expected)
    assert set(observed) == set(expected)


def test_exact_codon_distance_zero():
    space = CodingSpace('MIKEY')
    muts = MutationSpace(space, 'ATGATTAAAGAATAT')

    muts.set_distance_constraints(min_codons=0, max_codons=0)

    assert [*muts.enumerate()] == ['ATGATTAAAGAATAT']


def test_exact_codon_distance_one():
    space = CodingSpace('MIKEY')
    muts = MutationSpace(space, 'ATGATTAAAGAATAT')

    muts.set_distance_constraints(min_codons=1, max_codons=1)

    expected = [
        'ATGATTAAAGAATAC',
        'ATGATTAAAGAGTAT',
        'ATGATTAAGGAATAT',
        'ATGATCAAAGAATAT',
        'ATGATAAAAGAATAT'
    ]

    observed = [*muts.enumerate()]

    assert len(observed) == len(expected)
    assert set(observed) == set(expected)


def test_distance_constraints_reduce_count():
    space = CodingSpace('MIKEY')
    cds = space[0]

    unconstrained = MutationSpace(space, cds)
    constrained = MutationSpace(space, cds, max_nts=2)

    assert unconstrained.n_valid_variants == 24
    assert constrained.n_valid_variants == 15


def test_distance_constraints_allow_contains_only_constrained_cdss():
    space = CodingSpace('MIKEY')
    all_seqs = [*space]
    reference = all_seqs[0]

    muts = MutationSpace(space, reference)
    muts.set_distance_constraints(max_nts=3)

    for seq in all_seqs:
        if get_nt_diffs(seq, reference) <= 3:
            assert seq in muts
        else:
            assert seq not in muts


def test_distance_constraints_affect_sampling_and_enumeration():
    space = CodingSpace('SASSAFRAS')
    assert space.n_valid_sequences == 995328

    reference = 'TCTGCTTCTTCTGCTTTTCGTGCTTCT'
    muts = MutationSpace(space, reference)

    muts.set_distance_constraints(max_nts=0)
    assert muts.n_valid_variants == 1
    assert len([*muts]) == 1
    for _ in range(100):
        assert muts.sample() == reference

    muts.set_distance_constraints(min_nts=1, max_nts=1)
    assert muts.n_valid_variants == 25
    assert len([*muts]) == 25
    for _ in range(100):
        assert get_nt_diffs(reference, muts.sample()) == 1

    muts.set_distance_constraints(max_nts=3)
    assert muts.n_valid_variants == 2208
    assert len([*muts]) == 2208
    for _ in range(100):
        assert get_nt_diffs(reference, muts.sample()) <= 3

    muts.set_distance_constraints(max_codons=0)
    assert muts.n_valid_variants == 1
    assert len([*muts]) == 1
    for _ in range(100):
        assert muts.sample() == reference

    muts.set_distance_constraints(min_codons=1, max_codons=1)
    assert muts.n_valid_variants == 35
    assert len([*muts]) == 35
    for _ in range(100):
        assert get_codon_diffs(reference, muts.sample()) == 1

    muts.set_distance_constraints(max_codons=3)
    assert muts.n_valid_variants == 5276
    assert len([*muts]) == 5276
    for _ in range(100):
        assert get_codon_diffs(reference, muts.sample()) <= 3


def test_clear_distance_constraints_restores_original_count():
    space = CodingSpace('MIKEY')
    cds = space[0]

    muts = MutationSpace(space, cds)
    assert muts.n_valid_variants == 24

    muts.set_distance_constraints(max_nts=2)
    assert muts.n_valid_variants == 15

    muts.clear_distance_constraints()
    assert muts.n_valid_variants == 24


def test_banned_sequences_work_with_distance_constraints():
    space = CodingSpace('SASSAFRAS', context_l='AAAAAA', context_r='GGGGGG')
    assert space.n_valid_sequences == 995328

    banned = ['TTTT', 'AAATCT', 'TCTGGG']
    space = CodingSpace('SASSAFRAS', context_l='AAAAAA', context_r='GGGGGG', forbidden_motifs=banned)
    assert space.n_valid_sequences == 604800

    reference = 'TCCGCTTCTTCTGCTTTCCGTGCTTCC'

    muts = space.mutants(reference)
    assert muts.n_valid_variants == 604800

    muts.set_distance_constraints(max_nts=5)
    assert muts.n_valid_variants == 23155

    muts.set_distance_constraints(min_codons=5, max_codons=5)
    assert muts.n_valid_variants == 60866


def has_banned_sequence(cds, banned_sequences=(), context_l='', context_r=''):
    full_seq = context_l + cds + context_r
    return any(banned in full_seq for banned in banned_sequences)


def enumerate_naive_mutants(
    aa_seq,
    reference,
    banned_sequences=(),
    context_l='',
    context_r='',
    min_nts=None,
    max_nts=None,
    min_codons=None,
    max_codons=None,
):
    seqs = []

    for cds in enumerate_naive(aa_seq):
        nt_diffs = get_nt_diffs(cds, reference)
        codon_diffs = get_codon_diffs(cds, reference)

        if min_nts is not None and nt_diffs < min_nts:
            continue
        if max_nts is not None and nt_diffs > max_nts:
            continue
        if min_codons is not None and codon_diffs < min_codons:
            continue
        if max_codons is not None and codon_diffs > max_codons:
            continue
        if has_banned_sequence(cds, banned_sequences, context_l, context_r):
            continue

        seqs.append(cds)

    return seqs


SHORT_AA_SEQUENCES = ['K', 'KK', 'MIKEY', 'SALLY', 'SASSY']

CONTEXTS = [
    ('', ''),
    ('AAAAAA', ''),
    ('', 'GGGGGG'),
    ('AAAAAA', 'GGGGGG'),
]

BANNED_SEQUENCE_SETS = [
    [],
    ['AAA'],
    ['TTTT'],
    ['AAATCT'],
    ['TCTGGG'],
    ['TTTT', 'AAATCT', 'TCTGGG'],
]

DISTANCE_CONSTRAINTS = [
    {},
    {'max_nts': 0},
    {'max_nts': 1},
    {'max_nts': 3},
    {'max_nts': 10},
    {'min_nts': 1, 'max_nts': 3},
    {'min_nts': 5, 'max_nts': 20},
    {'max_codons': 0},
    {'max_codons': 1},
    {'max_codons': 3},
    {'max_codons': 10},
    {'min_codons': 5, 'max_codons': 20},
]


@pytest.mark.parametrize('aa_seq', SHORT_AA_SEQUENCES)
@pytest.mark.parametrize('context', CONTEXTS)
@pytest.mark.parametrize('banned_sequences', BANNED_SEQUENCE_SETS)
@pytest.mark.parametrize('distance_kwargs', DISTANCE_CONSTRAINTS)
def test_mutation_space_matches_naive_combinatorial(
    aa_seq,
    context,
    banned_sequences,
    distance_kwargs,
):
    context_l, context_r = context

    space = CodingSpace(
        aa_seq,
        context_l=context_l,
        context_r=context_r,
        forbidden_motifs=banned_sequences,
    )

    if space.n_valid_sequences == 0:
        return

    reference = space[0]
    muts = MutationSpace(space, reference, **distance_kwargs)

    expected = enumerate_naive_mutants(
        aa_seq,
        reference,
        banned_sequences=banned_sequences,
        context_l=context_l,
        context_r=context_r,
        **distance_kwargs,
    )

    assert [*muts] == expected
    assert muts.n_valid_variants == len(expected)

    expected_set = set(expected)
    for cds in enumerate_naive(aa_seq):
        assert (cds in muts) == (cds in expected_set)


REAL_CONTEXTS = [
    ('', ''),

    # Cloning / stop
    ('GCCACC', ''),
    ('', 'TAA'),
    ('GCCACC', 'TAA'),

    # Restriction-sites
    ('GAATTC', 'AAGCTT'),    # EcoRI / HindIII
    ('GGATCC', 'CTCGAG'),    # BamHI / XhoI
]


REAL_BANNED_CASES = [
    [],
    [RestrictionSite.EcoRI, RestrictionSite.BamHI],
    [RestrictionSite.EcoRI, RestrictionSite.BamHI, RestrictionSite.XhoI, RestrictionSite.HindIII],
]

# REAL_PROTEINS = {**NORMAL_PROTEINS, **DIFFICULT_PROTEINS, **ANTIBODIES, **LARGE_PROTEINS}
REAL_PROTEINS = {**NORMAL_PROTEINS, **LARGE_PROTEINS}


@pytest.mark.parametrize('name, aa_seq', REAL_PROTEINS.items())
@pytest.mark.parametrize('context', REAL_CONTEXTS)
@pytest.mark.parametrize('banned_sequences', REAL_BANNED_CASES)
@pytest.mark.parametrize('distance_kwargs', DISTANCE_CONSTRAINTS)
def test_mutation_space_samples_real_proteins(
    name,
    aa_seq,
    context,
    banned_sequences,
    distance_kwargs,
):
    context_l, context_r = context

    space = CodingSpace(
        aa_seq,
        context_l=context_l,
        context_r=context_r,
        forbidden_motifs=banned_sequences,
    )

    if space.n_valid_sequences == 0:
        return

    reference = space.sample()
    muts = MutationSpace(space, reference, **distance_kwargs)

    if muts.n_valid_variants == 0:
        return

    for _ in range(100):
        cds = muts.sample()

        assert cds in muts
        assert cds in space

        nt_diffs = get_nt_diffs(cds, reference)
        codon_diffs = get_codon_diffs(cds, reference)

        if muts.min_nts is not None:
            assert nt_diffs >= muts.min_nts
        if muts.max_nts is not None:
            assert nt_diffs <= muts.max_nts
        if muts.min_codons is not None:
            assert codon_diffs >= muts.min_codons
        if muts.max_codons is not None:
            assert codon_diffs <= muts.max_codons


def helper_mutation_counts_by_block(seqs, ref, block_size):
    counts = Counter()

    for seq in seqs:
        for codon_i in range(len(ref) // 3):
            start = codon_i * 3
            if seq[start:start + 3] != ref[start:start + 3]:
                block = codon_i // block_size
                counts[block] += 1

    return counts


@pytest.mark.parametrize('aa_seq', (
    'MIKEY' * 10,
    'MIKEY' * 100,
))
@pytest.mark.parametrize('banned', (
    (),
    ('GAATTC', 'GGATCC'),
))
@pytest.mark.parametrize('distance_constraints', (
    dict(min_nts=2, max_nts=20),
    dict(min_nts=10, max_nts=20),
    dict(min_codons=2, max_codons=20),
    dict(min_codons=8, max_codons=15),
    dict(min_nts=10, max_nts=20, min_codons=8, max_codons=15),
))
def test_mutation_sampling_is_even_across_sequence(aa_seq, banned, distance_constraints):
    space = CodingSpace(aa_seq, context_l='aaa', context_r='ttt', forbidden_motifs=banned)

    ref = space[0]

    muts = space.mutants(ref)
    muts.set_distance_constraints(**distance_constraints)

    n = 1000
    seqs = [muts.sample() for _ in range(n)]

    counts = helper_mutation_counts_by_block(seqs, ref, block_size=5)

    n_blocks = len(aa_seq) // 5
    total = sum(counts.values())
    expected = total / n_blocks

    first_half = sum(counts[i] for i in range(n_blocks // 2))
    second_half = sum(counts[i] for i in range(n_blocks // 2, n_blocks))

    assert abs(first_half - second_half) <= 0.05 * total

    for i in range(n_blocks):
        assert abs(counts[i] - expected) <= 0.5 * expected


def helper_codon_counts_by_block_position(seqs, block_size=5):
    counts = {within_block: [] for within_block in range(block_size)}

    n_codons = len(seqs[0]) // 3
    n_blocks = n_codons // block_size

    for within_block in range(block_size):
        counts[within_block] = [Counter() for _ in range(n_blocks)]

    for seq in seqs:
        for codon_i in range(n_codons):
            block = codon_i // block_size
            within_block = codon_i % block_size
            start = codon_i * 3
            counts[within_block][block][seq[start:start + 3]] += 1

    return counts


@pytest.mark.parametrize('aa_seq', (
    'MIKEY' * 10,
    'MIKEY' * 100,
))
@pytest.mark.parametrize('banned', (
    (),
    ('GAATTC', 'GGATCC'),
))
@pytest.mark.parametrize('distance_constraints', (
    dict(min_nts=2, max_nts=20),
    dict(min_nts=10, max_nts=20),
    dict(min_codons=2, max_codons=20),
    dict(min_codons=8, max_codons=15),
    dict(min_nts=10, max_nts=20, min_codons=8, max_codons=15),
))
def test_mutation_codon_distributions_are_stable_across_sequence(aa_seq, banned, distance_constraints):

    space = CodingSpace(aa_seq, context_l='aaa', context_r='ttt', forbidden_motifs=banned, seed=8765309)

    ref = space[0]

    muts = space.mutants(ref)
    muts.set_distance_constraints(**distance_constraints)

    n = 1000
    seqs = [muts.sample() for _ in range(n)]

    counts = helper_codon_counts_by_block_position(seqs)
    pvalues = []

    # skip M, only one codon
    for within_block in range(1, 5):
        block_counts = counts[within_block]
        codons = sorted({
            codon
            for count in block_counts
            for codon in count
        })

        table = [
            [count.get(codon, 0) for codon in codons]
            for count in block_counts
        ]

        _, pvalue, _, _ = chi2_contingency(table)
        pvalues.append(pvalue)

    assert sum(p >= 0.001 for p in pvalues) / len(pvalues) >= 0.99
    assert min(pvalues, default=1.0) >= 1e-8
