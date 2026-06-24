import pickle
import random

import pytest

from Bio.Seq import Seq

from codeine.translation.weights import CodonWeights
from codeine.space.coding import CodingSpace
from codeine.motifs.restriction import RestrictionSite


@pytest.mark.parametrize('aa_seq', ('MIKEY', 'MILDRED', 'STEVEN', 'WILLIAM'))
def test_coding_space_sequences_translate_correctly(aa_seq):
    space = CodingSpace(aa_seq=aa_seq)
    for _ in range(1000):
        cds = space.sample()
        translated = Seq(cds).translate()
        assert translated == aa_seq


def test_coding_space_fixed_codons_are_fixed():
    space = CodingSpace('MIKEY', codon_restrictions={2: 'ATA'})

    for _ in range(1000):
        cds = space.sample()
        codons = [cds[i:i + 3] for i in range(0, len(cds), 3)]

        assert codons[1] == 'ATA'
        assert str(Seq(cds).translate()) == 'MIKEY'


def test_coding_space_generates_different_sequences():
    space = CodingSpace('MIKEY')
    generated = [space.sample() for _ in range(1000)]
    assert len(generated) > 1


def test_coding_space_generates_all_cdss_for_small_aa_seq():
    space = CodingSpace('MF')
    generated = {space.sample() for _ in range(1000)}
    expected = {'ATGTTT', 'ATGTTC'}
    assert generated == expected


def test_coding_space_pinned_codons_are_fixed():
    space = CodingSpace('MIKEY')
    space.pin_codons({3: 'AAA'})
    for _ in range(100):
        cds = space.sample()
        codons = [cds[i:i + 3] for i in range(0, len(cds), 3)]
        assert codons[2] == 'AAA'


def test_coding_space_unpin_codons_restores_sampling():
    space = CodingSpace('MIKEY')
    space.pin_codons({3: 'AAA'})
    space.unpin_codons([3])
    sampled = set()
    for _ in range(100):
        cds = space.sample()
        codons = [cds[i:i + 3] for i in range(0, len(cds), 3)]
        sampled.add(codons[2])
    assert sampled == {'AAA', 'AAG'}


def test_coding_space_clear_pins_restores_sampling():
    space = CodingSpace('MIKEY')

    space.pin_codons({3: 'AAA'})
    space.clear_pins()

    sampled = set()
    for _ in range(100):
        cds = space.sample()
        codons = [cds[i:i + 3] for i in range(0, len(cds), 3)]
        sampled.add(codons[2])

    assert sampled == {'AAA', 'AAG'}


def test_coding_space_rejects_invalid_pin():
    space = CodingSpace('MIKEY')
    with pytest.raises(ValueError):
        space.pin_codons({3: 'GCT'})


def test_coding_space_rejects_out_of_range_pin():
    space = CodingSpace('MIKEY')

    with pytest.raises(ValueError):
        space.pin_codons({0: 'ATG'})

    with pytest.raises(ValueError):
        space.pin_codons({6: 'ATG'})


def test_coding_space_sample_excludes_context_by_default():
    space = CodingSpace(
        aa_seq='MF',
        context_l='AAAA',
        context_r='CCCC',
    )

    cds = space.sample()

    assert cds in {'ATGTTT', 'ATGTTC'}
    assert not cds.startswith('AAAA')
    assert not cds.endswith('CCCC')


def test_coding_space_mutants_raises_if_seq_is_invalid():
    space = CodingSpace('MIKEY')
    
    with pytest.raises(ValueError):
        _ = space.mutants('', [1, 2])

    with pytest.raises(ValueError):
        _ = space.mutants('ATG', [1, 2])

    with pytest.raises(ValueError):
        _ = space.mutants('ATTATTAAAGAATAT', [1, 2])

    with pytest.raises(ValueError):
        _ = space.mutants('ATGATTAAAGAATATATG', [1, 2])


def test_coding_space_mutants_returns_mutation_space():
    space = CodingSpace('MIKEY')
    muts = space.mutants('ATGATTAAAGAATAT', free_positions=[2])

    assert muts.cds == 'ATGATTAAAGAATAT'
    assert muts.free_positions == {2}


def test_coding_space_mutants_does_not_change_base_space():
    aa_seq = 'MIKEY'
    ref_cds = 'ATGATTAAAGAATAT'

    space = CodingSpace(aa_seq)

    sampled_seqs = [space.sample() for _ in range(1000)]
    assert not all(s[6:] == ref_cds[6:] for s in sampled_seqs)
    assert not all(s[3:6] == ref_cds[3:6] for s in sampled_seqs)

    mut = space.mutants(ref_cds, [2])
    sampled_seqs = [mut.sample() for _ in range(1000)]
    assert all(s[6:] == ref_cds[6:] for s in sampled_seqs)
    assert not all(s[3:6] == ref_cds[3:6] for s in sampled_seqs)

    sampled_seqs = [space.sample() for _ in range(1000)]
    assert not all(s[6:] == ref_cds[6:] for s in sampled_seqs)
    assert not all(s[3:6] == ref_cds[3:6] for s in sampled_seqs)


def test_coding_space_mutants_passes_constraints():
    space = CodingSpace('MIKEY')
    muts = space.mutants(
        'ATGATTAAAGAATAT',
        free_positions=[2, 3],
        min_nts=1,
        max_nts=3,
        min_codons=1,
        max_codons=2,
    )

    assert set(muts.free_positions) == set([2, 3])
    assert muts.min_nts == 1
    assert muts.max_nts == 3
    assert muts.min_codons == 1
    assert muts.max_codons == 2


def test_sequence_space_getitem():
    space = CodingSpace('MM')
    assert space[0] == 'ATGATG'


def test_view_iter():
    space = CodingSpace('MIKEY')
    seqs = [*space]
    assert len(seqs) == len(set(seqs)) == 24


def test_coding_space_enumerate():
    space = CodingSpace('F')
    assert list(space.enumerate()) == ['TTT', 'TTC']


def test_coding_space_contains():
    space = CodingSpace('F')
    assert space.contains('TTT')
    assert space.contains('TTC')
    assert not space.contains('ATG')


def test_coding_space_contains():
    space = CodingSpace('MIKEY')
    for _ in range(100):
        seq = space.sample()
        assert seq in space
        assert seq + 'ATG' not in space


def test_forbidden_motifs_are_stored():
    space = CodingSpace('MIKEY', forbidden_motifs=[RestrictionSite.EcoRI, 'AAAA'])
    assert space.forbidden_sequences == ('AAAA', 'GAATTC')


def test_max_homopolymer_is_stored():
    space = CodingSpace('MIKEY', max_homopolymer=4)
    assert space.max_homopolymer == 4


def test_max_homopolymer_is_expanded_correctly():
    space = CodingSpace('MIKEY', max_homopolymer=None)
    assert not space.forbidden_sequences

    space = CodingSpace('MIKEY', max_homopolymer=4)
    assert all(nt * 5 in space.forbidden_sequences for nt in 'ACGT')


def test_mixed_restrictions():
    space = CodingSpace('MIKEY', max_homopolymer=4, forbidden_motifs=[RestrictionSite.BsaI, 'GGTTCC'])
    assert set(space.forbidden_sequences) == {'GAGACC', 'GGTCTC', 'GGTTCC', 'AAAAA', 'CCCCC', 'GGGGG', 'TTTTT'}


def test_forbidden_motifs_repr():
    space = CodingSpace('MIKEY', forbidden_motifs=[RestrictionSite.EcoRI, 'AAAA'],)

    text = repr(space)
    assert 'Forbidden motifs:' in text
    assert 'EcoRI' in text
    assert 'GAATTC' in text
    assert 'AAAA' in text


def test_max_homopolymer_repr():
    space = CodingSpace('MIKEY', max_homopolymer=4)

    text = repr(space)
    assert 'Maximum homopolymer' in text
    assert '4' in text


def test_coding_space_seed_and_rng_cannot_both_be_provided():
    with pytest.raises(ValueError):
        CodingSpace('MIKEY', seed=420, rng=random.Random(69))


def test_coding_space_seed_consistent():
    samples_by_rep = []

    for _ in range(10):
        space = CodingSpace('MIKEY', seed=8675309)
        samples = [space.sample() for _ in range(100)]
        samples_by_rep.append(samples)

    assert all(samples == samples_by_rep[0] for samples in samples_by_rep)


def test_coding_space_rng_consistent():
    samples_by_rep = []

    for _ in range(10):
        rng = random.Random(8675309)
        space = CodingSpace('MIKEY', rng=rng)
        samples = [space.sample() for _ in range(100)]
        samples_by_rep.append(samples)

    assert all(samples == samples_by_rep[0] for samples in samples_by_rep)


def test_coding_space_no_seed_not_consistent():
    samples_by_rep = []

    for _ in range(10):
        space = CodingSpace('MIKEY')
        samples = [space.sample() for _ in range(100)]
        samples_by_rep.append(samples)

    assert not all(samples == samples_by_rep[0] for samples in samples_by_rep)


def test_coding_space_pickle_preserves_random_state():
    space = CodingSpace('MIKEY', seed=8675309)

    _ = [space.sample() for _ in range(100)]

    loaded = pickle.loads(pickle.dumps(space))

    assert [loaded.sample() for _ in range(100)] == [space.sample() for _ in range(100)]


def test_coding_space_pickle_preserves_pins():
    space = CodingSpace('MIKEY', seed=8675309)
    space.pin_codons({2: 'ATC'})

    loaded = pickle.loads(pickle.dumps(space))

    assert loaded.view.pinned_codons == space.view.pinned_codons
    assert loaded.n_valid_sequences == space.n_valid_sequences
    assert [*loaded.enumerate()] == [*space.enumerate()]


def test_coding_space_pickle_preserves_constraints():
    space = CodingSpace(
        'MIKEY',
        codon_restrictions={2: 'ATC'},
        forbidden_motifs=['AAAA'],
        max_homopolymer=4,
        seed=8675309,
    )

    loaded = pickle.loads(pickle.dumps(space))

    assert loaded.forbidden_motifs == space.forbidden_motifs
    assert loaded.max_homopolymer == space.max_homopolymer
    assert loaded.forbidden_sequences == space.forbidden_sequences
    assert loaded.n_valid_sequences == space.n_valid_sequences


def test_coding_space_pickle_preserves_future_sampling():
    space = CodingSpace('MIKEY', seed=8675309)

    _ = [space.sample() for _ in range(100)]

    loaded = pickle.loads(pickle.dumps(space))

    expected = [space.sample() for _ in range(1000)]
    observed = [loaded.sample() for _ in range(1000)]

    assert observed == expected


def test_coding_space_pickle_preserves_future_sampling_with_pins():
    space = CodingSpace('MIKEY', seed=8675309)
    space.pin_codons({2: 'ATC'})

    _ = [space.sample() for _ in range(100)]

    loaded = pickle.loads(pickle.dumps(space))

    expected = [space.sample() for _ in range(1000)]
    observed = [loaded.sample() for _ in range(1000)]

    assert observed == expected


def test_coding_space_pickle_preserves_initial_sampling():
    space = CodingSpace('MIKEY', seed=8675309)

    loaded = pickle.loads(pickle.dumps(space))

    expected = [space.sample() for _ in range(1000)]
    observed = [loaded.sample() for _ in range(1000)]

    assert observed == expected


def test_coding_space_save_load(tmp_path):
    path = tmp_path / 'space.pkl'

    space = CodingSpace('MIKEY', seed=8675309)
    space.save(path)

    loaded = CodingSpace.load(path)

    assert loaded.n_valid_sequences == space.n_valid_sequences
    assert [*loaded.enumerate()] == [*space.enumerate()]


def test_coding_space_save_load_preserves_random_state(tmp_path):
    path = tmp_path / 'space.pkl'

    space = CodingSpace('MIKEY', seed=8675309)

    _ = [space.sample() for _ in range(100)]

    space.save(path)
    loaded = CodingSpace.load(path)

    expected = [space.sample() for _ in range(1000)]
    observed = [loaded.sample() for _ in range(1000)]

    assert observed == expected


def test_coding_space_save_load_preserves_pins(tmp_path):
    path = tmp_path / 'space.pkl'

    space = CodingSpace('MIKEY', seed=8675309)
    space.pin_codons({2: 'ATC'})

    space.save(path)
    loaded = CodingSpace.load(path)

    assert loaded.view.pinned_codons == space.view.pinned_codons
    assert loaded.n_valid_sequences == space.n_valid_sequences
    assert [*loaded.enumerate()] == [*space.enumerate()]


def test_coding_space_exposes_graph_properties():
    cw = CodonWeights.ecoli()
    space = CodingSpace(
        'MIKEY',
        codon_restrictions={2: 'ATC'},
        codon_weights=cw,
        context_l='AAA',
        context_r='CCC',
    )

    assert space.aa_seq == 'MIKEY'
    assert space.translation_table is space.view.translation_table
    assert space.codon_weights is space.view.codon_weights
    assert space.codon_weights is cw
    assert space.codon_restrictions == space.view.codon_restrictions
    assert space.context_l == 'AAA'
    assert space.context_r == 'CCC'


def test_coding_space_exposes_pins():
    space = CodingSpace('MIKEY')

    space.pin_codons({2: 'ATC'})
    assert space.pinned_codons == {2: ['ATC']}

    space.clear_pins()
    assert space.pinned_codons == {}


def test_coding_space_can_set_forbidden_motifs():
    space = CodingSpace('MIKEY')

    space.set_forbidden_motifs(['AAA'])

    assert space.forbidden_motifs == ['AAA']
    assert set(space.forbidden_sequences) == {'AAA'}
    assert set(space.view.banned_sequences) == {'AAA'}

    space.clear_forbidden_motifs()

    assert space.forbidden_motifs is None
    assert set(space.forbidden_sequences) == set()
    assert set(space.view.banned_sequences) == set()


def test_coding_space_can_set_max_homopolymer():
    space = CodingSpace('MIKEY')
    space.set_max_homopolymer(2)

    assert space.max_homopolymer == 2
    assert set(space.forbidden_sequences) == {'AAA', 'CCC', 'GGG', 'TTT'}
    assert set(space.view.banned_sequences) == {'AAA', 'CCC', 'GGG', 'TTT'}

    space.clear_max_homopolymer()

    assert space.max_homopolymer is None
    assert set(space.forbidden_sequences) == set()
    assert set(space.view.banned_sequences) == set()


def test_coding_space_forbidden_motifs_and_max_homopolymer_combine():
    space = CodingSpace('MIKEY')

    space.set_forbidden_motifs(['GAG'])
    space.set_max_homopolymer(2)

    assert set(space.forbidden_sequences) == {'AAA', 'CCC', 'GGG', 'TTT', 'GAG'}
    assert set(space.view.banned_sequences) == {'AAA', 'CCC', 'GGG', 'TTT', 'GAG'}

    space.clear_forbidden_motifs()
    assert set(space.forbidden_sequences) == {'AAA', 'CCC', 'GGG', 'TTT'}
    assert set(space.view.banned_sequences) == {'AAA', 'CCC', 'GGG', 'TTT'}

    space.clear_max_homopolymer()
    assert set(space.forbidden_sequences) == set()
    assert set(space.view.banned_sequences) == set()


def test_forbidden_sequences_is_read_only():
    space = CodingSpace('MIKEY')

    space.set_forbidden_motifs('GGG')
    assert isinstance(space.forbidden_sequences, tuple)
    assert space.forbidden_sequences == ('GGG',)

    with pytest.raises(AttributeError):
        space.forbidden_sequences = ('AAA',)


