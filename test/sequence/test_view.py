import pickle
import pytest

from itertools import product
from unittest.mock import MagicMock

from codeine.sequence.graph import CodonGraph


def test_view_can_pin_codons():
    view = CodonGraph('MIKEY').view()
    view.pin_codons({3: 'AAA'})
    assert view.pinned_codons[3] == ['AAA']


def test_view_can_unpin_codons():
    view = CodonGraph('MIKEY').view()
    view.pin_codons({3: 'AAA'})
    view.unpin_codons([3])
    assert 3 not in view.pinned_codons


def test_pin_codons_validates_positions():
    view = CodonGraph('MIKEY').view()
    with pytest.raises(ValueError):
        view.pin_codons({999: 'ATG'})


def test_view_can_clear_pins():
    view = CodonGraph('MIKEY').view()
    view.pin_codons({3: 'AAA', 5: 'TAT'})
    view.clear_pins()
    assert view.pinned_codons == {}


def test_set_pinned_codons_replaces_existing_pins():
    view = CodonGraph('MIKEY').view()
    view.pin_codons({1: 'ATG'})
    view.set_pinned_codons({5: 'TAT'})
    assert view.pinned_codons == {5: ['TAT']}


def test_set_pinned_codons_compiles_once():
    view = CodonGraph('MIKEY').view()
    view.compile = MagicMock()
    view.set_pinned_codons({1: 'ATG'})
    view.compile.assert_called_once()


def test_set_pinned_codons_validates_positions():
    view = CodonGraph('MIKEY').view()
    with pytest.raises(ValueError):
        view.set_pinned_codons({999: 'ATG'})


def test_view_rejects_out_of_range_pin():
    view = CodonGraph('MIKEY').view()

    with pytest.raises(ValueError):
        view.pin_codons({0: 'ATG'})

    with pytest.raises(ValueError):
        view.pin_codons({6: 'ATG'})


def test_view_rejects_pin_outside_codon_restrictions():
    view = CodonGraph('MIKEY', codon_restrictions={3: ['AAA']}).view()
    with pytest.raises(ValueError):
        view.pin_codons({3: 'AAG'})


def test_view_getitem():
    view = CodonGraph('MF').view()
    assert view[0] == 'ATGTTT'
    assert view[1] == 'ATGTTC'


def test_get_works_for_very_large_sequences():
    view = CodonGraph('MIKEY' * 1000).view()
    _ = view[100]
    _ = view[1000000]
    _ = view[10**40]


def test_view_iter():
    view = CodonGraph('MIKEY').view()
    seqs = [*view]
    assert len(seqs) == len(set(seqs)) == 24


def test_view_contains():
    view = CodonGraph('MIKEY').view()
    for _ in range(100):
        seq = view.sample()
        assert seq in view
        assert seq + 'ATG' not in view


@pytest.mark.parametrize('aa_seq',
                         (
                                 'M',
                                 'MIKEY',
                                 'MILDRED',
                                 'ELEPHANT',
                                 'REGINALD',
                         )
                         )
def test_contains_passes_on_valid_sequences(aa_seq, standard_codon_table):
    view = CodonGraph(aa_seq).view()
    expected_all_seqs = helper_enumerate_sequences(aa_seq, standard_codon_table)
    for seq in expected_all_seqs:
        assert view.contains(seq)


def test_contains_fails_on_wrong_length_sequences():
    view = CodonGraph('MIKEY').view()
    assert not view.contains('')
    assert not view.contains('ATG')
    assert not view.contains('ATG' * 10)
    assert not view.contains('ATGA')  # not multiple of 3


@pytest.mark.parametrize(
    'aa_seq, invalid_seq',
    (
        ('M', 'ATT'),
        ('MIKEY', 'ATGATCAAAGAGTAA'),
    ),
)
def test_contains_fails_on_invalid_sequences(aa_seq, invalid_seq):
    view = CodonGraph(aa_seq).view()
    assert not view.contains(invalid_seq)


def test_contains_respects_pinning():
    view = CodonGraph('MS').view()
    assert view.contains('ATGTCT')
    assert view.contains('ATGTCC')

    view.pin_codons({2: 'TCT'})
    assert view.contains('ATGTCT')
    assert not view.contains('ATGTCC')

    view.clear_pins()
    assert view.contains('ATGTCT')
    assert view.contains('ATGTCC')


@pytest.fixture
def standard_codon_table():
    return {
        'F': ['TTT', 'TTC'],
        'L': ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'],
        'S': ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'],
        'Y': ['TAT', 'TAC'],
        'C': ['TGT', 'TGC'],
        'W': ['TGG'],
        'P': ['CCT', 'CCC', 'CCA', 'CCG'],
        'H': ['CAT', 'CAC'],
        'Q': ['CAA', 'CAG'],
        'R': ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
        'I': ['ATT', 'ATC', 'ATA'],
        'M': ['ATG'],
        'T': ['ACT', 'ACC', 'ACA', 'ACG'],
        'N': ['AAT', 'AAC'],
        'K': ['AAA', 'AAG'],
        'V': ['GTT', 'GTC', 'GTA', 'GTG'],
        'A': ['GCT', 'GCC', 'GCA', 'GCG'],
        'D': ['GAT', 'GAC'],
        'E': ['GAA', 'GAG'],
        'G': ['GGT', 'GGC', 'GGA', 'GGG']
    }


def helper_enumerate_sequences(aa_seq, aa_to_codons):
    codon_choices = [aa_to_codons[aa] for aa in aa_seq]
    seqs = [''.join(choices) for choices in product(*codon_choices)]
    return seqs


@pytest.mark.parametrize('aa_seq',
                         (
                                 'MIKEY',
                                 'MIKEY',
                                 'M' * 1000,
                                 'SSSSSS',
                                 'M',
                                 'MILDRED',
                                 'ELEPHANT',
                                 'REGINALD',
                         )
                         )
def test_n_valid_sequences_no_restrictions(aa_seq, standard_codon_table):
    view = CodonGraph(aa_seq).view()
    expected_n_all_seqs = len(helper_enumerate_sequences(aa_seq, standard_codon_table))
    assert view.n_valid_sequences == expected_n_all_seqs


def test_n_valid_sequences_fixed_codon(standard_codon_table):
    aa_seq = 'MIKEY'
    codon_restrictions = {2: 'ATC'}
    view = CodonGraph(aa_seq, codon_restrictions=codon_restrictions).view()

    sequences_all = helper_enumerate_sequences(aa_seq, standard_codon_table)
    sequences_restricted = [s for s in sequences_all if s[3:6] == 'ATC']

    assert len(sequences_restricted) != len(sequences_all)
    assert len(sequences_restricted) == view.n_valid_sequences


def test_n_valid_sequences_pinning_and_unpinning(standard_codon_table):
    aa_seq = 'MIKEY'
    view = CodonGraph(aa_seq).view()

    sequences_all = helper_enumerate_sequences(aa_seq, standard_codon_table)
    assert len(sequences_all) == view.n_valid_sequences

    codon_restrictions = {2: 'ATC'}
    sequences_restricted = [s for s in sequences_all if s[3:6] == 'ATC']
    assert len(sequences_restricted) != len(sequences_all)

    view.pin_codons(codon_restrictions)
    assert len(sequences_restricted) == view.n_valid_sequences

    view.clear_pins()
    assert len(sequences_all) == view.n_valid_sequences


@pytest.mark.parametrize('aa_seq',
                         (
                                 'MIKEY',
                                 'MIKEY',
                                 'M' * 1000,
                                 'SSSSSS',
                                 'M',
                                 'MILDRED',
                                 'ELEPHANT',
                                 'REGINALD',
                         )
                         )
def test_enumerate_sequences(aa_seq, standard_codon_table):
    view = CodonGraph(aa_seq).view()

    generated_all_seqs = [*view.enumerate()]
    expected_all_seqs = helper_enumerate_sequences(aa_seq, standard_codon_table)

    assert view.n_valid_sequences == len(generated_all_seqs) == len(expected_all_seqs)
    assert len(generated_all_seqs) == len(expected_all_seqs) == len(set(expected_all_seqs))
    assert set(generated_all_seqs) == set(expected_all_seqs)


def test_enumerate_pinned_sequences(standard_codon_table):
    aa_seq = 'MIKEY'
    view = CodonGraph(aa_seq).view()

    generated_all_seqs = [*view.enumerate()]
    expected_all_seqs = helper_enumerate_sequences(aa_seq, standard_codon_table)

    assert 24 == view.n_valid_sequences == len(generated_all_seqs) == len(expected_all_seqs)
    assert 24 == len(generated_all_seqs) == len(expected_all_seqs) == len(set(expected_all_seqs))
    assert set(generated_all_seqs) == set(expected_all_seqs)

    view.pin_codons({2: 'ATC'})
    generated_pinned_seqs = [*view.enumerate()]
    assert 8 == view.n_valid_sequences
    assert 8 == len(generated_pinned_seqs)
    assert all(seq[3:6] == 'ATC' for seq in generated_pinned_seqs)

    view.clear_pins()
    generated_unpinned_seqs = [*view.enumerate()]

    assert 24 == view.n_valid_sequences == len(generated_unpinned_seqs) == len(expected_all_seqs)
    assert 24 == len(generated_unpinned_seqs) == len(expected_all_seqs) == len(set(expected_all_seqs))
    assert set(generated_unpinned_seqs) == set(expected_all_seqs)


def test_view_seed_consistent():
    samples_by_rep = []

    for _ in range(10):
        view = CodonGraph('MIKEY').view(seed=8675309)
        samples = [view.sample() for _ in range(100)]
        samples_by_rep.append(samples)

    assert all(samples == samples_by_rep[0] for samples in samples_by_rep)


def test_view_no_seed_not_consistent():
    samples_by_rep = []

    for _ in range(10):
        view = CodonGraph('MIKEY').view()
        samples = [view.sample() for _ in range(100)]
        samples_by_rep.append(samples)

    assert not all(samples == samples_by_rep[0] for samples in samples_by_rep)


def test_view_seed_survives_recompile():
    view = CodonGraph('MIKEY').view(seed=8675309)

    before = [view.sample() for _ in range(20)]

    view.pin_codons({2: 'ATC'})
    pinned = [view.sample() for _ in range(20)]

    view.clear_pins()
    after = [view.sample() for _ in range(20)]

    repeat = CodonGraph('MIKEY').view(seed=8675309)

    expected_before = [repeat.sample() for _ in range(20)]
    repeat.pin_codons({2: 'ATC'})
    expected_pinned = [repeat.sample() for _ in range(20)]
    repeat.clear_pins()
    expected_after = [repeat.sample() for _ in range(20)]

    assert before == expected_before
    assert pinned == expected_pinned
    assert after == expected_after


def test_view_seed_and_rng_cannot_both_be_provided():
    import random

    with pytest.raises(ValueError):
        CodonGraph('MIKEY').view(seed=123, rng=random.Random(123))


def test_view_rng_is_used():
    import random

    rng = random.Random(8675309)
    view = CodonGraph('MIKEY').view(rng=rng)

    samples = [view.sample() for _ in range(100)]

    control_rng = random.Random(8675309)
    control_view = CodonGraph('MIKEY').view(rng=control_rng)
    expected = [control_view.sample() for _ in range(100)]

    assert samples == expected


def test_view_copy_copies_pins_not_rng_state():
    view = CodonGraph('MIKEY').view(seed=8675309)
    view.pin_codons({2: 'ATC'})

    copied = view.copy()

    assert copied.pinned_codons == view.pinned_codons
    assert copied.n_valid_sequences == view.n_valid_sequences
    assert [*copied.enumerate()] == [*view.enumerate()]


def test_codon_graph_view_pickle_preserves_random_state():
    view = CodonGraph('MIKEY').view(seed=8675309)
    _ = [view.sample() for _ in range(100)]

    loaded = pickle.loads(pickle.dumps(view))

    assert [loaded.sample() for _ in range(100)] == [view.sample() for _ in range(100)]


def test_codon_graph_view_pickle_preserves_pins():
    view = CodonGraph('MIKEY').view(seed=8675309)
    view.pin_codons({2: 'ATC'})

    loaded = pickle.loads(pickle.dumps(view))

    assert loaded.pinned_codons == view.pinned_codons
    assert loaded.n_valid_sequences == view.n_valid_sequences
    assert [*loaded.enumerate()] == [*view.enumerate()]