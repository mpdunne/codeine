import pickle
import pytest
import random

from itertools import product
from unittest.mock import MagicMock

from codeine.translation.tables import TranslationTable
from codeine.graph.graph import CodonGraph
from codeine.graph.nodes import CodonNode
from codeine.graph.tracking import BannedSequenceTracker


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


def test_getitem_works_for_very_large_sequences():
    view = CodonGraph('MIKEY' * 1000).view()
    _ = view[100]
    _ = view[1000000]
    _ = view[10**40]


def test_getitem_respects_pins():
    view = CodonGraph('MIKEY').view()
    view.pin_codons({2: 'ATC'})

    assert all(view[i][3:6] == 'ATC' for i in range(view.n_valid_sequences))


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


def test_contains_respects_contexts():
    view = CodonGraph('M', context_l='AAA', context_r='CCC').view()
    assert 'ATG' in view


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


def test_view_copy_copies_banned_sequences():
    view = CodonGraph('MIKEY').view()
    view.set_banned_sequences(['AAA'])

    copied = view.copy()

    assert copied.banned_sequences == view.banned_sequences
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


def test_pinning_reuses_banned_tracker():
    view = CodonGraph('MIKEY').view()
    view.set_banned_sequences(['GAATTC'])

    tracker = view._banned_tracker
    view.pin_codons({1: ['ATG']})
    assert view._banned_tracker is tracker

    tracker = view._banned_tracker
    view.clear_pins()
    assert view._banned_tracker is tracker


def test_setting_banned_sequences_rebuilds_tracker():
    view = CodonGraph('MIKEY').view()
    view.set_banned_sequences(['GAATTC'])

    tracker = view._banned_tracker

    view.set_banned_sequences(['GAATTC'])
    assert view._banned_tracker is not tracker

    view.set_banned_sequences(['GAATTC', 'ccgatt'])
    assert view._banned_tracker is not tracker


def test_view_doesnt_compile_immediately():
    view = CodonGraph('MIKEY').view()
    assert view._requires_compile


def test_n_valid_sequences_compiles_view():
    view = CodonGraph('MIKEY').view()

    assert view.n_valid_sequences == 24
    assert not view._requires_compile


def test_set_pinned_codons_marks_view_for_compile():
    view = CodonGraph('MIKEY').view()
    _ = view.n_valid_sequences

    view.compile = MagicMock(wraps=view.compile)
    view.set_pinned_codons({1: 'ATG'})

    view.compile.assert_not_called()
    assert view._requires_compile


def test_pin_codons_marks_view_for_compile():
    view = CodonGraph('MIKEY').view()
    _ = view.n_valid_sequences

    view.pin_codons({2: 'ATC'})

    assert view._requires_compile
    assert view.n_valid_sequences == 8
    assert not view._requires_compile


def test_unpin_codons_marks_view_for_compile():
    view = CodonGraph('MIKEY').view()
    view.pin_codons({2: 'ATC'})
    _ = view.n_valid_sequences

    view.unpin_codons([2])

    assert view._requires_compile
    assert view.n_valid_sequences == 24


def test_clear_pins_marks_view_for_compile():
    view = CodonGraph('MIKEY').view()
    view.pin_codons({2: 'ATC'})
    _ = view.n_valid_sequences

    view.clear_pins()

    assert view._requires_compile
    assert view.n_valid_sequences == 24


def test_public_methods_compile_if_required():
    view = CodonGraph('MIKEY').view()
    view.pin_codons({2: 'ATC'})

    assert view.sample()[3:6] == 'ATC'
    assert view[0][3:6] == 'ATC'
    assert 'ATGATCAAAGAGTAT' in view
    assert not view._requires_compile


def test_copy_preserves_compile_state():
    view = CodonGraph('MIKEY').view()
    view.pin_codons({2: 'ATC'})
    _ = view.n_valid_sequences

    copied = view.copy()

    assert not copied._requires_compile
    assert copied._compiled is view._compiled
    assert copied.n_valid_sequences == view.n_valid_sequences
    assert [*copied.enumerate()] == [*view.enumerate()]


def test_copy_preserves_uncompiled_state():
    view = CodonGraph('MIKEY').view()
    view.pin_codons({2: 'ATC'})

    copied = view.copy()

    assert copied._requires_compile
    assert copied.pinned_codons == view.pinned_codons
    assert copied.n_valid_sequences == 8


SHORT_AA_SEQUENCES = (
    'M',
    'MIKEY',
    'MILDRED',
    'ELEPHANT',
    'REGINALD',
)

MEDIUM_AA_SEQUENCES = (
    'MIKEY' * 10,
    'MIKEY' * 20,
)

# TODO Improve behaviour for long sequences
LONG_AA_SEQUENCES = (
    'MIKEY' * 100,
    'MIKEY' * 250,
#    'MIKEY' * 500,
#    'MIKEY' * 1000,
)

CONTEXTS_L = ('', 'a', 'aa', 'aaa', 'aaggaaggaagg')
CONTEXTS_R = ('', 't', 'tt', 'ttt', 'ttccttccttcc')


def helper_ban_sequences_and_check_enumerate(
        aa_seq,
        banned_sequences,
        context_l='',
        context_r='',
):
    tt = TranslationTable()

    unconstrained_graph = CodonGraph(
        aa_seq,
        context_l=context_l,
        context_r=context_r,
        translation_table=tt,
    )
    unconstrained_view = unconstrained_graph.view()
    unconstrained_seqs = set(unconstrained_view)

    graph = CodonGraph(
        aa_seq,
        context_l=context_l,
        context_r=context_r,
        translation_table=tt,
    )
    view = graph.view()
    view.set_banned_sequences(banned_sequences=banned_sequences)
    observed_seqs = set(view)

    normalised_banned_sequences = [
        banned_sequence.upper()
        for banned_sequence in banned_sequences
    ]

    context_l = context_l.upper()
    context_r = context_r.upper()

    expected_removed_seqs = {
        seq for seq in unconstrained_seqs
        if any(
            banned_sequence in f'{context_l}{seq.upper()}{context_r}'
            for banned_sequence in normalised_banned_sequences
        )
    }

    expected_remaining_seqs = unconstrained_seqs - expected_removed_seqs
    observed_removed_seqs = unconstrained_seqs - observed_seqs

    assert observed_removed_seqs == expected_removed_seqs
    assert observed_seqs == expected_remaining_seqs
    assert view.n_valid_sequences == len(expected_remaining_seqs)

    for seq in expected_removed_seqs:
        assert seq not in view

    for seq in expected_remaining_seqs:
        assert seq in view


def helper_ban_sequences_and_check_sample(
        aa_seq,
        banned_sequences,
        context_l='',
        context_r='',
        n_samples=100,
):
    tt = TranslationTable()

    unconstrained_graph = CodonGraph(
        aa_seq,
        context_l=context_l,
        context_r=context_r,
        translation_table=tt,
    )
    unconstrained_view = unconstrained_graph.view()

    graph = CodonGraph(
        aa_seq,
        context_l=context_l,
        context_r=context_r,
        translation_table=tt,
    )
    view = graph.view()
    view.set_banned_sequences(banned_sequences=banned_sequences)

    normalised_banned_sequences = [
        banned_sequence.upper()
        for banned_sequence in banned_sequences
    ]

    assert view.n_valid_sequences >= 0

    if view.n_valid_sequences == 0:
        with pytest.raises(ValueError):
            view.sample()
        return

    for _ in range(n_samples):
        seq = view.sample()

        assert seq in unconstrained_view
        assert seq in view

        full_seq = f'{context_l.upper()}{seq.upper()}{context_r.upper()}'

        assert all(
            banned_sequence not in full_seq
            for banned_sequence in normalised_banned_sequences
        )


#TODO: improve it.
_ = '''def helper_check_sampling_matches_expected_distributions(
        aa_seq,
        banned_sequences,
        context_l='',
        context_r='',
        n_samples=50_000,
):
    tt = TranslationTable()

    unconstrained_view = CodonGraph(
        aa_seq,
        context_l=context_l,
        context_r=context_r,
        translation_table=tt,
    ).view(seed=8675309)

    banned_view = CodonGraph(
        aa_seq,
        context_l=context_l,
        context_r=context_r,
        translation_table=tt,
    ).view(seed=8675309)

    banned_view.set_banned_sequences(banned_sequences)

    banned_sequences = [seq.upper() for seq in banned_sequences]

    fancy_counts = Counter(
        banned_view.sample()
        for _ in range(n_samples)
    )

    rejection_counts = Counter()
    while sum(rejection_counts.values()) < n_samples:
        seq = unconstrained_view.sample()
        if all(banned not in seq.upper() for banned in banned_sequences):
            rejection_counts[seq] += 1

    assert set(fancy_counts) == set(rejection_counts) == set(banned_view)

    for seq in banned_view:
        fancy_p = fancy_counts[seq] / n_samples
        rejection_p = rejection_counts[seq] / n_samples

        assert abs(fancy_p - rejection_p) < 0.02
'''

def test_view_sampler_keys_include_tracker_state():
    view = CodonGraph('MIKEY').view()

    for key in view.samplers:
        node, state = key
        assert isinstance(node, CodonNode)
        assert state == frozenset()


def test_sample_respects_pins():
    view = CodonGraph('MIKEY').view(seed=8675309)
    view.pin_codons({2: 'ATT'})

    for _ in range(100):
        assert view.sample()[3:6] == 'ATT'


@pytest.mark.parametrize('aa_seq', SHORT_AA_SEQUENCES)
@pytest.mark.parametrize('context_l', CONTEXTS_L)
@pytest.mark.parametrize('context_r', CONTEXTS_R)
def test_view_sampling_works_without_banned_sequences(aa_seq, context_l, context_r):
    helper_ban_sequences_and_check_sample(aa_seq, [], context_l, context_r, n_samples=1000)


@pytest.mark.parametrize('aa_seq', SHORT_AA_SEQUENCES)
@pytest.mark.parametrize('context_l', CONTEXTS_L)
@pytest.mark.parametrize('context_r', CONTEXTS_R)
def test_view_enumerate_works_without_banned_sequences(aa_seq, context_l, context_r):
    helper_ban_sequences_and_check_enumerate(aa_seq, [], context_l, context_r)


# TODO investigate this.
_ = '''
@pytest.mark.parametrize('aa_seq', SHORT_AA_SEQUENCES)
@pytest.mark.parametrize('context_l', CONTEXTS_L)
@pytest.mark.parametrize('context_r', CONTEXTS_R)
def test_banned_sampling_matches_rejection_sampling(aa_seq, context_l, context_r):
    helper_check_sampling_matches_expected_distributions(aa_seq, [], context_l=context_l, context_r=context_r)'''


def test_banned_sequence_entirely_in_left_context_gives_empty_space():
    graph = CodonGraph('MIKEY', context_l='GAATTC')
    view = graph.view()
    view.set_banned_sequences(['GAATTC'])
    assert view.n_valid_sequences == 0

    graph = CodonGraph('MIKEY', context_l='GAATTC')
    view = graph.view()
    view.set_banned_sequences(['AATTC'])
    assert view.n_valid_sequences == 0

    graph = CodonGraph('MIKEY', context_l='GAATTC')
    view = graph.view()
    view.set_banned_sequences(['GAATT'])
    assert view.n_valid_sequences == 0


def test_banned_sequence_entirely_in_right_context_gives_empty_space():
    graph = CodonGraph('MIKEY', context_l='GAATTC')
    view = graph.view()
    view.set_banned_sequences(['GAATTC'])
    assert view.n_valid_sequences == 0

    graph = CodonGraph('MIKEY', context_l='GAATTC')
    view = graph.view()
    view.set_banned_sequences(['AATTC'])
    assert view.n_valid_sequences == 0

    graph = CodonGraph('MIKEY', context_r='GAATTC')
    view = graph.view()
    view.set_banned_sequences(['AATT'])
    assert view.n_valid_sequences == 0


def test_banned_sequence_spanning_left_context_and_first_codon_is_respected():
    view = CodonGraph('ELEPHANT', context_l='AAGGATGATG').view()
    view.set_banned_sequences(['AAGGATGATGAA'])

    for seq in view:
        assert 'AAGGATGATGAA' not in f'AAGGATGATG{seq}'


def test_banned_sequence_spanning_last_codon_and_right_context_is_respected():
    view = CodonGraph('ELEPHANT', context_r='AAGGATGATG').view()
    view.set_banned_sequences(['CGAAGGATGATG'])

    for seq in view:
        assert 'CGAAGGATGATG' not in f'{seq}AAGGATGATG'


def test_banned_sequence_spanning_both_contexts_is_respected():
    view = CodonGraph('ELEPHANT', context_l='TTAA', context_r='AAGG').view()
    view.set_banned_sequences(['AA' + 'GAGCTTGAGCCGCATGCCAATACG' + 'AA'])
    assert 'GAGCTTGAGCCGCATGCCAATACG' not in view


def helper_arbitrary_coding_sequence(aa_seq, translation_table):
    return ''.join(translation_table.aa_to_codons[aa][0] for aa in aa_seq)


def test_tracker_finds_ban_inside_left_context():
    graph = CodonGraph('REGINALD', context_l='aaggaaggaagg')
    banned_seqs = (
        'ATG',
        'TAAAAG',
        'AAGGAA',
        'ATTAAGG',
        'GAATAC',
    )
    tracker = BannedSequenceTracker(graph, banned_seqs)
    assert tracker.paths


@pytest.mark.parametrize('aa_seq', SHORT_AA_SEQUENCES)
@pytest.mark.parametrize('context_l', CONTEXTS_L)
@pytest.mark.parametrize('context_r', CONTEXTS_R)
def test_regression_banned_sequences_short_aa_sequence(aa_seq, context_l, context_r):
    banned_seqs = (
        'ATG',
        'TAAAAG',
        'AAGGAA',
        'ATTAAGG',
        'GAATAC',
    )
    helper_ban_sequences_and_check_enumerate(aa_seq, banned_seqs,
                                             context_l=context_l, context_r=context_r)


@pytest.mark.parametrize('aa_seq', SHORT_AA_SEQUENCES)
@pytest.mark.parametrize('context_l', CONTEXTS_L)
@pytest.mark.parametrize('context_r', CONTEXTS_R)
def test_regression_banned_sequences_short_aa_sequence_overlapping_banned_sequences(aa_seq, context_l, context_r):
    banned_seqs = (
        'TAAAAG',
        'AAAGGA',
        'AAGGAA',
        'GGAATA',
    )
    helper_ban_sequences_and_check_enumerate(aa_seq, banned_seqs,
                                             context_l=context_l, context_r=context_r)


@pytest.mark.parametrize('aa_seq', SHORT_AA_SEQUENCES + MEDIUM_AA_SEQUENCES + LONG_AA_SEQUENCES)
@pytest.mark.parametrize('context_l', CONTEXTS_L)
@pytest.mark.parametrize('context_r', CONTEXTS_R)
def test_regression_banned_whole_sequences(aa_seq, context_l, context_r):
    tt = TranslationTable()
    cds = helper_arbitrary_coding_sequence(aa_seq, tt)

    unconstrained_graph = CodonGraph(
        aa_seq,
        context_l=context_l,
        context_r=context_r,
        translation_table=tt,
    )
    unconstrained_view = unconstrained_graph.view()

    assert cds in unconstrained_view
    unconstrained_n_sequences = unconstrained_view.n_valid_sequences
    assert unconstrained_n_sequences >= 0

    graph = CodonGraph(
        aa_seq,
        context_l=context_l,
        context_r=context_r,
        translation_table=tt,
    )

    view = graph.view()
    view.set_banned_sequences([cds])

    assert cds not in view
    assert view.n_valid_sequences == unconstrained_n_sequences - 1

    # Same thing but with multiple sequences.
    seqs = [unconstrained_view[i] for i in range(min(unconstrained_n_sequences, 5))]

    graph = CodonGraph(
        aa_seq,
        context_l=context_l,
        context_r=context_r,
        translation_table=tt,
    )

    view = graph.view()
    view.set_banned_sequences(seqs)

    for seq in seqs:
        assert seq not in view

    assert view.n_valid_sequences == unconstrained_n_sequences - len(seqs)


@pytest.mark.parametrize('aa_seq', SHORT_AA_SEQUENCES + MEDIUM_AA_SEQUENCES + LONG_AA_SEQUENCES)
@pytest.mark.parametrize('context_l', CONTEXTS_L)
@pytest.mark.parametrize('context_r', CONTEXTS_R)
def test_regression_banned_sequences_that_arent_present_anyway(aa_seq, context_l, context_r):
    tt = TranslationTable()

    unconstrained_graph = CodonGraph(
        aa_seq,
        context_l=context_l,
        context_r=context_r,
        translation_table=tt,
    )

    unconstrained_view = unconstrained_graph.view()

    graph = CodonGraph(
        aa_seq,
        context_l=context_l,
        context_r=context_r,
        translation_table=tt,
    )

    view = graph.view()
    banned_sequences = ['CCCCCCCCCCCC', 'GGGGGGGGGGGG', 'TTTTTTTTTTTT']
    view.set_banned_sequences(banned_sequences)

    assert view.n_valid_sequences == unconstrained_view.n_valid_sequences


@pytest.mark.parametrize('aa_seq', MEDIUM_AA_SEQUENCES + LONG_AA_SEQUENCES)
@pytest.mark.parametrize('context_l', CONTEXTS_L)
@pytest.mark.parametrize('context_r', CONTEXTS_R)
def test_regression_banned_sequences_long_aa_sequence(aa_seq, context_l, context_r):
    banned_seqs = (
        'TAAAAG',
        'AAAGGA',
        'AAGGAA',
        'GGAATA',
    )

    helper_ban_sequences_and_check_sample(aa_seq, banned_seqs,
                                          context_l, context_r, n_samples=1000)


@pytest.mark.parametrize('aa_seq', MEDIUM_AA_SEQUENCES + LONG_AA_SEQUENCES)
@pytest.mark.parametrize('context_l', CONTEXTS_L)
@pytest.mark.parametrize('context_r', CONTEXTS_R)
def test_regression_banned_sequences_long_aa_sequence_overlapping_banned_sequences(aa_seq, context_l, context_r):
    tt = TranslationTable()

    seq = helper_arbitrary_coding_sequence(aa_seq, tt)

    banned_sequences = [
        seq[0:12],
        seq[3:15],
        seq[6:18],
        seq[9:21],
    ]

    helper_ban_sequences_and_check_sample(aa_seq, banned_sequences,
                                          context_l, context_r, n_samples=1000)


@pytest.mark.parametrize('aa_seq', MEDIUM_AA_SEQUENCES + LONG_AA_SEQUENCES)
@pytest.mark.parametrize('context_l', CONTEXTS_L)
@pytest.mark.parametrize('context_r', CONTEXTS_R)
def test_regression_banned_sequences_long_aa_sequence_nested_banned_sequences(aa_seq, context_l, context_r):
    tt = TranslationTable()

    seq = helper_arbitrary_coding_sequence(aa_seq, tt)

    banned_sequences = (
        seq[0:12],
        seq[3:12],
        seq[5:10],
    )

    helper_ban_sequences_and_check_sample(aa_seq, banned_sequences,
                                                 context_l, context_r, n_samples=1000)


@pytest.mark.parametrize('aa_seq', MEDIUM_AA_SEQUENCES + LONG_AA_SEQUENCES)
@pytest.mark.parametrize('context_l', CONTEXTS_L)
@pytest.mark.parametrize('context_r', CONTEXTS_R)
def test_regression_banned_sequences_long_aa_sequence_long_banned_sequences(aa_seq, context_l, context_r):

    tt = TranslationTable()

    unconstrained_graph = CodonGraph(
        aa_seq,
        context_l=context_l,
        context_r=context_r,
        translation_table=tt,
    )

    unconstrained_view = unconstrained_graph.view()

    # Jenny?????
    rng = random.Random(8675309)

    # Grab 50 random 50nt-long sequences from here.
    # (I personally think that's "long", don't know about you!)
    banned_seqs = []
    for _ in range(50):
        ix = rng.randrange(unconstrained_view.n_valid_sequences)
        seq = unconstrained_view[ix]
        start = rng.randrange(len(seq) - 49)
        banned_seqs.append(seq[start:start + 50])

    helper_ban_sequences_and_check_sample(
        aa_seq,
        banned_seqs,
        context_l=context_l,
        context_r=context_r,
    )


@pytest.mark.parametrize('aa_seq', SHORT_AA_SEQUENCES)
@pytest.mark.parametrize('context_l', CONTEXTS_L)
@pytest.mark.parametrize('context_r', CONTEXTS_R)
def test_regression_banned_sequences_short_aa_sequence_many_banned_sequences(aa_seq, context_l, context_r):
    banned_seqs = (
        'ATG',
        'TAA',
        'AAG',
        'GAA',
        'TAC',
        'ATT',
        'CTG',
        'GGT',
        'GGC',
        'AAA',
        'AAC',
        'GAG',
    )

    helper_ban_sequences_and_check_enumerate(
        aa_seq,
        banned_seqs,
        context_l=context_l,
        context_r=context_r,
    )


@pytest.mark.parametrize('aa_seq', MEDIUM_AA_SEQUENCES + LONG_AA_SEQUENCES)
@pytest.mark.parametrize('context_l', CONTEXTS_L)
@pytest.mark.parametrize('context_r', CONTEXTS_R)
def test_regression_banned_sequences_long_aa_sequence_many_banned_sequences(aa_seq, context_l, context_r):
    tt = TranslationTable()

    unconstrained_graph = CodonGraph(
        aa_seq,
        context_l=context_l,
        context_r=context_r,
        translation_table=tt,
    )

    unconstrained_view = unconstrained_graph.view()

    rng = random.Random(8675309)

    banned_seqs = []
    for _ in range(500):
        ix = rng.randrange(unconstrained_view.n_valid_sequences)
        seq = unconstrained_view[ix]

        start = rng.randrange(len(seq) - 11)
        banned_seqs.append(seq[start:start + 12])

    helper_ban_sequences_and_check_sample(
        aa_seq,
        banned_seqs,
        context_l=context_l,
        context_r=context_r,
        n_samples=1000,
    )
