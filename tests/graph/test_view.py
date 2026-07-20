import pickle
import pytest
import random

from collections import Counter
from itertools import product
from unittest.mock import MagicMock
from scipy.stats import chisquare, chi2_contingency

from codeine.constraints.base import Constraint, DEAD_STATE, SAFE_STATE
from codeine.constraints.gc import GCConstraint
from codeine.constraints.banned import BannedSequenceConstraint
from codeine.translation.tables import TranslationTable
from codeine.translation.weights import CodonWeights
from codeine.graph.base import CodonGraph
from codeine.graph.view import CodonGraphView, COMPILED, COMPILE_SHALLOW, COMPILE_DEEP, COMPILE_EXTEND

from tests.data import NORMAL_PROTEINS


def test_view_exposes_translation_table():
    graph = CodonGraph('MIKEY')
    view = graph.view()
    assert view.translation_table is graph.tt


def test_view_exposes_codon_restrictions():
    graph = CodonGraph('MIKEY', codon_restrictions={2: 'ATC'})
    view = graph.view()
    assert view.codon_restrictions == graph.codon_restrictions


def test_view_exposes_contexts():
    graph = CodonGraph('MIKEY', context_l='AAA', context_r='CCC')
    view = graph.view()
    assert view.context_l == 'AAA'
    assert view.context_r == 'CCC'


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
    seqs = [*view.enumerate()]

    for i in range(len(seqs)):
        assert view[i] == seqs[i]


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


class RejectAllConstraint(Constraint):

    @property
    def initial_state(self):
        return 0

    def advance(self, state, pos, choice):
        return DEAD_STATE

    def link(self, graph):
        pass


def test_constraint_can_reject_all_sequences():
    view = CodonGraph('MIKEY', context_l='aaa', context_r='ttt').view()

    view.set_constraints([RejectAllConstraint()])

    assert view.n_valid_sequences == 0
    assert [*view.enumerate()] == []

    with pytest.raises(ValueError):
        view.sample()


def test_clear_constraints_restores_default_behaviour():
    view = CodonGraph('MIKEY', context_l='aaa', context_r='ttt').view()

    expected_sequences = [*view.enumerate()]
    expected_count = view.n_valid_sequences

    view.set_constraints([RejectAllConstraint()])
    assert view.n_valid_sequences == 0

    view.clear_constraints()

    assert view.n_valid_sequences == expected_count
    assert [*view.enumerate()] == expected_sequences


def test_constraints_are_copied_with_view():
    view = CodonGraph('MIKEY', context_l='aaa', context_r='ttt').view()
    view.set_constraints([RejectAllConstraint()])

    copied = view.copy()

    assert copied.n_valid_sequences == 0
    assert [*copied.enumerate()] == []


class RejectChoiceConstraint(Constraint):

    def __init__(self, rejected_choice):
        self.rejected_choice = rejected_choice

    @property
    def initial_state(self):
        return 0

    def advance(self, state, pos, choice):
        if state < 0:
            return state

        if choice == self.rejected_choice:
            return DEAD_STATE

        return 0

    def link(self, graph):
        pass


class SafeAfterFirstConstraint(Constraint):

    @property
    def initial_state(self):
        return 0

    def advance(self, state, pos, choice):
        if state < 0:
            return state

        if pos >= 1:
            return SAFE_STATE

        return 0

    def link(self, graph):
        pass


def test_view_applies_multiple_constraints():
    graph = CodonGraph('MIKEY', context_l='aaa', context_r='ttt')
    all_sequences = set(graph.view().enumerate())

    view = graph.view()
    view.set_constraints([RejectChoiceConstraint('ATA'), RejectChoiceConstraint('AAG')])

    expected = {sequence for sequence in all_sequences
                if sequence[3:6] != 'ATA' and sequence[6:9] != 'AAG'}

    assert set(view.enumerate()) == expected


def test_view_rejects_choice_when_any_constraint_is_dead():
    graph = CodonGraph('MIKEY')
    all_sequences = set(graph.view().enumerate())

    view = graph.view()
    view.set_constraints([SafeAfterFirstConstraint(), RejectChoiceConstraint('GAG')]),

    expected = {sequence for sequence in all_sequences if sequence[9:12] != 'GAG'}

    assert set(view.enumerate()) == expected


def test_safe_constraint_does_not_disable_other_constraints():
    graph = CodonGraph('MIKEY')
    all_sequences = set(graph.view().enumerate())

    view = graph.view()
    view.set_constraints([SafeAfterFirstConstraint(), RejectChoiceConstraint('TAC')])

    expected = {sequence for sequence in all_sequences if sequence[12:15] != 'TAC'}

    assert set(view.enumerate()) == expected


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


def test_enumerate_range_matches_full_enumeration_slice():
    view = CodonGraph('MIKEY').view()
    all_seqs = [*view.enumerate()]

    assert [*view.enumerate_range(0, 0)] == []
    assert [*view.enumerate_range(5, 5)] == []
    assert [*view.enumerate_range(0, 5)] == all_seqs[0:5]
    assert [*view.enumerate_range(5, 12)] == all_seqs[5:12]
    assert [*view.enumerate_range(12, 24)] == all_seqs[12:24]
    assert [*view.enumerate_range(0)] == all_seqs


def test_sequences_at_matches_full_enumeration_slice():
    view = CodonGraph('MIKEY').view()
    all_seqs = [*view.enumerate()]

    assert view[5:5] == []
    assert view[:0] == all_seqs[:0]
    assert view[:5] == all_seqs[:5]
    assert view[5:12] == all_seqs[5:12]
    assert view[12:] == all_seqs[12:]
    assert view[:] == all_seqs


def test_sequences_at_with_step_matches_expected_slice():
    view = CodonGraph('MIKEY').view()
    all_seqs = [*view.enumerate()]

    assert view[::2] == all_seqs[::2]
    assert view[1::3] == all_seqs[1::3]
    assert view[20:2:-4] == all_seqs[20:2:-4]


def test_sequence_at_matches_single_item_slice():
    view = CodonGraph('MIKEY').view()

    for index in range(view.n_valid_sequences):
        assert view.sequence_at(index) == view[index:index + 1][0]


def test_large_sequence_at_matches_large_slice():
    view = CodonGraph('MIKEY' * 1000).view()

    index = 10**69

    assert view.sequence_at(index) == view[index:index + 1][0]


def test_large_enumerate_range_matches_repeated_sequence_at():
    view = CodonGraph('MIKEY' * 1000).view()

    start = 10**69
    stop = start + 10

    observed = [*view.enumerate_range(start, stop)]
    expected = [view.sequence_at(index) for index in range(start, stop)]

    assert observed == expected


def test_large_slice_matches_large_enumerate_range():
    view = CodonGraph('MIKEY' * 1000).view()

    start = 10**69
    stop = start + 10

    assert view[start:stop] == [*view.enumerate_range(start, stop)]


def test_large_stepped_slice_matches_repeated_sequence_at():
    view = CodonGraph('MIKEY' * 1000).view()

    start = 10**69
    stop = start + 42
    step = 7

    observed = view[start:stop:step]
    expected = [view.sequence_at(index) for index in range(start, stop, step)]

    assert observed == expected


def test_enumerate_range_raise_if_bounds_are_bad():
    view = CodonGraph('MIKEY').view()
    n = view.n_valid_sequences

    with pytest.raises(IndexError):
        _ = [*view.enumerate_range(-1, 1)]

    with pytest.raises(IndexError):
        _ = [*view.enumerate_range(2, 1)]

    with pytest.raises(IndexError):
        _ = [*view.enumerate_range(0, n + 1)]


def test_large_enumerate_range_respects_pins_and_constraints():
    view = CodonGraph('MIKEY' * 100).view()
    view.pin_codons({2: 'ATC'})
    view.set_constraints([RejectChoiceConstraint('TAC')])

    start = 10**69
    stop = start + 10

    observed = [*view.enumerate_range(start, stop)]
    expected = [view.sequence_at(index) for index in range(start, stop)]

    assert observed == expected
    assert all(seq[3:6] == 'ATC' for seq in observed)
    assert all('TAC' not in seq for seq in observed)


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


def test_view_copy_copies_pins():
    view = CodonGraph('MIKEY').view(seed=8675309)
    view.pin_codons({2: 'ATC'})

    copied = view.copy()

    assert copied.pinned_codons == view.pinned_codons
    assert copied.n_valid_sequences == view.n_valid_sequences
    assert [*copied.enumerate()] == [*view.enumerate()]


def test_view_copy_copies_constraints():
    view = CodonGraph('MIKEY').view()
    view.set_constraints([RejectChoiceConstraint('ATA')])

    copied = view.copy()

    assert copied.constraints == view.constraints
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


def test_view_doesnt_compile_immediately():
    view = CodonGraph('MIKEY').view()
    assert view._compile_status


def test_n_valid_sequences_compiles_view():
    view = CodonGraph('MIKEY').view()

    assert view.n_valid_sequences == 24
    assert not view._compile_status


def test_set_pinned_codons_marks_view_for_compile():
    view = CodonGraph('MIKEY').view()
    _ = view.n_valid_sequences

    view.compile = MagicMock(wraps=view.compile)
    view.set_pinned_codons({1: 'ATG'})

    view.compile.assert_not_called()
    assert view._compile_status


def test_pin_codons_marks_view_for_compile():
    view = CodonGraph('MIKEY').view()
    _ = view.n_valid_sequences

    view.pin_codons({2: 'ATC'})

    assert view._compile_status
    assert view.n_valid_sequences == 8
    assert not view._compile_status


def test_unpin_codons_marks_view_for_compile():
    view = CodonGraph('MIKEY').view()
    view.pin_codons({2: 'ATC'})
    _ = view.n_valid_sequences

    view.unpin_codons([2])

    assert view._compile_status
    assert view.n_valid_sequences == 24


def test_clear_pins_marks_view_for_compile():
    view = CodonGraph('MIKEY').view()
    view.pin_codons({2: 'ATC'})
    _ = view.n_valid_sequences

    view.clear_pins()

    assert view._compile_status
    assert view.n_valid_sequences == 24


def test_public_methods_compile_if_required():
    view = CodonGraph('MIKEY').view()
    view.pin_codons({2: 'ATC'})
    assert view._compile_status

    view = CodonGraph('MIKEY').view()
    assert view._compile_status
    _ = view.sample()
    assert not view._compile_status

    view = CodonGraph('MIKEY').view()
    assert view._compile_status
    _ = view[0]
    assert not view._compile_status

    view = CodonGraph('MIKEY').view()
    assert view._compile_status
    assert 'ATGATCAAAGAGTAT' in view
    assert not view._compile_status


def test_changing_copied_view_pins_leaves_original_untouched():
    view = CodonGraph('MIKEY').view()
    copied = view.copy()

    copied.pin_codons({1: 'ATG'})
    assert view.pinned_codons == {}
    assert copied.pinned_codons == {1: ['ATG']}

    view.pin_codons({2: 'ATT'})
    assert view.pinned_codons == {2: ['ATT']}
    assert copied.pinned_codons == {1: ['ATG']}


def test_changing_copied_view_changing_constraints_leaves_original_untouched():
    view = CodonGraph('MIKEY').view()
    view.set_constraints([RejectChoiceConstraint('AAA')])
    copied = view.copy()

    copied.set_constraints([RejectChoiceConstraint('TTT')])
    assert len(view.constraints) == 1
    assert view.constraints[0].rejected_choice == 'AAA'
    assert len(copied.constraints) == 1
    assert copied.constraints[0].rejected_choice == 'TTT'

    view.clear_constraints()
    assert len(view.constraints) == 0
    assert len(copied.constraints) == 1
    assert copied.constraints[0].rejected_choice == 'TTT'


def test_copied_view_recompiles_independently_after_change():
    view = CodonGraph('M').view()
    copied = view.copy()

    assert view._compile_status
    assert copied._compile_status

    assert view.n_valid_sequences == copied.n_valid_sequences
    assert not view._compile_status
    assert not copied._compile_status

    copied.set_constraints([RejectChoiceConstraint('ATG')])
    assert not view._compile_status
    assert copied._compile_status

    assert copied.n_valid_sequences == 0
    assert view.n_valid_sequences > 0
    assert not view._compile_status
    assert not copied._compile_status


def test_copy_preserves_compile_state():
    view = CodonGraph('MIKEY').view()
    view.pin_codons({2: 'ATC'})
    _ = view.n_valid_sequences

    copied = view.copy()

    assert not copied._compile_status
    assert copied._compiled is view._compiled
    assert copied.n_valid_sequences == view.n_valid_sequences
    assert [*copied.enumerate()] == [*view.enumerate()]


def test_copy_preserves_uncompiled_state():
    view = CodonGraph('MIKEY').view()
    view.pin_codons({2: 'ATC'})

    copied = view.copy()

    assert copied._compile_status
    assert copied.pinned_codons == view.pinned_codons
    assert copied.n_valid_sequences == 8


def test_copied_view_copies_rng_state():
    view = CodonGraph('MIKEY').view(seed=123)
    copied = view.copy()

    assert copied._rng is not view._rng
    assert copied.sample() == view.sample()
    assert copied.sample() == view.sample()


def test_view_defaults_to_dna_when_only_weights_are_given():
    weights = CodonWeights.ecoli()

    view = CodonGraph('MIKEY').view(weights=weights)

    assert view.translation_table.rna is False
    assert view.codon_weights is weights


def test_view_preserves_matching_weights():
    tt = TranslationTable()
    weights = CodonWeights.ecoli()

    view = CodonGraph('MIKEY', translation_table=tt).view(weights=weights)

    assert view.translation_table is tt
    assert view.codon_weights is weights


def test_view_converts_dna_weights_for_rna_table():
    tt = TranslationTable(rna=True)
    weights = CodonWeights.ecoli()

    view = CodonGraph('MIKEY', translation_table=tt).view(weights=weights)

    assert view.translation_table is tt
    assert view.codon_weights is not weights
    assert 'AUG' in view.codon_weights.weights
    assert 'ATG' not in view.codon_weights.weights
    assert 'ATG' in weights.weights
    assert 'AUG' not in weights.weights


def test_view_converts_rna_weights_for_dna_table():
    rna_table = TranslationTable(rna=True)
    weights = CodonWeights.uniform(table=rna_table)

    view = CodonGraph('MIKEY', translation_table=TranslationTable()).view(weights=weights)

    assert view.translation_table.rna is False
    assert view.codon_weights is not weights
    assert 'ATG' in view.codon_weights.weights
    assert 'AUG' not in view.codon_weights.weights
    assert 'AUG' in weights.weights
    assert 'ATG' not in weights.weights


def test_view_rejects_incompatible_codon_weights():
    weights = CodonWeights.uniform()

    with pytest.raises(ValueError):
        CodonGraph('MIKEY', translation_table=TranslationTable(table_id=2)).view(weights=weights)


def test_view_defaults_to_uniform_weights():
    graph = CodonGraph('MIKEY')
    view = graph.view()

    assert view.codon_weights.weights == CodonWeights.uniform(table=graph.tt).weights


def test_view_uses_provided_weights():
    graph = CodonGraph('MIKEY')
    weights = CodonWeights.ecoli()

    view = graph.view(weights=weights)

    assert view.codon_weights is weights


def test_view_converts_weights_to_graph_molecule_type():
    graph = CodonGraph('MIKEY', translation_table=TranslationTable(rna=True))
    weights = CodonWeights.ecoli()

    view = graph.view(weights=weights)

    assert view.codon_weights is not weights
    assert 'AUG' in view.codon_weights.weights
    assert 'ATG' not in view.codon_weights.weights


def test_set_weights_updates_weights_and_requires_compile():
    graph = CodonGraph('MIKEY')
    view = graph.view()

    view.compile()
    assert not view._compile_status

    weights = CodonWeights.ecoli()
    view.set_weights(weights)

    assert view.codon_weights is weights
    assert view._compile_status


def test_clear_weights_restores_uniform_weights():
    graph = CodonGraph('MIKEY')
    view = graph.view(weights=CodonWeights.ecoli())

    view.clear_weights()

    assert view.codon_weights.weights == CodonWeights.uniform(
        table=graph.tt,
    ).weights


def test_copy_preserves_weights():
    graph = CodonGraph('MIKEY')
    weights = CodonWeights.ecoli()
    view = graph.view(weights=weights)

    copied = view.copy()

    assert copied.codon_weights is weights


def test_sample_respects_pins():
    view = CodonGraph('MIKEY').view(seed=8675309)
    view.pin_codons({2: 'ATT'})

    for _ in range(100):
        assert view.sample()[3:6] == 'ATT'


def test_sample_many_returns_correct_number():
    view = CodonGraph('MIKEY').view(seed=8675309)
    seqs = view.sample(n=10)

    assert len(seqs) == 10
    assert all(isinstance(seq, str) for seq in seqs)


def test_sample_many_zero():
    view = CodonGraph('MIKEY').view()
    assert view.sample(n=0) == []


def test_sample_many_negative_n_raises():
    view = CodonGraph('MIKEY').view()
    with pytest.raises(ValueError):
        view.sample(n=-1)


def test_sample_many_sequences_are_valid():
    view = CodonGraph('MIKEY').view()
    seqs = view.sample(n=100)
    assert all(seq in view for seq in seqs)


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

    bsc = BannedSequenceConstraint(banned_sequences)
    view.set_constraints([bsc])

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

    graph = CodonGraph(
        aa_seq,
        context_l=context_l,
        context_r=context_r,
        translation_table=tt,
    )
    view = graph.view()
    bsc = BannedSequenceConstraint(banned_sequences)
    view.set_constraints([bsc])

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

        assert tt.translate(seq) == aa_seq

        full_seq = f'{context_l.upper()}{seq.upper()}{context_r.upper()}'

        assert all(
            banned_sequence not in full_seq
            for banned_sequence in normalised_banned_sequences
        )


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


def test_banned_sequences_can_be_set_and_are_normalised():
    graph = CodonGraph('MIKEY')
    view = CodonGraphView(graph, constraints=[BannedSequenceConstraint(['aaa', 'AAA', 'ttt'])])
    assert len(view.constraints) == 1
    assert isinstance(view.constraints[0], BannedSequenceConstraint)
    assert set(view.constraints[0].banned_sequences) == {'AAA', 'TTT'}

    view = CodonGraphView(graph, constraints=[BannedSequenceConstraint(['ccc', 'CCC', 'ggg'])])
    assert len(view.constraints) == 1
    assert isinstance(view.constraints[0], BannedSequenceConstraint)
    assert set(view.constraints[0].banned_sequences) == {'CCC', 'GGG'}


def test_clear_banned_sequences_removes_bans_and_marks_stale():
    graph = CodonGraph('MIKEY')
    bsc = BannedSequenceConstraint(['aaa', 'AAA', 'ttt'])
    view = CodonGraphView(graph, constraints=[bsc])
    assert view.constraints == (bsc,)
    assert view._compile_status

    view.compile()
    assert not view._compile_status

    view.clear_constraints()
    assert view.constraints == ()
    assert view._compile_status


def test_view_exposes_graph_properties():
    graph = CodonGraph('MIKEY')
    view = graph.view()

    assert view.aa_seq == graph.aa_seq
    assert view.translation_table is graph.tt
    assert view.codon_restrictions is graph.codon_restrictions
    assert view.context_l == graph.context_l
    assert view.context_r == graph.context_r


@pytest.mark.parametrize(
    'banned_sequence',
    (
        'GAATTC',
        'AATTC',
        'GAATT',
    ),
)
def test_regression_banned_sequence_entirely_in_left_context_gives_empty_space(banned_sequence):
    graph = CodonGraph('MIKEY', context_l='GAATTC')
    view = graph.view()
    bsc = BannedSequenceConstraint([banned_sequence])
    view.set_constraints([bsc])
    assert view.n_valid_sequences == 0


@pytest.mark.parametrize(
    'banned_sequence',
    (
        'GAATTC',
        'AATTC',
        'GAATT',
    ),
)
def test_regression_banned_sequence_entirely_in_right_context_gives_empty_space(banned_sequence):
    graph = CodonGraph('MIKEY', context_r='GAATTC')
    view = graph.view()
    bsc = BannedSequenceConstraint([banned_sequence])
    view.set_constraints([bsc])
    assert view.n_valid_sequences == 0


def test_regression_banned_sequence_spanning_left_context_and_first_codon_is_respected():
    view = CodonGraph('ELEPHANT', context_l='AAGGATGATG').view()
    bsc = BannedSequenceConstraint(['AAGGATGATGAA'])
    view.set_constraints([bsc])

    for seq in view:
        assert 'AAGGATGATGAA' not in f'AAGGATGATG{seq}'


def test_regression_banned_sequence_spanning_last_codon_and_right_context_is_respected():
    view = CodonGraph('ELEPHANT', context_r='AAGGATGATG').view()
    bsc = BannedSequenceConstraint(['CGAAGGATGATG'])
    view.set_constraints([bsc])

    for seq in view:
        assert 'CGAAGGATGATG' not in f'{seq}AAGGATGATG'


def test_regression_banned_sequence_spanning_both_contexts_is_respected():
    view = CodonGraph('ELEPHANT', context_l='TTAA', context_r='AAGG').view()
    bsc = BannedSequenceConstraint(['AA' + 'GAGCTTGAGCCGCATGCCAATACG' + 'AA'])
    view.set_constraints([bsc])
    assert 'GAGCTTGAGCCGCATGCCAATACG' not in view


def helper_arbitrary_coding_sequence(aa_seq, translation_table):
    return ''.join(translation_table.aa_to_codons[aa][0] for aa in aa_seq)


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

    bsc = BannedSequenceConstraint([cds])
    view.set_constraints([bsc])

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
    bsc = BannedSequenceConstraint(seqs)
    view.set_constraints([bsc])

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
    bsc = BannedSequenceConstraint(banned_sequences)
    view.set_constraints([bsc])

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
    n_valid_sequences = unconstrained_view.n_valid_sequences
    for _ in range(50):
        ix = rng.randrange(n_valid_sequences)
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


def helper_chi_square_codon_test(observed_counts, expected_counts):
    """
    Check that a set of sampled items looks like it's drawn from a given distirbution.
    """
    codons = sorted({*expected_counts, *observed_counts})

    obs = [observed_counts.get(codon, 0) for codon in codons]
    exp = [expected_counts.get(codon, 0) for codon in codons]

    return chisquare(obs, exp)


def helper_codon_counts_by_position(seqs):
    n_codons = len(seqs[0]) // 3
    counts = {pos: Counter() for pos in range(1, n_codons + 1)}

    for seq in seqs:
        for pos in counts:
            start = (pos - 1) * 3
            counts[pos][seq[start:start + 3]] += 1

    return counts


@pytest.mark.parametrize('aa_seq', (
    'M',
    'MIKEY',
    'SASSAFRAS',
    'MIKEY' * 100,
    'MIKEY' * 500,
))
@pytest.mark.parametrize('codon_weights', (
        CodonWeights.uniform,
        CodonWeights.ecoli,
        CodonWeights.mouse,
))
def test_codon_distributions_roughly_match_weights(aa_seq, codon_weights):
    cw = codon_weights()
    view = CodonGraph(aa_seq, context_l='aaa', context_r='ttt').view(seed=8675309, weights=cw)

    n = 1000

    seqs = [view.sample() for _ in range(n)]
    observed_counts = helper_codon_counts_by_position(seqs)

    pvalues = []

    for i, aa in enumerate(aa_seq):
        pos = i + 1
        expected = {codon: weight * n for codon, weight in view.codon_weights.by_aa(aa).items()}
        observed = observed_counts[pos]

        if len(expected) == 1:
            assert len(observed) == 1
        else:
            result = helper_chi_square_codon_test(observed, expected)
            pvalues.append(result.pvalue)

    if len(pvalues) > 0:

        # Most should pass
        assert sum(p >= 0.001 for p in pvalues) / len(pvalues) >= 0.99

        # No truly terrible results please.
        assert min(pvalues, default=1.0) >= 1e-6


def helper_chi_square_two_sample_test(counts_a, counts_b):
    """
    Check that two sets of samples look like they're drawn from
    the same base distribution.
    """
    codons = sorted({*counts_a.keys(), *counts_b.keys()})

    table = [
        [counts_a.get(codon, 0) for codon in codons],
        [counts_b.get(codon, 0) for codon in codons],
    ]

    return chi2_contingency(table)


@pytest.mark.parametrize('name,aa_seq', NORMAL_PROTEINS.items())
@pytest.mark.parametrize('banned', (
        ('GAATTC', 'GGATCC'),
        ('CTCGAG', 'AAGCTT'),
        ('GAATTC', 'GGATCC', 'CTCGAG', 'AAGCTT')
))
def test_codon_distributions_roughly_match_weights_banned_sequences(name, aa_seq, banned):
    cw = CodonWeights.ecoli()
    view = CodonGraph(aa_seq, context_l='aaa', context_r='ttt').view(seed=8675309, weights=cw)

    n = 10000

    rejection_sampled_seqs = []
    while len(rejection_sampled_seqs) < n:
        seq = view.sample()
        if not any(b in seq for b in banned):
            rejection_sampled_seqs.append(seq)

    bsc = BannedSequenceConstraint(banned)
    view.set_constraints([bsc])

    smart_seqs = [view.sample() for _ in range(n)]

    rejection_counts = helper_codon_counts_by_position(rejection_sampled_seqs)
    smart_counts = helper_codon_counts_by_position(smart_seqs)

    pvalues = []

    for pos in range(1, len(aa_seq) + 1):
        aa = aa_seq[pos - 1]

        if len(view.codon_weights.by_aa(aa)) == 1:
            assert len(rejection_counts[pos]) == len(smart_counts[pos]) == 1
        else:
            result = helper_chi_square_two_sample_test(smart_counts[pos], rejection_counts[pos])
            _, pvalue, _, _ = result
            pvalues.append(pvalue)

    assert sum(p >= 0.001 for p in pvalues) / len(pvalues) >= 0.99
    assert min(pvalues, default=1.0) >= 1e-6


def test_sequence_at_raises_on_unexpected_dead_end():
    view = CodonGraph('M').view()
    view.compile()

    choice_results = list(view.choice_results_by_state_id)
    choice_results[view.initial_state_id] = ()
    view.choice_results_by_state_id = tuple(choice_results)

    with pytest.raises(RuntimeError):
        view.sequence_at(0)


def test_set_weights_marks_view_for_sampler_update():
    view = CodonGraph('MIKEY').view()
    view.compile()

    view.set_weights(CodonWeights.ecoli())

    assert view._compile_status == COMPILE_SHALLOW


def test_set_pinned_codons_marks_view_for_results_update():
    view = CodonGraph('MIKEY').view()
    view.compile()

    view.set_pinned_codons({2: 'ATC'})

    assert view._compile_status == COMPILE_SHALLOW


def test_set_constraints_marks_view_for_topology_update():
    view = CodonGraph('MIKEY').view()
    view.compile()

    view.set_constraints([GCConstraint(min_perc=40)])

    assert view._compile_status == COMPILE_DEEP


def test_compile_requirement_is_not_downgraded():
    view = CodonGraph('MIKEY').view()

    view.set_constraints([GCConstraint(min_perc=40)])
    assert view._compile_status == COMPILE_DEEP

    view.set_weights(CodonWeights.ecoli())
    assert view._compile_status == COMPILE_DEEP


def test_compile_requirement_is_upgraded():
    view = CodonGraph('MIKEY').view()
    view.compile()

    view.set_weights(CodonWeights.ecoli())
    assert view._compile_status == COMPILE_SHALLOW

    view.set_pinned_codons({2: 'ATC'})
    assert view._compile_status == COMPILE_SHALLOW


def test_pin_codons_preserves_deep_topology():
    graph = CodonGraph('MIKEY')
    view = graph.view()

    view.compile()

    states = view._compiled.states
    child_results = view._compiled.child_results_by_state_id

    view.pin_codons({1: 'ATG'})
    view.compile()

    assert view._compiled.states == states
    assert view._compiled.child_results_by_state_id == child_results



def test_unpin_codons_restores_sequence_count():
    graph = CodonGraph('FF')
    view = graph.view()

    original_count = view.n_valid_sequences

    view.pin_codons({1: 'TTT'})
    pinned_count = view.n_valid_sequences

    view.clear_pins()

    assert pinned_count < original_count
    assert view.n_valid_sequences == original_count


def test_set_weights_preserves_deep_topology():
    graph = CodonGraph('MIKEY')
    view = graph.view()

    view.compile()

    states = view._compiled.states
    child_results = view._compiled.child_results_by_state_id

    view.set_weights(CodonWeights.ecoli())
    view.compile()

    assert view._compiled.states == states
    assert view._compiled.child_results_by_state_id == child_results


def test_shallow_compile_clears_cached_samplers():
    graph = CodonGraph('MIKEY')
    view = graph.view()

    view.sample()
    assert any(sampler is not None for sampler in view.samplers_by_state_id)

    view.set_weights(CodonWeights.ecoli())
    view.compile()

    assert all(sampler is None for sampler in view.samplers_by_state_id)


def test_pins_are_applied_on_top_of_constraints():
    graph = CodonGraph('MIKEY')
    view = graph.view(constraints=[GCConstraint(min_perc=40)])

    unconstrained_count = view.n_valid_sequences

    view.pin_codons({2: 'ATT'})

    assert view.n_valid_sequences < unconstrained_count
    assert all(sequence[3:6] == 'ATT' for sequence in view)


class RejectChoicesConstraint(Constraint):
    """
    Test constraint that rejects specified graph choices.
    """

    def __init__(self, rejected_choices):
        self.rejected_choices = frozenset(rejected_choices)

    @property
    def initial_state(self):
        return ()

    @property
    def is_trivial(self):
        return not self.rejected_choices

    def link(self, graph):
        pass

    def advance(self, state, pos, choice):
        if state == SAFE_STATE:
            return SAFE_STATE

        if choice in self.rejected_choices:
            return DEAD_STATE

        return state


def test_add_constraint_before_first_compile_matches_full_compile():
    graph = CodonGraph('MIKEY')
    constraint = RejectChoicesConstraint({'ATA'})

    view = graph.view()
    view.add_constraints(constraint)
    view.compile()

    full_view = graph.view(constraints=[constraint])
    full_view.compile()

    assert view._compiled == full_view._compiled
    assert view._compile_status == COMPILED
    assert view._pending_constraints == ()


def test_add_constraint_extends_existing_compile():
    graph = CodonGraph('MIKEY')
    constraint = RejectChoicesConstraint({'ATA'})

    view = graph.view()
    view.compile()
    base_compiled = view._compiled

    view.add_constraints(constraint)

    assert view.constraints == (constraint,)
    assert view._pending_constraints == (constraint,)
    assert view._compile_status == COMPILE_EXTEND

    view.compile()

    full_view = graph.view(constraints=[constraint])
    full_view.compile()

    assert view._compiled == full_view._compiled
    assert view._compiled is not base_compiled
    assert view._compile_status == COMPILED
    assert view._pending_constraints == ()


def test_add_constraints_extends_with_all_new_constraints():
    graph = CodonGraph('MIKEY')
    constraints = (
        RejectChoicesConstraint({'ATA'}),
        RejectChoicesConstraint({'AAG'}),
    )

    view = graph.view()
    view.compile()
    view.add_constraints(constraints)
    view.compile()

    full_view = graph.view(constraints=constraints)
    full_view.compile()

    assert view._compiled == full_view._compiled
    assert view.constraints == constraints


def test_add_constraints_accumulates_before_compile():
    graph = CodonGraph('MIKEY')
    first = RejectChoicesConstraint({'ATA'})
    second = RejectChoicesConstraint({'AAG'})

    view = graph.view()
    view.compile()
    view.add_constraints(first)
    view.add_constraints(second)

    assert view._pending_constraints == (first, second)
    assert view._compile_status == COMPILE_EXTEND

    view.compile()

    full_view = graph.view(constraints=[first, second])
    full_view.compile()

    assert view._compiled == full_view._compiled


def test_set_constraints_discards_pending_constraints():
    graph = CodonGraph('MIKEY')
    pending = RejectChoicesConstraint({'ATA'})
    replacement = RejectChoicesConstraint({'AAG'})

    view = graph.view()
    view.compile()
    view.add_constraints(pending)
    view.set_constraints([replacement])

    assert view.constraints == (replacement,)
    assert view._pending_constraints == ()
    assert view._compile_status == COMPILE_DEEP

    view.compile()

    full_view = graph.view(constraints=[replacement])
    full_view.compile()

    assert view._compiled == full_view._compiled


def test_add_no_constraints_does_not_require_compile():
    view = CodonGraph('MIKEY').view()
    view.compile()

    compiled = view._compiled
    view.add_constraints([])

    assert view._compiled is compiled
    assert view._compile_status == COMPILED
    assert view._pending_constraints == ()


def test_copy_preserves_pending_constraints():
    constraint = RejectChoicesConstraint({'ATA'})

    view = CodonGraph('MIKEY').view()
    view.compile()
    view.add_constraints(constraint)

    copied = view.copy()

    assert copied.constraints == (constraint,)
    assert copied._pending_constraints == (constraint,)
    assert copied._compile_status == COMPILE_EXTEND

    copied.compile()

    full_view = copied.graph.view(constraints=[constraint])
    full_view.compile()

    assert copied._compiled == full_view._compiled
