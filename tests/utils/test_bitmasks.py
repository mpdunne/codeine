import pytest

from codeine import TranslationTable
from codeine.graph.base import CodonGraph
from codeine.utils.bitmasks import choices_to_nt_bitmasks, sequence_to_nt_bitmasks


@pytest.mark.parametrize(
    ('sequence', 'expected'),
    [
        ('', []),
        ('A', [1]),
        ('C', [2]),
        ('G', [4]),
        ('T', [8]),
        ('U', [8]),
        ('ACGT', [1, 2, 4, 8]),
        ('ACGU', [1, 2, 4, 8]),
        ('AAAA', [1, 1, 1, 1]),
        ('TGCA', [8, 4, 2, 1]),
    ],
)
def test_sequence_to_nt_bitmasks(sequence, expected):
    assert sequence_to_nt_bitmasks(sequence) == expected


@pytest.mark.parametrize(
    ('choices', 'expected'),
    [
        ((('A',),), [1]),
        ((('C',),), [2]),
        ((('G',),), [4]),
        ((('T',),), [8]),
        ((('U',),), [8]),
        ((('',),), []),
        ((('ACGT',),), [1, 2, 4, 8]),
        ((('ACGU',),), [1, 2, 4, 8]),
    ],
)
def test_choices_to_nt_bitmasks_single_choice(choices, expected):
    assert choices_to_nt_bitmasks(choices) == expected


def test_choices_to_nt_bitmasks_alternative_choices():
    assert choices_to_nt_bitmasks((('AAA', 'AAC', 'AAG', 'AAT'),)) == [1, 1, 15]
    assert choices_to_nt_bitmasks((('ATA', 'ATC', 'ATT'),)) == [1, 8, 11]


def test_choices_to_nt_bitmasks_multiple_variable_positions():
    assert choices_to_nt_bitmasks((('AAA', 'CCG'),)) == [3, 3, 5]
    assert choices_to_nt_bitmasks((('TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'),)) == [10, 8, 15]


def test_choices_to_nt_bitmasks_duplicate_choices():
    assert choices_to_nt_bitmasks((('AAA', 'AAA', 'AAG'),)) == [1, 1, 5]


@pytest.mark.parametrize(
    'choices',
    [
        (('A', 'AT'),),
        (('ATG', 'AA'),),
        (('AAA',), ('CCC', 'CC')),
    ],
)
def test_unequal_choice_lengths(choices):
    with pytest.raises(ValueError, match='Choices within each choice set must have equal length'):
        choices_to_nt_bitmasks(choices)


def test_choices_to_nt_bitmasks_all_bases():
    assert choices_to_nt_bitmasks((('A', 'C', 'G', 'T'),)) == [15]


def test_choices_to_nt_bitmasks_dna_and_rna():
    assert choices_to_nt_bitmasks((('A', 'C', 'G', 'T', 'U'),)) == [15]


def test_choices_to_nt_bitmaskss_multiple_choice_sets():
    assert choices_to_nt_bitmasks((('ATG',), ('AAA', 'AAG'), ('TTC', 'TTT'))) == [
        1, 8, 4,
        1, 1, 5,
        8, 8, 10,
    ]


def test_choices_to_nt_bitmasks_nonempty_contexts():
    choices = (
        ('ATCGACCA',),
        ('ATA', 'ATC', 'ATT'),
        ('AAA', 'AAG'),
        ('CCGTTAG',),
    )

    assert choices_to_nt_bitmasks(choices) == [
        1, 8, 2, 4, 1, 2, 2, 1,
        1, 8, 11,
        1, 1, 5,
        2, 2, 4, 8, 8, 1, 4,
    ]


def test_choices_to_nt_bitmasks_empty_contexts():
    assert choices_to_nt_bitmasks((('',), ('ATG',), ('',))) == [1, 8, 4]


def test_choices_to_nt_bitmasks_long_contexts():
    sequence = 'ATCGACCAATTGTTCCATCG'

    assert choices_to_nt_bitmasks(((sequence,),)) == sequence_to_nt_bitmasks(sequence)


def test_fixed_choice_matches_sequence_masks():
    sequence = 'ATCGACCAATTG'

    assert choices_to_nt_bitmasks(((sequence,),)) == sequence_to_nt_bitmasks(sequence)


def test_each_choice_is_subset_of_choice_masks():
    choices = (('TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'),)
    nt_masks = choices_to_nt_bitmasks(choices)

    for sequence in choices[0]:
        sequence_masks = sequence_to_nt_bitmasks(sequence)

        assert all(sequence_mask & choice_mask == sequence_mask
                   for sequence_mask, choice_mask in zip(sequence_masks, nt_masks))


def test_generated_sequences_are_subsets_of_choice_masks():
    graph = CodonGraph('MIKEY')
    view = graph.view()

    choices = tuple(tuple(node.transitions) for node in graph.codon_nodes)
    nt_masks = choices_to_nt_bitmasks(choices)

    for sequence in view.enumerate():
        sequence_masks = sequence_to_nt_bitmasks(sequence)

        assert len(sequence_masks) == len(nt_masks)
        assert all(sequence_mask & choice_mask == sequence_mask
                   for sequence_mask, choice_mask in zip(sequence_masks, nt_masks))
