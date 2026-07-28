import pytest

from codeine.constraints.tandem import TandemRepeatConstraint
from codeine.graph.base import CodonGraph
from codeine.translation.tables import TranslationTable

from tests.data import NORMAL_PROTEINS, DIFFICULT_PROTEINS, ANTIBODIES


###############################
# Helpers
###############################


def helper_has_tandem_repeat(
    sequence,
    repeat_length,
    min_copies,
):
    """
    Return whether a sequence contains the specified exact tandem repeat.
    """
    repeat_span = repeat_length * min_copies

    for start in range(len(sequence) - repeat_span + 1):
        repeat = sequence[start:start + repeat_length]

        if sequence[start:start + repeat_span] == repeat * min_copies:
            return True

    return False


@pytest.mark.parametrize(
    'sequence,repeat_length,min_copies,expected',
    [
        ('AAAA', 1, 2, True),
        ('AAAA', 1, 4, True),
        ('AAAT', 1, 4, False),
        ('ATAT', 2, 2, True),
        ('ATATAT', 2, 2, True),
        ('ATATAT', 2, 3, True),
        ('ATATAC', 2, 3, False),
        ('TATAT', 2, 2, True),
        ('ACTACT', 3, 2, True),
        ('ACTACTACT', 3, 2, True),
        ('ACTACTACT', 3, 3, True),
        ('ACTACTACC', 3, 3, False),
        ('AATAATAAT', 3, 3, True),
        ('CACTACTA', 3, 2, True),
        ('ACGTACGT', 4, 2, True),
        ('ACGTACGTACGT', 4, 3, True),
        ('ACGTACGA', 4, 2, False),
        ('AT', 2, 2, False),
        ('', 1, 2, False),
    ],
)
def test_has_tandem_repeat(sequence, repeat_length, min_copies, expected):
    assert helper_has_tandem_repeat(sequence, repeat_length, min_copies) is expected


def helper_enumerate_expected(
    aa_seq,
    repeat_length,
    min_copies,
    context_l='',
    context_r='',
):
    """
    Brute-force the sequences that should remain after applying the constraint.
    """
    view = CodonGraph(aa_seq, context_l=context_l, context_r=context_r).view()

    return {
        sequence for sequence in view.enumerate()
        if not helper_has_tandem_repeat(context_l + sequence + context_r, repeat_length, min_copies)
    }


###############################
# Construction and validation
###############################


def test_tandem_repeat_constraint_stores_parameters():
    constraint = TandemRepeatConstraint(repeat_length=4, min_copies=3)

    assert constraint.repeat_length == 4
    assert constraint.min_copies == 3


def test_tandem_repeat_constraint_defaults_to_two_copies():
    constraint = TandemRepeatConstraint(repeat_length=4)

    assert constraint.min_copies == 2


@pytest.mark.parametrize('repeat_length', [0, -1])
def test_repeat_length_must_be_positive(repeat_length):
    with pytest.raises(ValueError, match='repeat_length must be at least 1'):
        TandemRepeatConstraint(repeat_length=repeat_length)


@pytest.mark.parametrize('repeat_length', [None, 1.5, '4', True])
def test_repeat_length_must_be_an_integer(repeat_length):
    with pytest.raises(TypeError, match='repeat_length must be an integer'):
        TandemRepeatConstraint(repeat_length=repeat_length)


@pytest.mark.parametrize('min_copies', [0, 1, -1])
def test_min_copies_must_be_at_least_two(min_copies):
    with pytest.raises(ValueError, match='min_copies must be at least 2'):
        TandemRepeatConstraint(repeat_length=4, min_copies=min_copies)


@pytest.mark.parametrize('min_copies', [None, 2.5, '2', True])
def test_min_copies_must_be_an_integer(min_copies):
    with pytest.raises(TypeError, match='min_copies must be an integer'):
        TandemRepeatConstraint(repeat_length=4, min_copies=min_copies)


###############################
# Functionality
###############################

@pytest.mark.parametrize(
    'aa_seq,repeat_length,min_copies',
    [
        # Single-nucleotide repeat units.
        ('F', 1, 2),
        ('FF', 1, 3),

        # Codon-aligned repeats.
        ('FF', 3, 2),
        ('FFF', 3, 3),
        ('FFFF', 3, 4),

        # Repeat units shorter than a codon.
        ('FF', 2, 2),
        ('FFF', 2, 3),
        ('FFFF', 2, 4),

        # Repeat units crossing codon boundaries.
        ('FIY', 4, 2),
        ('FYFI', 5, 2),
        ('FIFI', 7, 2),

        # Repeat units longer than one codon.
        ('FFFF', 4, 2),
        ('FFFF', 6, 2),
        ('FFFFF', 8, 2),

        # Different codon degeneracies.
        ('II', 2, 2),
        ('LL', 3, 2),
        ('LIY', 5, 2),

        # No possible repeats.
        ('ME', 2, 2),
        ('MKW', 3, 2),

        # Repeat span longer than the coding sequence.
        ('F', 2, 2),
        ('FI', 4, 2),
        ('FI', 3, 3),
    ],
)
def test_enumerated_sequences_do_not_contain_tandem_repeats(aa_seq, repeat_length, min_copies):
    expected = helper_enumerate_expected(aa_seq, repeat_length, min_copies)

    view = CodonGraph(aa_seq).view()

    constraint = TandemRepeatConstraint(repeat_length=repeat_length, min_copies=min_copies)
    view.set_constraints([constraint])

    observed = set(view.enumerate())

    assert observed == expected
    assert view.n_valid_sequences == len(expected)


###############################
# Large spaces
###############################

@pytest.mark.parametrize(
    'aa_seq,repeat_length,min_copies',
    [
        ('F' * 20, 1, 6),
        ('F' * 20, 2, 4),
        ('F' * 20, 3, 3),
        ('F' * 20, 4, 3),
        ('F' * 20, 5, 2),
        ('F' * 20, 6, 2),
        ('IY' * 20, 1, 6),
        ('IY' * 20, 2, 4),
        ('IY' * 20, 3, 3),
        ('IY' * 20, 4, 3),
        ('IY' * 20, 5, 2),
        ('IY' * 20, 6, 2),
        ('FIL' * 15, 7, 2),
        ('FIL' * 15, 8, 2),
    ],
)
def test_sampled_sequences_do_not_contain_tandem_repeats(aa_seq, repeat_length, min_copies):
    view = CodonGraph(aa_seq).view(seed=8675309)

    constraint = TandemRepeatConstraint(repeat_length=repeat_length, min_copies=min_copies)
    view.set_constraints([constraint])

    for sequence in view.sample(n=1000):
        assert not helper_has_tandem_repeat(sequence, repeat_length, min_copies)


###############################
# Contexts
###############################


@pytest.mark.parametrize(
    'aa_seq,context_l,context_r,repeat_length,min_copies',
    [
        # Repeat begins in the left context.
        ('F', 'T', '', 2, 2),       # T + TTT -> TT x 2
        ('F', 'TTC', '', 3, 2),     # TTC + TTC -> TTC x 2
        ('FI', 'AT', '', 4, 2),
        ('FFFF', 'AT', '', 7, 2),

        # Repeat ends in the right context.
        ('F', '', 'T', 2, 2),       # TTT + T -> TT x 2
        ('F', '', 'TTC', 3, 2),     # TTC + TTC -> TTC x 2
        ('FI', '', 'AT', 4, 2),
        ('FFFF', '', 'AT', 7, 2),

        # Different synonymous encodings violate opposite boundaries.
        ('F', 'T', 'ACA', 2, 2),    # T + TTT -> TT x 2 /  TTC + ACA -> CA x 2

        # Repeat units cross codon boundaries.
        ('FF', 'CT', '', 4, 2),     # CT + TTCTTT -> CTTT x 2
        ('FF', '', 'CT', 4, 2),     # TTCTTT + CT -> TTCT x 2
        ('FIY', 'A', 'T', 5, 2),
        ('FIY', 'AT', 'GC', 7, 2),

        # Single-nucleotide repeats cross context boundaries.
        ('F', 'T', '', 1, 4),
        ('F', '', 'T', 1, 4),

        # Three-copy repeats cross context boundaries.
        ('F', 'TT', '', 2, 3),
        ('F', '', 'TT', 2, 3),

        # Contexts shorter and longer than codons.
        ('FI', 'A', 'C', 2, 2),
        ('FI', 'AT', 'CG', 3, 2),
        ('FI', 'ATGC', 'CGTA', 4, 2),
        ('FIYF', 'AT', 'GC', 5, 2),
    ],
)
def test_tandem_repeats_in_contexts(
    aa_seq,
    context_l,
    context_r,
    repeat_length,
    min_copies,
):
    expected = helper_enumerate_expected(
        aa_seq,
        repeat_length,
        min_copies,
        context_l=context_l,
        context_r=context_r,
    )

    view = CodonGraph(aa_seq, context_l=context_l, context_r=context_r).view()

    constraint = TandemRepeatConstraint(repeat_length=repeat_length, min_copies=min_copies)
    view.set_constraints([constraint])

    observed = set(view.enumerate())

    assert observed == expected
    assert view.n_valid_sequences == len(expected)


@pytest.mark.parametrize(
    'context_l,context_r',
    [
        ('ATAT', ''),
        ('', 'ATAT'),
        ('TATA', ''),
        ('', 'TATA'),
        ('ACAC', ''),
        ('', 'ACAC'),
        ('CGCG', ''),
        ('', 'CGCG'),
        ('ATATAT', ''),
        ('', 'ATATAT'),
        ('GCGCGC', ''),
        ('', 'GCGCGC'),
        ('GCTCATATATGCTC', ''),
        ('', 'GCTCATATATGCTC'),
        ('ATTAGCGCGCATTA', ''),
        ('', 'ATTAGCGCGCATTA'),
    ],
)
def test_repeat_entirely_within_context_makes_space_empty(context_l, context_r):
    view = CodonGraph('MIKEY', context_l=context_l, context_r=context_r).view()

    constraint = TandemRepeatConstraint(repeat_length=2, min_copies=2)
    view.set_constraints([constraint])

    assert view.n_valid_sequences == 0
    assert list(view.enumerate()) == []


##################
# Real life.
##################

@pytest.mark.parametrize(
    'aa_seq,repeat_length,min_copies',
    [
        (NORMAL_PROTEINS['ubiquitin'], 3, 3),
        (NORMAL_PROTEINS['gfp'], 6, 2),
        (NORMAL_PROTEINS['mcherry'], 4, 3),
        (DIFFICULT_PROTEINS['collagen'], 3, 3),
        (DIFFICULT_PROTEINS['also_hard'], 6, 2),
        (ANTIBODIES['caplacizumab'], 5, 2),
    ],
)
def test_real_protein_samples_do_not_contain_tandem_repeats(
    aa_seq,
    repeat_length,
    min_copies,
):
    view = CodonGraph(aa_seq).view(seed=8675309)

    constraint = TandemRepeatConstraint(repeat_length=repeat_length, min_copies=min_copies)
    view.set_constraints([constraint])

    for sequence in view.sample(n=100):
        assert not helper_has_tandem_repeat(sequence, repeat_length, min_copies)


@pytest.mark.parametrize(
    'aa_seq,context_l,context_r,repeat_length,min_copies',
    [
        (NORMAL_PROTEINS['ubiquitin'], 'ATAT', '', 4, 2),
        (ANTIBODIES['caplacizumab'], '', 'GCGCGC', 3, 2),
    ],
)
def test_real_protein_samples_with_contexts(
    aa_seq,
    context_l,
    context_r,
    repeat_length,
    min_copies,
):
    view = CodonGraph(aa_seq, context_l=context_l, context_r=context_r).view(seed=8675309)

    constraint = TandemRepeatConstraint(repeat_length=repeat_length, min_copies=min_copies)
    view.set_constraints([constraint])

    for sequence in view.sample(n=100):
        full_sequence = context_l + sequence + context_r

        assert not helper_has_tandem_repeat(full_sequence, repeat_length, min_copies)


########################
# Multiple constraints
########################


def test_multiple_tandem_repeat_constraints():
    aa_seq = 'FIYFI'

    constraints = [
        TandemRepeatConstraint(2, 3),
        TandemRepeatConstraint(3, 2),
        TandemRepeatConstraint(5, 2),
    ]

    expected = {
        sequence
        for sequence in CodonGraph(aa_seq).view().enumerate()
        if all(
            not helper_has_tandem_repeat(sequence, constraint.repeat_length, constraint.min_copies)
            for constraint in constraints
        )
    }

    view = CodonGraph(aa_seq).view()
    view.set_constraints(constraints)

    observed = set(view.enumerate())

    assert observed == expected
    assert view.n_valid_sequences == len(expected)


@pytest.mark.parametrize(
    'aa_seq',
    [
        NORMAL_PROTEINS['ubiquitin'],
        NORMAL_PROTEINS['gfp'],
        NORMAL_PROTEINS['mcherry'],
        DIFFICULT_PROTEINS['collagen'],
        ANTIBODIES['caplacizumab'],
    ],
)
def test_sampled_sequences_satisfy_multiple_tandem_repeat_constraints(aa_seq):
    constraints = [
        TandemRepeatConstraint(1, 6),
        TandemRepeatConstraint(2, 4),
        TandemRepeatConstraint(3, 3),
        TandemRepeatConstraint(4, 3),
        TandemRepeatConstraint(5, 2),
        TandemRepeatConstraint(6, 2),
    ]

    view = CodonGraph(aa_seq).view(seed=8675309)
    view.set_constraints(constraints)

    for sequence in view.sample(n=100):
        for constraint in constraints:
            assert not helper_has_tandem_repeat(sequence, constraint.repeat_length, constraint.min_copies)


#################
# Longer repeats
#################

@pytest.mark.parametrize(
    'seq,repeat_length,min_copies',
    [
        ('ATATATATATATATAT' + 'TATATATATATATATA' + 'T', 15, 2),
        ('C' + 'ATATATATATATATAT' + 'TATATATATATATATA', 15, 2),
        ('ATC' * 6 + 'ATC' * 6 + 'ATC' * 6 + 'ATC', 18, 3),
        ('T' + 'ATC' * 6 + 'ATC' * 6 + 'ATC' * 6 + 'ATCTT', 18, 3),
    ],
)
def test_long_and_possible_repeat_unit(seq, repeat_length, min_copies):
    aa_seq = TranslationTable().translate(seq)

    view = CodonGraph(aa_seq).view()
    unconstrained_n = view.n_valid_sequences

    view.set_constraints([TandemRepeatConstraint(repeat_length=repeat_length, min_copies=min_copies)])
    view.compile()
    constrained_n = view.n_valid_sequences

    assert constrained_n < unconstrained_n


@pytest.mark.parametrize(
    'seq,repeat_length,min_copies',
    [
        ('AT' * 99 + 'CGA', 100, 2),
        ('ATC' * 100 + 'ATC' * 100 + 'ATC' * 99 + 'CCC', 100, 3),
        ('ATGC' * 1000 + 'ATGC' * 999 + 'AAAAA', 4000, 3),
        ('ATCGACGTAGGCATGCCGTAC' * 9 + 'A' * 21, 21, 10),
        ('ATCGACGTAGGCATGCCGTAC' * 99 + 'A' * 21, 21, 100),
    ],
)
def test_long_but_impossible_repeat_unit(seq, repeat_length, min_copies):
    aa_seq = TranslationTable().translate(seq)

    view = CodonGraph(aa_seq).view()
    unconstrained_n = view.n_valid_sequences

    view.set_constraints([TandemRepeatConstraint(repeat_length=repeat_length, min_copies=min_copies)])
    view.compile()
    constrained_n = view.n_valid_sequences

    assert constrained_n == unconstrained_n
