import random
import pytest

from codeine.constraints.repeats import DirectRepeatConstraint, InvertedRepeatConstraint, RepeatConstraint
from codeine.graph.base import CodonGraph
from codeine.translation.tables import TranslationTable
from codeine.tools.repeats import contains_direct_repeat, contains_inverted_repeat


###############################
# Validation & construction
###############################


class ExampleRepeatConstraint(RepeatConstraint):
    pass


def test_repeat_constraint_stores_args_correctly():
    constraint = ExampleRepeatConstraint(repeat_length=10, min_distance=3, max_distance=20, inverted=True)

    assert constraint.repeat_length == 10
    assert constraint.min_distance == 3
    assert constraint.max_distance == 20
    assert constraint.inverted is True


def test_repeat_constraint_uses_defaults():
    constraint = ExampleRepeatConstraint(repeat_length=10)

    assert constraint.repeat_length == 10
    assert constraint.min_distance == 0
    assert constraint.max_distance is None
    assert constraint.inverted is False


@pytest.mark.parametrize('repeat_length', [None, 1.0, '1', [], True, False])
def test_repeat_length_must_be_integer(repeat_length):
    with pytest.raises(TypeError, match='repeat_length must be an integer'):
        ExampleRepeatConstraint(repeat_length)


@pytest.mark.parametrize('repeat_length', [0, -1, -100])
def test_repeat_length_must_be_positive(repeat_length):
    with pytest.raises(ValueError, match='repeat_length must be at least 1'):
        ExampleRepeatConstraint(repeat_length)


@pytest.mark.parametrize('min_distance', [None, 1.0, '1', [], True, False])
def test_min_distance_must_be_integer(min_distance):
    with pytest.raises(TypeError, match='min_distance must be an integer'):
        ExampleRepeatConstraint(repeat_length=10, min_distance=min_distance)


@pytest.mark.parametrize('min_distance', [-1, -100])
def test_min_distance_must_be_non_negative(min_distance):
    with pytest.raises(ValueError, match='min_distance must be at least 0'):
        ExampleRepeatConstraint(repeat_length=10, min_distance=min_distance)


@pytest.mark.parametrize('max_distance', [1.0, '10', [], True, False])
def test_max_distance_must_be_integer_or_none(max_distance):
    with pytest.raises(TypeError, match='max_distance must be an integer or None'):
        ExampleRepeatConstraint(repeat_length=10, max_distance=max_distance)


def test_max_distance_must_not_be_less_than_minimum():
    with pytest.raises(ValueError, match='max_distance must be at least min_distance'):
        ExampleRepeatConstraint(repeat_length=10, min_distance=5, max_distance=4)


@pytest.mark.parametrize('inverted', [None, 0, 1, 'yes'])
def test_inverted_must_be_boolean(inverted):
    with pytest.raises(TypeError, match='inverted must be a boolean'):
        ExampleRepeatConstraint(repeat_length=10, inverted=inverted)


def test_direct_repeat_stores_values_correctly():
    constraint = DirectRepeatConstraint(repeat_length=10, min_distance=3, max_distance=20)

    assert constraint.repeat_length == 10
    assert constraint.min_distance == 3
    assert constraint.max_distance == 20
    assert constraint.inverted is False


def test_inverted_repeat_stores_values_correctly():
    constraint = InvertedRepeatConstraint(repeat_length=10, min_distance=3, max_distance=20)

    assert constraint.repeat_length == 10
    assert constraint.min_distance == 3
    assert constraint.max_distance == 20
    assert constraint.inverted is True


###################################################
# Trivial constraints and repeat-free sequences
###################################################


def make_repeat_safe_sequence(length, repeat_length, alphabet='AC'):
    """
    Construct a sequence containing no repeated repeat_length-mers.
    """
    if not alphabet:
        raise ValueError('alphabet must not be empty')

    if length > len(alphabet) ** repeat_length + repeat_length - 1:
        raise ValueError('requested sequence is too long for this repeat length')

    sequence = [alphabet[0]] * (repeat_length - 1)
    seen = set()

    while len(sequence) < length:
        for nt in alphabet[::-1]:
            kmer = ''.join(sequence[-repeat_length + 1:] + [nt])

            if kmer not in seen:
                seen.add(kmer)
                sequence.append(nt)
                break
        else:
            raise RuntimeError('Could not construct repeat-safe sequence')

    return ''.join(sequence)


def reverse_complement(sequence):
    return sequence.translate(str.maketrans('ACGT', 'TGCA'))[::-1]


@pytest.mark.parametrize(
    ('repeat_length', 'context_l_size', 'context_r_size'),
    [
        (12, 0, 100),
        (24, 50, 50),
        (120, 100, 200),
    ],
)
def test_repeat_safe_contexts_have_no_direct_repeats(
    repeat_length,
    context_l_size,
    context_r_size,
):
    sequence = make_repeat_safe_sequence(context_l_size + context_r_size, repeat_length)

    context_l = sequence[:context_l_size]
    context_r = sequence[context_l_size:]

    graph = CodonGraph('M', context_l=context_l, context_r=context_r)

    constraint = RepeatConstraint(repeat_length=repeat_length)
    constraint.link(graph)

    assert constraint.is_trivial


def test_repeated_coding_region_can_form_direct_repeat():
    constraint = RepeatConstraint(repeat_length=15)
    constraint.link(CodonGraph('MIKEYAAAAAMIKEY'))

    assert not constraint.is_trivial


def test_fixed_codons_kill_direct_repeat():
    constraint = RepeatConstraint(repeat_length=15)
    graph = CodonGraph('MIKEYAAAAAMIKEY')
    constraint.link(graph)

    assert not constraint.is_trivial

    constraint = RepeatConstraint(repeat_length=15)
    graph = CodonGraph('MIKEYAAAAAMIKEY', codon_restrictions={2: 'ATA', 12: 'ATT'})
    constraint.link(graph)

    assert constraint.is_trivial


def test_exact_direct_repeat_in_context():
    repeat = 'ATGATCAAAGAATAC'

    graph = CodonGraph('M', context_l=repeat + 'ACACACAC' + repeat)

    constraint = RepeatConstraint(repeat_length=len(repeat))
    constraint.link(graph)

    assert not constraint.is_trivial


def test_single_change_kills_direct_repeat_in_context():
    repeat = 'ATGATCAAAGAATAC'
    modified_repeat = 'ATGATCAAGGAATAC'

    graph = CodonGraph('M', context_l=repeat + 'ACACACAC' + modified_repeat)

    constraint = RepeatConstraint(repeat_length=len(repeat))
    constraint.link(graph)

    assert constraint.is_trivial


@pytest.mark.parametrize(
    'min_distance,max_distance',
    [
        (6, 12),
        (12, 24),
        (18, 36),
    ],
)
def test_direct_repeat_distance(min_distance, max_distance):
    constraint = RepeatConstraint(repeat_length=15, min_distance=min_distance, max_distance=max_distance)

    repeat_free_cds = 'AAAAAAAAAAAAAACCCCCCCCCCCCCCCACCCCCCCCCCCCCAACCC'
    repeat_free_aa = 'KKKKNPPPPPPPPPQP'

    def make_graph(linker_aa_len):
        linker_aa = repeat_free_aa[:linker_aa_len]
        linker_cds = repeat_free_cds[:linker_aa_len * 3]
        aa_seq = 'CCCCC' + linker_aa + 'CCCCC'
        codon_restrictions = {6 + ix: linker_cds[ix * 3:(ix + 1) * 3] for ix in range(linker_aa_len)}
        return CodonGraph(aa_seq, codon_restrictions=codon_restrictions)

    # Too close :(
    graph = make_graph(min_distance // 3 - 1)
    constraint.link(graph)
    assert constraint.is_trivial

    # Just right! :D
    graph = make_graph(min_distance // 3)
    constraint.link(graph)
    assert not constraint.is_trivial

    # Too far :O
    graph = make_graph(max_distance // 3 + 1)
    constraint.link(graph)
    assert constraint.is_trivial


@pytest.mark.parametrize(
    ('aa_seq', 'repeat_length'),
    [
        ('MIKEYAAAAAMIKEY', 15),
        ('MKTLEFQAAAAAMKTLEFQ', 21),
        ('GCPRYKKLLLLLGCPRYKK', 21),
    ],
)
def test_repeated_coding_regions_are_not_trivial(aa_seq, repeat_length):
    graph = CodonGraph(aa_seq)

    constraint = RepeatConstraint(repeat_length=repeat_length)
    constraint.link(graph)

    assert not constraint.is_trivial


@pytest.mark.parametrize(
    ('repeat_length', 'aa_length', 'context_l_size', 'context_r_size'),
    [
        (6, 10, 0, 0),
        (24, 50, 20, 30),
        (120, 100, 100, 50),
    ],
)
def test_ag_only_coding_spaces_have_no_inverted_repeats(
    repeat_length,
    aa_length,
    context_l_size,
    context_r_size,
):
    # K and E are encoded only by A/G codons:
    # K = AAA/AAG and E = GAA/GAG.
    # Their reverse complements contain only T/C, so an A/G-only coding
    # space cannot contain an inverted repeat.
    for seed in range(10):
        rng = random.Random(seed)

        aa_seq = ''.join(rng.choice('KE') for _ in range(aa_length))
        context_l = ''.join(rng.choice('AG') for _ in range(context_l_size))
        context_r = ''.join(rng.choice('AG') for _ in range(context_r_size))

        graph = CodonGraph(aa_seq, context_l=context_l, context_r=context_r)

        constraint = RepeatConstraint(repeat_length=repeat_length, inverted=True,)
        constraint.link(graph)

        assert constraint.is_trivial


@pytest.mark.parametrize(
    'repeat',
    [
        'ATGATCAAAGAATAC',
        'ATCGACCAATTG',
    ],
)
def test_exact_inverted_repeat(repeat):
    graph = CodonGraph('M', context_l=repeat + 'AGAGAGAG' + reverse_complement(repeat))

    constraint = RepeatConstraint(repeat_length=len(repeat), inverted=True)
    constraint.link(graph)

    assert not constraint.is_trivial


def test_fixed_codons_kill_inverted_repeat():
    aa_seq = 'KEEKEEFLFLFLKEEEKEEKEEKEEKEE'  # <- FL's encodings can form inverted repeats with the KE encodings.

    constraint = RepeatConstraint(repeat_length=15, inverted=True)
    graph = CodonGraph(aa_seq)
    constraint.link(graph)

    assert not constraint.is_trivial

    constraint = RepeatConstraint(repeat_length=15, inverted=True)
    graph = CodonGraph(aa_seq, codon_restrictions={8: 'TTG'})  # This breaks the inverted repeats.
    constraint.link(graph)

    assert constraint.is_trivial


def test_single_change_kills_inverted_repeat_in_context():
    repeat = 'ATGATCAAAGAATAC'
    modified_repeat = 'ATGATCAAGGAATAC'

    graph = CodonGraph('M', context_l=repeat + 'AGAGAGAG' + reverse_complement(modified_repeat))

    constraint = RepeatConstraint(repeat_length=len(repeat), inverted=True)
    constraint.link(graph)

    assert constraint.is_trivial


@pytest.mark.parametrize(
    'min_distance,max_distance',
    [
        (6, 12),
        (12, 24),
        (18, 36),
    ],
)
def test_inverted_repeat_distance(min_distance, max_distance):
    constraint = RepeatConstraint(
        repeat_length=15,
        min_distance=min_distance,
        max_distance=max_distance,
        inverted=True,
    )

    repeat_free_cds = 'AAAAAAAAAAAAAACCCCCCCCCCCCCCCACCCCCCCCCCCCCAACCC'
    repeat_free_aa = 'KKKKNPPPPPPPPPQP'

    def make_graph(linker_aa_len):
        linker_aa = repeat_free_aa[:linker_aa_len]
        linker_cds = repeat_free_cds[:linker_aa_len * 3]
        aa_seq = 'KEEKE' + linker_aa + 'FFFFF'
        codon_restrictions = {6 + ix: linker_cds[ix * 3:(ix + 1) * 3] for ix in range(linker_aa_len)}
        return CodonGraph(aa_seq, codon_restrictions=codon_restrictions)

    # Too close :(
    linker_aa_len = min_distance // 3 - 1
    constraint.link(make_graph(linker_aa_len))
    repeat_ixs = (0, 15 + linker_aa_len * 3)
    found_ixs = [(start_l, start_r) for start_l, start_r, _requirements in constraint.repeats]
    assert repeat_ixs not in found_ixs

    # Just right! :D
    linker_aa_len = min_distance // 3
    constraint.link(make_graph(linker_aa_len))
    repeat_ixs = (0, 15 + linker_aa_len * 3)
    found_ixs = [(start_l, start_r) for start_l, start_r, _requirements in constraint.repeats]
    assert repeat_ixs in found_ixs

    # Too far :(
    linker_aa_len = max_distance // 3 + 1
    constraint.link(make_graph(linker_aa_len))
    repeat_ixs = (0, 15 + linker_aa_len * 3)
    found_ixs = [(start_l, start_r) for start_l, start_r, _requirements in constraint.repeats]
    assert repeat_ixs not in found_ixs


@pytest.mark.parametrize(
    'repeat',
    [
        'ATGATCAAAGAATAC',
        'ATCGACCAATTG',
        'GTCAGGATCCGATGCAAT',
    ],
)
def test_inverted_repeats_are_not_trivial(repeat):
    graph = CodonGraph('M', context_l=repeat + 'AGAGAGAG' + reverse_complement(repeat))

    constraint = RepeatConstraint(repeat_length=len(repeat), inverted=True)
    constraint.link(graph)

    assert not constraint.is_trivial


##################
# Repeats exist!
##################


@pytest.mark.parametrize(
    'aa_seq,repeat_length',
    [
        ('II', 3),
        ('IYIY', 6),
        ('MIKEYMIKEY', 15),
        ('MKTLEFQMKTLEFQ', 21),
    ],
)
def test_finds_direct_repeats(aa_seq, repeat_length):
    constraint = RepeatConstraint(repeat_length)
    constraint.link(CodonGraph(aa_seq))

    assert constraint.repeats


@pytest.mark.parametrize(
    'repeat',
    [
        'ATGATC',
        'ATCGACCAATTG',
        'ATGATCAAAGAATAC',
    ],
)
def test_finds_inverted_repeats(repeat):
    constraint = RepeatConstraint(len(repeat), inverted=True)
    constraint.link(CodonGraph('M', context_l=repeat + 'AGAGAGAG' + reverse_complement(repeat)))

    assert constraint.repeats


def test_finds_long_direct_repeat():
    constraint = RepeatConstraint(1500)
    constraint.link(CodonGraph('MIKEY' * 200))

    assert len(constraint.repeats) == 1


def test_finds_long_inverted_repeat():
    constraint = RepeatConstraint(1500, inverted=True)
    constraint.link(CodonGraph('KE' * 250 + 'FL' * 250))

    assert len(constraint.repeats) == 1


def test_finds_distant_direct_repeat():
    repeat_length = 15
    spacer = make_repeat_safe_sequence(3000, repeat_length)
    tt = TranslationTable()

    aa_seq = 'MIKEY' + tt.translate(spacer) + 'MIKEY'
    codon_restrictions = {6 + ix: spacer[ix * 3:(ix + 1) * 3] for ix in range(len(spacer) // 3)}
    graph = CodonGraph(aa_seq, codon_restrictions=codon_restrictions)

    constraint = RepeatConstraint(repeat_length)
    constraint.link(graph)

    assert len(constraint.repeats) == 1


def test_finds_distant_inverted_repeat():
    repeat_length = 15
    spacer = 'AC' * 1500
    tt = TranslationTable()

    aa_seq = 'KEEKE' + tt.translate(spacer) + 'FFFFF'
    codon_restrictions = {6 + ix: spacer[ix * 3:(ix + 1) * 3] for ix in range(len(spacer) // 3)}
    graph = CodonGraph(aa_seq, codon_restrictions=codon_restrictions)

    constraint = RepeatConstraint(repeat_length, inverted=True)
    constraint.link(graph)

    assert len(constraint.repeats) == 1


@pytest.mark.parametrize(
    'context_l,cds,context_r,repeat_length,expected',
    [
        ('ATG', 'AACGCA', 'T', 3, [(0, 7)]),
        ('AAAAATG', 'AACGCA', 'T', 3, [(4, 11)]),
        ('', 'ATGAAC', 'GCAT', 3, [(0, 7)]),
        ('AAAAA', 'ATGAAC', 'GCAT', 3, [(5, 12)]),
        ('ATG', 'AACGCA', 'TAAAAA', 3, [(0, 7)]),
        ('AAAAATG', 'AACGCA', 'TAAAAA', 3, [(4, 11)]),
    ],
)
def test_finds_inverted_repeats_across_full_sequence(context_l, cds, context_r, repeat_length, expected):

    # This test exists to catch a previous bug in which we were missing matches whose rev. comp.
    # appeared before the first part of the match, in the rev. comped sequence.
    #
    # reference:       A A A A A T G A A C G C A T
    #                          A T G
    #
    # rev comp:        A T G C G T T C A T T T T T  ---> caught by a right-moving scan
    #                  A T G
    #
    # We need to scan both ways, unlike direct repeats :)
    #
    # Compare with this example.
    #
    # reference:       A T G A A C G C A T A A A A
    #                  A T G
    #
    # rev comp:        T T T T A T G C G T T C A T  <--- caught by a left-moving scan
    #                          A T G

    tt = TranslationTable()
    aa_seq = tt.translate(cds)
    codon_restrictions = {ix + 1: cds[ix * 3:(ix + 1) * 3] for ix in range(len(cds) // 3)}

    graph = CodonGraph(aa_seq, context_l=context_l, context_r=context_r, codon_restrictions=codon_restrictions)

    constraint = RepeatConstraint(repeat_length, inverted=True)
    constraint.link(graph)

    observed = [(start_l, start_r) for start_l, start_r, _requirements in constraint.repeats]
    assert observed == expected


##########################################
# Enumeration tests & empty spaces
##########################################

@pytest.mark.parametrize(
    'context_l,aa_seq,context_r,repeat_length',
    [
        # In-frame repeats in CDS
        ('GCTAACGTAC', 'MIKEYMIKEY', 'TCGATGCAGT', 15),
        ('ATCGACCAATTG', 'IYIY', 'GCTAACGTAC', 6),
        ('GCTAACGTAC', 'IYAIYAIY', 'TCGATGCAGT', 6),

        # Out-of-frame repeats in CDS
        ('', 'RNYKQT', '', 4),
        ('', 'IHERQW', '', 4),
        ('', 'MTKPHY', '', 4),

        # Repeats overlapping contexts
        ('GGACGTT', 'KRALE', 'TCTCTTCC', 4),  # overlaps context_l
        ('TATCC', 'DNRVG', 'CAAG', 4),  # overlaps context_r
        ('CACTGCG', 'RIWGS', 'AGGGGA', 5),  # overlaps context_r
    ],
)
def test_direct_repeat_constraint_enumeration(context_l, aa_seq, context_r, repeat_length):
    view = CodonGraph(aa_seq, context_l=context_l, context_r=context_r).view()
    unconstrained = [*view.enumerate()]

    view.set_constraints([DirectRepeatConstraint(repeat_length)])
    constrained = [*view.enumerate()]

    filtered = [cds for cds in unconstrained if not contains_direct_repeat(context_l + cds + context_r, repeat_length)]

    assert 0 != len(filtered) == len(constrained) != len(unconstrained)
    assert set(filtered) == set(constrained)


@pytest.mark.parametrize(
    'context_l,aa_seq,context_r,repeat_length',
    [
        # In-frame inverse repeats in CDS
        ('GCTAACGTAC', 'KEEKEFFFFF', 'TCGATGCAGT', 15),
        ('ATCGACCAATTG', 'KEFL', 'GCTAACGTAC', 6),
        ('GCTAACGTAC', 'KEAFLAFL', 'TCGATGCAGT', 6),

        # Out-of-frame repeats in CDS
        ('', 'GYWCKP', '', 4),
        ('', 'VSNAYC', '', 4),
        ('', 'CGVTMK', '', 4),

        # Inverted repeats overlapping contexts
        ('TAATCCT', 'HQTMI', 'CCAGCAT', 4),  # overlaps context_l
        ('TGAGCC', 'PTHGA', 'GCTTTT', 4),  # overlaps context_l
        ('GATACGGG', 'FQYCA', 'TGAA', 4),  # overlaps context_r
    ],
)
def test_inverted_repeat_constraint_enumeration(context_l, aa_seq, context_r, repeat_length):
    view = CodonGraph(aa_seq, context_l=context_l, context_r=context_r).view()
    unconstrained = [*view.enumerate()]

    view.set_constraints([InvertedRepeatConstraint(repeat_length)])
    constrained = [*view.enumerate()]

    filtered = [cds for cds in unconstrained if not contains_inverted_repeat(context_l + cds + context_r, repeat_length)]

    assert 0 != len(filtered) == len(constrained) != len(unconstrained)
    assert set(filtered) == set(constrained)


@pytest.mark.parametrize(
    'context_l,aa_seq,context_r,repeat_length',
    [
        # In-frame repeats in CDS
        ('', 'MMMMMM', '', 9),
        ('', 'MMMMMMMM', '', 12),

        # Out-of-frame repeats in CDS
        ('', 'MMWM', '', 4),
        ('', 'WMMM', '', 4),

        # In contexts
        ('ATCGACCAATTG', 'MIKEY', 'GCTAACGTACATCGACCAATTG', 12),
        ('GCTAAC', 'KEFL', 'TCGATGCAGTGCTAAC', 6),

        # Overlapping contexts
        ('GCTAACG', 'MIKEY', 'TCGATGCAGTGCTAAC', 6),   # overlaps context_l
        ('ATCGACCA', 'KEFL', 'GCTAACATCGACCA', 8),     # overlaps context_r
    ],
)
def test_direct_repeat_constraint_empty_space(context_l, aa_seq, context_r, repeat_length):
    view = CodonGraph(aa_seq, context_l=context_l, context_r=context_r).view()
    view.set_constraints([DirectRepeatConstraint(repeat_length)])

    assert not [*view.enumerate()]


@pytest.mark.parametrize(
    'context_l,aa_seq,context_r,repeat_length',
    [
        # In-frame inverted repeats in CDS
        ('', 'MWH', '', 3),
        ('', 'MHM', '', 3),

        # Out-of-frame inverted repeats in CDS
        ('', 'MHMH', '', 4),
        ('', 'MPWH', '', 4),

        # In contexts
        ('ATCGACCAATTG', 'MIKEY', reverse_complement('ATCGACCAATTG'), 12),
        ('GCTAACGTAC', 'KEFL', reverse_complement('GCTAACGTAC'), 10),

        # Overlapping contexts
        ('GCTAACG', 'MIKEY', 'TCGATGCAGT' + reverse_complement('GCTAAC'), 6),  # overlaps context_l
        ('ATCGACCA', 'KEFL', 'GCTAAC' + reverse_complement('TCGACCA'), 7),  # overlaps context_r
    ],
)
def test_inverted_repeat_constraint_empty_space(context_l, aa_seq, context_r, repeat_length):
    view = CodonGraph(aa_seq, context_l=context_l, context_r=context_r).view()
    view.set_constraints([InvertedRepeatConstraint(repeat_length)])

    assert not [*view.enumerate()]


##########################################
# Sampling tests on larger spaces.
##########################################

@pytest.mark.parametrize(
    'aa_seq,repeat_length,min_distance,max_distance',
    [
        ('MIKEYQTVEGPILAIEWAEGTLPWPMIKEY', 9, 0, None),
        ('MIKEYQTVEGPILAIEWAEGTLPWPMIKEY', 15, 0, None),
        ('MIKEYMIKEYQTLAIEWAEGTLPWPMIKEY', 15, 16, None),
        ('MIKEYQTLAIEWMIKEYAEGTLPWPMIKEY', 15, 16, None),
        ('MIKEYMIKEYMIKEYMIKEY', 15, 0, None),
        ('MIKEYMIKEYMIKEYMIKEYMIKEYQTLAIEWAEGTLPWPMIKEY', 75, 75, None),
    ],
)
def test_direct_repeat_constraint_sampling(aa_seq, repeat_length, min_distance, max_distance):
    view = CodonGraph(aa_seq).view()

    constraint = DirectRepeatConstraint(repeat_length, min_distance=min_distance, max_distance=max_distance)
    view.set_constraints([constraint])

    for _ in range(100):
        cds = view.sample()
        assert not contains_direct_repeat(cds, repeat_length, min_distance=min_distance, max_distance=max_distance)


@pytest.mark.parametrize(
    'aa_seq,repeat_length,min_distance,max_distance',
    [
        ('KEEKEQTVEGPILAIEWAEGTLPWPFFFFF', 9, 0, None),
        ('KEEKEQTVEGPILAIEWAEGTLPWPFFFFF', 15, 0, None),
        ('KEEKEKEEKEQTLAIEWAEGTLPWPFFFFF', 15, 16, None),
        ('KEEKEKEEKEKEEKEFFFFFFFFFFFFFFF', 15, 0, None),
        ('KEKEKEKEKEKEKEKEKEKEKEKEKQTLAIEWAEGTLPWPFFFFFFFFFFFFFFFFFFFFFFFFF', 75, 75, None),
    ],
)
def test_inverted_repeat_constraint_sampling(aa_seq, repeat_length, min_distance, max_distance):
    view = CodonGraph(aa_seq).view()

    constraint = InvertedRepeatConstraint(repeat_length, min_distance=min_distance, max_distance=max_distance)
    view.set_constraints([constraint])

    for _ in range(100):
        cds = view.sample()
        assert not contains_inverted_repeat(cds, repeat_length, min_distance=min_distance, max_distance=max_distance)
