import pytest

from Bio.Seq import Seq

from codeine.sequence.space import CodingSpace
from codeine.motifs.restriction import RestrictionSite


@pytest.mark.parametrize('aa_seq', ('MIKEY', 'MILDRED', 'STEVEN', 'WILLIAM'))
def test_ss_sequences_translate_correctly(aa_seq):
    space = CodingSpace(aa_seq=aa_seq)
    for _ in range(1000):
        cds = space.sample()
        translated = Seq(cds).translate()
        assert translated == aa_seq


def test_ss_fixed_codons_are_fixed():
    space = CodingSpace('MIKEY', codon_restrictions={2: 'ATA'})

    for _ in range(1000):
        cds = space.sample()
        codons = [cds[i:i + 3] for i in range(0, len(cds), 3)]

        assert codons[1] == 'ATA'
        assert str(Seq(cds).translate()) == 'MIKEY'


def test_ss_generates_different_sequences():
    space = CodingSpace('MIKEY')
    generated = [space.sample() for _ in range(1000)]
    assert len(generated) > 1


def test_ss_generates_all_sequences_for_small_case():
    space = CodingSpace('MF')
    generated = {space.sample() for _ in range(1000)}
    expected = {'ATGTTT', 'ATGTTC'}
    assert generated == expected


def test_ss_pinned_codons_are_fixed():
    space = CodingSpace('MIKEY')
    space.pin_codons({3: 'AAA'})
    for _ in range(100):
        cds = space.sample()
        codons = [cds[i:i + 3] for i in range(0, len(cds), 3)]
        assert codons[2] == 'AAA'


def test_ss_unpin_codons_restores_sampling():
    space = CodingSpace('MIKEY')
    space.pin_codons({3: 'AAA'})
    space.unpin_codons([3])
    sampled = set()
    for _ in range(100):
        cds = space.sample()
        codons = [cds[i:i + 3] for i in range(0, len(cds), 3)]
        sampled.add(codons[2])
    assert sampled == {'AAA', 'AAG'}


def test_ss_clear_pins_restores_sampling():
    space = CodingSpace('MIKEY')

    space.pin_codons({3: 'AAA'})
    space.clear_pins()

    sampled = set()
    for _ in range(100):
        cds = space.sample()
        codons = [cds[i:i + 3] for i in range(0, len(cds), 3)]
        sampled.add(codons[2])

    assert sampled == {'AAA', 'AAG'}


def test_ss_rejects_invalid_pin():
    space = CodingSpace('MIKEY')
    with pytest.raises(ValueError):
        space.pin_codons({3: 'GCT'})


def test_ss_rejects_out_of_range_pin():
    space = CodingSpace('MIKEY')

    with pytest.raises(ValueError):
        space.pin_codons({0: 'ATG'})

    with pytest.raises(ValueError):
        space.pin_codons({6: 'ATG'})


def test_ss_sample_excludes_context_by_default():
    space = CodingSpace(
        aa_seq='MF',
        context_l='AAAA',
        context_r='CCCC',
    )

    cds = space.sample()

    assert cds in {'ATGTTT', 'ATGTTC'}
    assert not cds.startswith('AAAA')
    assert not cds.endswith('CCCC')


def test_mutation_space_raises_if_seq_is_invalid():
    space = CodingSpace('MIKEY')
    
    with pytest.raises(ValueError):
        _ = space.mutants('', [1, 2])

    with pytest.raises(ValueError):
        _ = space.mutants('ATG', [1, 2])

    with pytest.raises(ValueError):
        _ = space.mutants('ATTATTAAAGAATAT', [1, 2])

    with pytest.raises(ValueError):
        _ = space.mutants('ATGATTAAAGAATATATG', [1, 2])


def test_mutation_space_raises_if_positions_are_invalid():
    space = CodingSpace('MIKEY')

    with pytest.raises(ValueError):
        _ = space.mutants('ATGATTAAAGAATATATG', [0])

    with pytest.raises(ValueError):
        _ = space.mutants('ATGATTAAAGAATATATG', [-1])

    with pytest.raises(ValueError):
        _ = space.mutants('ATGATTAAAGAATATATG', [1, 6])


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


def test_base_space_remains_unchanged_after_making_mutation_space():
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


def test_mutation_space_n_valid_sequences():
    aa_seq = 'MIKEY'
    ref_cds = 'ATGATTAAAGAATAT'

    space = CodingSpace(aa_seq)
    assert space.n_valid_sequences == 24

    mut = space.mutants(ref_cds, [1])
    assert mut.n_valid_sequences == 1

    mut = space.mutants(ref_cds, [2])
    assert mut.n_valid_sequences == 3

    mut = space.mutants(ref_cds, [1, 2, 3])
    assert mut.n_valid_sequences == 6

    assert space.n_valid_sequences == 24


def test_sequence_space_getitem():
    space = CodingSpace('MM')
    assert space[0] == 'ATGATG'


def test_view_iter():
    space = CodingSpace('MIKEY')
    seqs = [*space]
    assert len(seqs) == len(set(seqs)) == 24


def test_sequence_space_enumerate():
    space = CodingSpace('F')
    assert list(space.enumerate()) == ['TTT', 'TTC']


def test_sequence_space_contains():
    space = CodingSpace('F')
    assert space.contains('TTT')
    assert space.contains('TTC')
    assert not space.contains('ATG')


def test_sequence_space_mutants_pins_non_mutated_positions():
    space = CodingSpace('FF')
    muts = space.mutants('TTTTTT', positions=[2])
    assert list(muts.enumerate()) == ['TTTTTT', 'TTTTTC']


def test_space_contains():
    space = CodingSpace('MIKEY')
    for _ in range(100):
        seq = space.sample()
        assert seq in space
        assert seq + 'ATG' not in space


def test_expand_forbidden_motifs_none():
    assert CodingSpace._expand_and_validate_forbidden_motifs(None, rna=False) == []


def test_expand_single_string_dna():
    motifs = 'gaattc'
    validated = CodingSpace._expand_and_validate_forbidden_motifs(motifs, rna=False)
    assert validated == ['GAATTC']


def test_expand_single_string_rna():
    motifs = 'GAATTC'
    validated = CodingSpace._expand_and_validate_forbidden_motifs(motifs, rna=True)
    assert validated == ['GAAUUC']


def test_expand_sequence_of_strings_deduplicates_and_sorts():
    motifs = ['tttt', 'UUUU', 'aaaa']
    validated = CodingSpace._expand_and_validate_forbidden_motifs(motifs, rna=False)
    expected = ['AAAA', 'TTTT']
    assert validated == expected


def test_expand_restriction_sites():
    validated = CodingSpace._expand_and_validate_forbidden_motifs([RestrictionSite.EcoRI], rna=False)
    expected = ['GAATTC']
    assert validated == expected

    validated = CodingSpace._expand_and_validate_forbidden_motifs([RestrictionSite.BsaI], rna=False)
    expected = ['GAGACC', 'GGTCTC']
    assert validated == expected

    validated = CodingSpace._expand_and_validate_forbidden_motifs([RestrictionSite.BsaI, 'GGTTCC'], rna=False)
    expected = ['GAGACC', 'GGTCTC', 'GGTTCC']
    assert validated == expected


def test_expand_mixed_forbidden_motifs():
    motifs = [RestrictionSite.EcoRI, 'AAAA']
    validated = CodingSpace._expand_and_validate_forbidden_motifs(motifs, rna=False)
    expected = ['AAAA', 'GAATTC']
    assert validated == expected


def test_empty_forbidden_motif_raises():
    motifs = ''

    with pytest.raises(ValueError, match='Forbidden motifs cannot be empty'):
        CodingSpace._expand_and_validate_forbidden_motifs(motifs, rna=False)


def test_invalid_forbidden_motif_type_raises():
    motifs = [420]
    with pytest.raises(TypeError, match='Forbidden motifs must be strings or codeine.RestrictionSite.'):
        CodingSpace._expand_and_validate_forbidden_motifs(motifs, rna=False)


def test_validate_max_homopolymer_none():
    assert CodingSpace._validate_max_homopolymer(None) is None


def test_validate_max_homopolymer_int():
    assert CodingSpace._validate_max_homopolymer(4) == 4


def test_validate_max_homopolymer_rejects_non_int():
    with pytest.raises(TypeError, match='max_homopolymer must be an integer'):
        CodingSpace._validate_max_homopolymer(4.5)


def test_validate_max_homopolymer_rejects_less_than_one():
    with pytest.raises(ValueError, match='max_homopolymer must be at least 1'):
        CodingSpace._validate_max_homopolymer(0)
