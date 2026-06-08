import pytest

from Bio.Seq import Seq

from codeine.sequence.space import SequenceSpace


@pytest.mark.parametrize('aa_seq', ('MIKEY', 'MILDRED', 'STEVEN', 'WILLIAM'))
def test_ss_sequences_translate_correctly(aa_seq):
    ss = SequenceSpace(aa_seq=aa_seq)
    for _ in range(1000):
        cds = ss.sample()
        translated = Seq(cds).translate()
        assert translated == aa_seq


def test_ss_fixed_codons_are_fixed():
    ss = SequenceSpace("MIKEY", codon_restrictions={2: "ATA"})

    for _ in range(1000):
        cds = ss.sample()
        codons = [cds[i:i + 3] for i in range(0, len(cds), 3)]

        assert codons[1] == "ATA"
        assert str(Seq(cds).translate()) == "MIKEY"


def test_ss_generates_different_sequences():
    ss = SequenceSpace("MIKEY")
    generated = [ss.sample() for _ in range(1000)]
    assert len(generated) > 1


def test_ss_generates_all_sequences_for_small_case():
    ss = SequenceSpace("MF")
    generated = {ss.sample() for _ in range(1000)}
    expected = {"ATGTTT", "ATGTTC"}
    assert generated == expected


def test_ss_pinned_codons_are_fixed():
    ss = SequenceSpace('MIKEY')
    ss.pin_codons({3: 'AAA'})
    for _ in range(100):
        cds = ss.sample()
        codons = [cds[i:i + 3] for i in range(0, len(cds), 3)]
        assert codons[2] == 'AAA'


def test_ss_unpin_codons_restores_sampling():
    ss = SequenceSpace('MIKEY')
    ss.pin_codons({3: 'AAA'})
    ss.unpin_codons([3])
    sampled = set()
    for _ in range(100):
        cds = ss.sample()
        codons = [cds[i:i + 3] for i in range(0, len(cds), 3)]
        sampled.add(codons[2])
    assert sampled == {'AAA', 'AAG'}


def test_ss_clear_pins_restores_sampling():
    ss = SequenceSpace('MIKEY')

    ss.pin_codons({3: 'AAA'})
    ss.clear_pins()

    sampled = set()
    for _ in range(100):
        cds = ss.sample()
        codons = [cds[i:i + 3] for i in range(0, len(cds), 3)]
        sampled.add(codons[2])

    assert sampled == {'AAA', 'AAG'}


def test_ss_rejects_invalid_pin():
    ss = SequenceSpace('MIKEY')
    with pytest.raises(ValueError):
        ss.pin_codons({3: 'GCT'})


def test_ss_rejects_out_of_range_pin():
    ss = SequenceSpace('MIKEY')

    with pytest.raises(ValueError):
        ss.pin_codons({0: 'ATG'})

    with pytest.raises(ValueError):
        ss.pin_codons({6: 'ATG'})

def test_ss_sample_excludes_context_by_default():
    ss = SequenceSpace(
        aa_seq="MF",
        context_l="AAAA",
        context_r="CCCC",
    )

    cds = ss.sample()

    assert cds in {"ATGTTT", "ATGTTC"}
    assert not cds.startswith("AAAA")
    assert not cds.endswith("CCCC")


def test_ss_sample_can_include_context():
    ss = SequenceSpace(
        aa_seq="MF",
        context_l="AAAA",
        context_r="CCCC",
    )

    generated = {ss.sample(include_context=True) for _ in range(1000)}

    assert generated == {
        "AAAAATGTTTCCCC",
        "AAAAATGTTCCCCC",
    }


def test_ss_sample_with_context_still_translates_cds_region():
    ss = SequenceSpace(
        aa_seq="MIKEY",
        context_l="AAAA",
        context_r="CCCC",
    )

    full_seq = ss.sample(include_context=True)
    cds = full_seq[len("AAAA"):-len("CCCC")]

    assert str(Seq(cds).translate()) == "MIKEY"


def test_mutation_space_raises_if_seq_is_invalid():
    ss = SequenceSpace('MIKEY')
    
    with pytest.raises(ValueError):
        _ = ss.mutants('', [1, 2])

    with pytest.raises(ValueError):
        _ = ss.mutants('ATG', [1, 2])

    with pytest.raises(ValueError):
        _ = ss.mutants('ATTATTAAAGAATAT', [1, 2])

    with pytest.raises(ValueError):
        _ = ss.mutants('ATGATTAAAGAATATATG', [1, 2])


def test_mutation_space_raises_if_positions_are_invalid():
    ss = SequenceSpace('MIKEY')

    with pytest.raises(ValueError):
        _ = ss.mutants('ATGATTAAAGAATATATG', [0])

    with pytest.raises(ValueError):
        _ = ss.mutants('ATGATTAAAGAATATATG', [-1])

    with pytest.raises(ValueError):
        _ = ss.mutants('ATGATTAAAGAATATATG', [1, 6])


@pytest.mark.parametrize('aa_seq,positions',
                         (
                                 ('MIKEY', [2, 3]),
                                 ('MILDRED', [2, 3, 4]),
                                 ('STEVEN', [1, 2]),
                                 ('WILLIAM', [2, 5, 7]),
                         ))
def test_mutation_space_mutates_only_specified_positions(aa_seq, positions):
    ss = SequenceSpace(aa_seq)
    ref_cds = ss.sample()

    mut = ss.mutants(ref_cds, positions)
    sampled_seqs = [mut.sample() for _ in range(1000)]
    assert all(ss.contains(s) for s in sampled_seqs)
    assert all(mut.contains(s) for s in sampled_seqs)

    fixed_positions = [pos for pos in range(1, len(aa_seq) + 1) if pos not in positions]

    values_at_fixed_positions = [tuple(seq[(pos - 1) * 3: pos * 3] for pos in fixed_positions) for seq in sampled_seqs]
    values_at_unfixed_positions = [tuple(seq[(pos - 1) * 3: pos * 3] for pos in positions) for seq in sampled_seqs]

    assert len(set(values_at_fixed_positions)) == 1
    assert len(set(values_at_unfixed_positions)) != 1


def test_base_space_remains_unchanged_after_making_mutation_space():
    aa_seq = 'MIKEY'
    ref_cds = 'ATGATTAAAGAATAT'

    ss = SequenceSpace(aa_seq)

    sampled_seqs = [ss.sample() for _ in range(1000)]
    assert not all(s[6:] == ref_cds[6:] for s in sampled_seqs)
    assert not all(s[3:6] == ref_cds[3:6] for s in sampled_seqs)

    mut = ss.mutants(ref_cds, [2])
    sampled_seqs = [mut.sample() for _ in range(1000)]
    assert all(s[6:] == ref_cds[6:] for s in sampled_seqs)
    assert not all(s[3:6] == ref_cds[3:6] for s in sampled_seqs)

    sampled_seqs = [ss.sample() for _ in range(1000)]
    assert not all(s[6:] == ref_cds[6:] for s in sampled_seqs)
    assert not all(s[3:6] == ref_cds[3:6] for s in sampled_seqs)


def test_mutation_space_n_valid_sequences():
    aa_seq = 'MIKEY'
    ref_cds = 'ATGATTAAAGAATAT'

    ss = SequenceSpace(aa_seq)
    assert ss.n_valid_sequences == 24

    mut = ss.mutants(ref_cds, [1])
    assert mut.n_valid_sequences == 1

    mut = ss.mutants(ref_cds, [2])
    assert mut.n_valid_sequences == 3

    mut = ss.mutants(ref_cds, [1, 2, 3])
    assert mut.n_valid_sequences == 6

    assert ss.n_valid_sequences == 24
