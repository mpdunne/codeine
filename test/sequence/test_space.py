import pytest

from Bio.Seq import Seq

from codeine.sequence.space import SequenceSpace


@pytest.mark.parametrize('aa_seq', ('MIKEY', 'MILDRED', 'STEVEN', 'WILLIAM'))
def test_sg_sequences_translate_correctly(aa_seq):
    sg = SequenceSpace(aa_seq=aa_seq)
    for _ in range(1000):
        cds = sg.sample()
        translated = Seq(cds).translate()
        assert translated == aa_seq


def test_sg_fixed_codons_are_fixed():
    sg = SequenceSpace("MIKEY", codon_restrictions={2: "ATA"})

    for _ in range(1000):
        cds = sg.sample()
        codons = [cds[i:i + 3] for i in range(0, len(cds), 3)]

        assert codons[1] == "ATA"
        assert str(Seq(cds).translate()) == "MIKEY"


def test_sg_generates_different_sequences():
    sg = SequenceSpace("MIKEY")
    generated = [sg.sample() for _ in range(1000)]
    assert len(generated) > 1


def test_sg_generates_all_sequences_for_small_case():
    sg = SequenceSpace("MF")
    generated = {sg.sample() for _ in range(1000)}
    expected = {"ATGTTT", "ATGTTC"}
    assert generated == expected


def test_sg_pinned_codons_are_fixed():
    sg = SequenceSpace('MIKEY')
    sg.pin_codons({3: 'AAA'})
    for _ in range(100):
        cds = sg.sample()
        codons = [cds[i:i + 3] for i in range(0, len(cds), 3)]
        assert codons[2] == 'AAA'


def test_sg_unpin_codons_restores_sampling():
    sg = SequenceSpace('MIKEY')
    sg.pin_codons({3: 'AAA'})
    sg.unpin_codons([3])
    sampled = set()
    for _ in range(100):
        cds = sg.sample()
        codons = [cds[i:i + 3] for i in range(0, len(cds), 3)]
        sampled.add(codons[2])
    assert sampled == {'AAA', 'AAG'}


def test_sg_clear_pins_restores_sampling():
    sg = SequenceSpace('MIKEY')

    sg.pin_codons({3: 'AAA'})
    sg.clear_pins()

    sampled = set()
    for _ in range(100):
        cds = sg.sample()
        codons = [cds[i:i + 3] for i in range(0, len(cds), 3)]
        sampled.add(codons[2])

    assert sampled == {'AAA', 'AAG'}


def test_sg_rejects_invalid_pin():
    sg = SequenceSpace('MIKEY')
    with pytest.raises(ValueError):
        sg.pin_codons({3: 'GCT'})


def test_sg_rejects_out_of_range_pin():
    sg = SequenceSpace('MIKEY')

    with pytest.raises(ValueError):
        sg.pin_codons({0: 'ATG'})

    with pytest.raises(ValueError):
        sg.pin_codons({6: 'ATG'})

def test_sg_sample_excludes_context_by_default():
    sg = SequenceSpace(
        aa_seq="MF",
        context_l="AAAA",
        context_r="CCCC",
    )

    cds = sg.sample()

    assert cds in {"ATGTTT", "ATGTTC"}
    assert not cds.startswith("AAAA")
    assert not cds.endswith("CCCC")


def test_sg_sample_can_include_context():
    sg = SequenceSpace(
        aa_seq="MF",
        context_l="AAAA",
        context_r="CCCC",
    )

    generated = {sg.sample(include_context=True) for _ in range(1000)}

    assert generated == {
        "AAAAATGTTTCCCC",
        "AAAAATGTTCCCCC",
    }


def test_sg_sample_with_context_still_translates_cds_region():
    sg = SequenceSpace(
        aa_seq="MIKEY",
        context_l="AAAA",
        context_r="CCCC",
    )

    full_seq = sg.sample(include_context=True)
    cds = full_seq[len("AAAA"):-len("CCCC")]

    assert str(Seq(cds).translate()) == "MIKEY"


def test_mutation_mode_raises_if_seq_is_invalid():
    with pytest.raises(ValueError):
        ss = SequenceSpace('MIKEY')
        ss.enter_mutation_mode('', [1, 2])

    with pytest.raises(ValueError):
        ss = SequenceSpace('MIKEY')
        ss.enter_mutation_mode('ATG', [1, 2])

    with pytest.raises(ValueError):
        ss = SequenceSpace('MIKEY')
        ss.enter_mutation_mode('ATTATTAAAGAATAT', [1, 2])

    with pytest.raises(ValueError):
        ss = SequenceSpace('MIKEY')
        ss.enter_mutation_mode('ATGATTAAAGAATATATG', [1, 2])


@pytest.mark.parametrize('aa_seq,positions',
                         (
                                 ('MIKEY', [2, 3]),
                                 ('MILDRED', [2, 3, 4]),
                                 ('STEVEN', [1, 2]),
                                 ('WILLIAM', [2, 5, 7]),
                         ))
def test_mutation_mode_mutates_only_specified_positions(aa_seq, positions):
    ss = SequenceSpace(aa_seq)
    ref_cds = ss.sample()

    ss.enter_mutation_mode(ref_cds, positions)
    sampled_seqs = [ss.sample() for _ in range(1000)]
    assert all(ss.contains(s) for s in sampled_seqs)

    fixed_positions = [pos for pos in range(1, len(aa_seq) + 1) if pos not in positions]

    values_at_fixed_positions = [tuple(seq[(pos - 1) * 3: pos * 3] for pos in fixed_positions) for seq in sampled_seqs]
    values_at_unfixed_positions = [tuple(seq[(pos - 1) * 3: pos * 3] for pos in positions) for seq in sampled_seqs]

    assert len(set(values_at_fixed_positions)) == 1
    assert len(set(values_at_unfixed_positions)) != 1


def test_can_exit_mutation_mode():
    aa_seq = 'MIKEY'
    ref_cds = 'ATGATTAAAGAATAT'

    ss = SequenceSpace(aa_seq)

    ss.enter_mutation_mode(ref_cds, [2])
    sampled_seqs = [ss.sample() for _ in range(1000)]
    assert all(s[6:] == ref_cds[6:] for s in sampled_seqs)
    assert not all(s[3:6] == ref_cds[3:6] for s in sampled_seqs)

    ss.exit_mutation_mode()
    sampled_seqs = [ss.sample() for _ in range(1000)]
    assert not all(s[6:] == ref_cds[6:] for s in sampled_seqs)
    assert not all(s[3:6] == ref_cds[3:6] for s in sampled_seqs)


def test_mutation_mode_updates_n_valid_sequences():
    aa_seq = 'MIKEY'
    ref_cds = 'ATGATTAAAGAATAT'

    ss = SequenceSpace(aa_seq)
    assert ss.graph.n_valid_sequences == 24

    ss.enter_mutation_mode(ref_cds, [1])
    assert ss.graph.n_valid_sequences == 1
    ss.exit_mutation_mode()

    ss.enter_mutation_mode(ref_cds, [2])
    assert ss.graph.n_valid_sequences == 3

    ss.exit_mutation_mode()
    ss.enter_mutation_mode(ref_cds, [1, 2, 3])
    assert ss.graph.n_valid_sequences == 6

    ss.exit_mutation_mode()
    assert ss.graph.n_valid_sequences == 24


