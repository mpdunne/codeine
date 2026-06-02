import pytest

from Bio.Seq import Seq

from codeine.sequence.generate import SequenceGenerator


@pytest.mark.parametrize('aa_seq', ('MIKEY', 'MILDRED', 'STEVEN', 'WILLIAM'))
def test_sg_sequences_translate_correctly(aa_seq):
    sg = SequenceGenerator(aa_seq=aa_seq)
    for _ in range(1000):
        cds = sg.generate()
        translated = Seq(cds).translate()
        assert translated == aa_seq


def test_sg_fixed_codons_are_fixed():
    sg = SequenceGenerator("MIKEY", fixed_codons={2: "ATA"})

    for _ in range(1000):
        cds = sg.generate()
        codons = [cds[i:i + 3] for i in range(0, len(cds), 3)]

        assert codons[1] == "ATA"
        assert str(Seq(cds).translate()) == "MIKEY"


def test_sg_generates_different_sequences():
    sg = SequenceGenerator("MIKEY")
    generated = [sg.generate() for _ in range(1000)]
    assert len(generated) > 1


def test_sg_generates_all_sequences_for_small_case():
    sg = SequenceGenerator("MF")

    generated = {sg.generate() for _ in range(1000)}
    expected = {"ATGTTT", "ATGTTC"}

    assert generated == expected


