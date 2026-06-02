import pytest

from codeine.sequence.graph import CodonGraph


def test_empty_sequence_raises():
    with pytest.raises(ValueError):
        CodonGraph("")


def test_invalid_fixed_codon_position_raises():
    with pytest.raises(ValueError):
        CodonGraph("MIKEY", fixed_codons={-1: "ATG"})

    with pytest.raises(ValueError):
        CodonGraph("MIKEY", fixed_codons={6: "ATG"})


def test_invalid_fixed_codon_identity_raises():
    with pytest.raises(ValueError):
        CodonGraph("MIKEY", fixed_codons={0: "ATT"})


def test_fixed_codons_are_uppercased():
    graph = CodonGraph("MIKEY", fixed_codons={3: "aaa"})

    assert graph.fixed_codons[3] == "AAA"


def test_lowercase_sequence_is_accepted():
    graph = CodonGraph("mikey", fixed_codons={3: "aaa"})

    assert graph.aa_seq == "MIKEY"
    assert graph.fixed_codons[3] == "AAA"


def test_initial_node_exists():
    graph = CodonGraph("MIKEY")

    assert graph.initial_node is not None
