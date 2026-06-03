import pytest

from codeine.sequence.graph import CodonGraph


def test_empty_sequence_raises():
    with pytest.raises(ValueError):
        CodonGraph("")


def test_invalid_codon_restriction_positions_raises():
    with pytest.raises(ValueError):
        CodonGraph("MIKEY", codon_restrictions={-1: "ATG"})

    with pytest.raises(ValueError):
        CodonGraph("MIKEY", codon_restrictions={6: "ATG"})


def test_invalid_codon_restriction_value_raises():
    with pytest.raises(ValueError):
        CodonGraph("MIKEY", codon_restrictions={0: "ATT"})

    with pytest.raises(ValueError):
        CodonGraph("MIKEY", codon_restrictions={0: ["ATT"]})


def test_codon_restrictions_are_uppercased():
    graph = CodonGraph("MIKEY", codon_restrictions={3: "aaa"})
    assert graph.codon_restrictions[3] == ["AAA"]

    graph = CodonGraph("MIKEY", codon_restrictions={3: ["aaa"]})
    assert graph.codon_restrictions[3] == ["AAA"]

    graph = CodonGraph('MIKEY', codon_restrictions={3: ['aaa', 'aag']})
    assert graph.codon_restrictions[3] == ['AAA', 'AAG']


def test_single_codon_restriction_is_applied():
    graph = CodonGraph('MIKEY', codon_restrictions={3: 'AAA'})
    assert graph.nodes[2].codons == ['AAA']


def test_multiple_codon_restriction_is_applied():
    graph = CodonGraph('MIKEY', codon_restrictions={3: ['AAA', 'AAG']})
    assert graph.nodes[2].codons == ['AAA', 'AAG']


def test_lowercase_sequence_is_accepted():
    graph = CodonGraph("mikey", codon_restrictions={3: "aaa"})
    assert graph.aa_seq == "MIKEY"
    assert graph.codon_restrictions[3] == ["AAA"]


def test_initial_node_exists():
    graph = CodonGraph("MIKEY")
    assert graph.initial_node is not None
