import pytest

from codeine.sequence.graph import CodonGraph, CodonNode



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


def test_node_can_pin_codon():
    node = CodonNode(pos=1, aa='M', codons=['ATG'])
    node.pin_codon('ATG')
    assert node.pinned_codon == 'ATG'
    assert node.sample_codon() == 'ATG'


def test_node_can_unpin_codon():
    node = CodonNode(pos=3, aa='K', codons=['AAA', 'AAG'])
    node.pin_codon('AAA')
    node.unpin_codon()
    assert node.pinned_codon is None
    sampled_codons = {node.sample_codon() for _ in range(100)}
    assert sampled_codons == {'AAA', 'AAG'}


def test_node_pin_codon_uppercases():
    node = CodonNode(pos=3, aa='K', codons=['AAA', 'AAG'])
    node.pin_codon('aaa')
    assert node.pinned_codon == 'AAA'
    assert node.sample_codon() == 'AAA'


def test_node_cannot_pin_invalid_codon():
    node = CodonNode(pos=3, aa='K', codons=['AAA', 'AAG'])
    with pytest.raises(ValueError):
        node.pin_codon('GCT')


def test_graph_can_pin_codons():
    graph = CodonGraph('MIKEY')
    graph.pin_codons({3: 'AAA'})
    assert graph.nodes[2].pinned_codon == 'AAA'
    assert graph.nodes[2].sample_codon() == 'AAA'


def test_graph_can_unpin_codons():
    graph = CodonGraph('MIKEY')
    graph.pin_codons({3: 'AAA'})
    graph.unpin_codons([3])
    assert graph.nodes[2].pinned_codon is None


def test_graph_can_clear_pins():
    graph = CodonGraph('MIKEY')
    graph.pin_codons({3: 'AAA', 5: 'TAT'})
    graph.clear_pins()
    assert all(node.pinned_codon is None for node in graph.nodes)


def test_graph_rejects_out_of_range_pin():
    graph = CodonGraph('MIKEY')

    with pytest.raises(ValueError):
        graph.pin_codons({0: 'ATG'})

    with pytest.raises(ValueError):
        graph.pin_codons({6: 'ATG'})


def test_graph_rejects_pin_outside_codon_restrictions():
    graph = CodonGraph('MIKEY', codon_restrictions={3: ['AAA']})
    with pytest.raises(ValueError):
        graph.pin_codons({3: 'AAG'})