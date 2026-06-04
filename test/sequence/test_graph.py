import pytest

from codeine.sequence.graph import CodonGraph, CodonNode, ContextNode


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
    assert graph.codon_nodes[2].codons == ['AAA']


def test_multiple_codon_restriction_is_applied():
    graph = CodonGraph('MIKEY', codon_restrictions={3: ['AAA', 'AAG']})
    assert graph.codon_nodes[2].codons == ['AAA', 'AAG']


def test_lowercase_sequence_is_accepted():
    graph = CodonGraph("mikey", codon_restrictions={3: "aaa"})
    assert graph.aa_seq == "MIKEY"
    assert graph.codon_restrictions[3] == ["AAA"]


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
    assert graph.codon_nodes[2].pinned_codon == 'AAA'
    assert graph.codon_nodes[2].sample_codon() == 'AAA'


def test_graph_can_unpin_codons():
    graph = CodonGraph('MIKEY')
    graph.pin_codons({3: 'AAA'})
    graph.unpin_codons([3])
    assert graph.codon_nodes[2].pinned_codon is None


def test_graph_can_clear_pins():
    graph = CodonGraph('MIKEY')
    graph.pin_codons({3: 'AAA', 5: 'TAT'})
    graph.clear_pins()
    assert all(node.pinned_codon is None for node in graph.codon_nodes)


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


def test_graph_has_initial_and_final_nodes():
    graph = CodonGraph('MIKEY')
    assert isinstance(graph.initial_node, ContextNode)
    assert isinstance(graph.final_node, ContextNode)


def test_initial_and_final_nodes_store_flanks():
    graph = CodonGraph('MIKEY', flank_l='AAA', flank_r='TTT')
    assert graph.initial_node.sequence == 'AAA'
    assert graph.final_node.sequence == 'TTT'


def test_initial_node_has_no_parents():
    graph = CodonGraph('MIKEY')
    assert graph.initial_node.parents == set()


def test_final_node_has_no_transitions():
    graph = CodonGraph('MIKEY')
    assert graph.final_node.transitions == {}


def test_last_codon_node_points_to_final_node():
    graph = CodonGraph('MIKEY')
    last_node = graph.codon_nodes[-1]
    assert set(last_node.transitions) == set(last_node.codons)
    assert all(target is graph.final_node for target in last_node.transitions.values())
    assert last_node in graph.final_node.parents


def test_codon_nodes_by_pos_excludes_initial_and_final_nodes():
    graph = CodonGraph('MIKEY')
    assert set(graph.codon_nodes_by_pos) == {1, 2, 3, 4, 5}


def test_codon_nodes_excludes_initial_and_final_nodes():
    graph = CodonGraph('MIKEY')
    assert all(isinstance(node, CodonNode) for node in graph.codon_nodes)
    assert graph.initial_node not in graph.codon_nodes
    assert graph.final_node not in graph.codon_nodes


def test_only_two_context_nodes():
    graph = CodonGraph('MIKEY')
    non_codon_nodes = [node for node in graph.nodes if not isinstance(node, CodonNode)]
    assert len(non_codon_nodes) == 2
