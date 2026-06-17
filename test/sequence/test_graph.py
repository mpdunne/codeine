import pytest
import pickle

from codeine.sequence.graph import CodonGraph, CodonNode, ContextNode, EndNode
from codeine.translation.tables import TranslationTable
from codeine.translation.weights import CodonWeights


def test_empty_sequence_raises():
    with pytest.raises(ValueError):
        CodonGraph('')


def test_invalid_codon_restriction_positions_raises():
    with pytest.raises(ValueError):
        CodonGraph('MIKEY', codon_restrictions={-1: 'ATG'})

    with pytest.raises(ValueError):
        CodonGraph('MIKEY', codon_restrictions={6: 'ATG'})


def test_invalid_codon_restriction_value_raises():
    with pytest.raises(ValueError):
        CodonGraph('MIKEY', codon_restrictions={1: 'ATT'})

    with pytest.raises(ValueError):
        CodonGraph('MIKEY', codon_restrictions={2: ['TTT']})


def test_codon_restrictions_are_uppercased():
    graph = CodonGraph('MIKEY', codon_restrictions={3: 'aaa'})
    assert graph.codon_restrictions[3] == ['AAA']

    graph = CodonGraph('MIKEY', codon_restrictions={3: ['aaa']})
    assert graph.codon_restrictions[3] == ['AAA']

    graph = CodonGraph('MIKEY', codon_restrictions={3: ['aaa', 'aag']})
    assert graph.codon_restrictions[3] == ['AAA', 'AAG']


def test_single_codon_restriction_is_applied():
    graph = CodonGraph('MIKEY', codon_restrictions={3: 'AAA'})
    [node] = graph.codon_nodes_by_pos[3]
    assert node.codons == ['AAA']


def test_multiple_codon_restriction_is_applied():
    graph = CodonGraph('MIKEY', codon_restrictions={3: ['AAA', 'AAG']})
    [node] = graph.codon_nodes_by_pos[3]
    assert node.codons == ['AAA', 'AAG']


def test_lowercase_sequence_is_accepted():
    graph = CodonGraph('mikey', codon_restrictions={3: 'aaa'})
    assert graph.aa_seq == 'MIKEY'
    assert graph.codon_restrictions[3] == ['AAA']


def test_graph_has_initial_and_final_nodes():
    graph = CodonGraph('MIKEY')
    assert isinstance(graph.initial_node, ContextNode)
    assert isinstance(graph.final_node, EndNode)


def test_context_nodes_store_flanks():
    graph = CodonGraph('MIKEY', context_l='AAA', context_r='TTT')
    assert graph.left_context_node.sequence == 'AAA'
    assert graph.right_context_node.sequence == 'TTT'


def test_initial_node_has_no_parents():
    graph = CodonGraph('MIKEY')
    assert graph.initial_node.parents == set()


def test_final_node_has_no_transitions():
    graph = CodonGraph('MIKEY')
    assert graph.final_node.transitions == {}


def test_left_context_node_points_to_first_codon_node():
    graph = CodonGraph('MIKEY', context_l='AAA')
    [first_node] = graph.codon_nodes_by_pos[1]

    assert graph.left_context_node.transitions == {'AAA': first_node}
    assert (graph.left_context_node, 'AAA') in first_node.parents


def test_last_codon_node_points_to_right_context_node():
    aa_seq = 'MIKEY'
    graph = CodonGraph(aa_seq)
    [last_node] = graph.codon_nodes_by_pos[len(aa_seq)]

    assert set(last_node.transitions) == set(last_node.codons)
    assert all(target is graph.right_context_node for target in last_node.transitions.values())
    assert len(graph.right_context_node.parents) == 2
    assert (last_node, 'TAC') in graph.right_context_node.parents
    assert (last_node, 'TAT') in graph.right_context_node.parents


def test_right_context_node_points_to_final_node():
    graph = CodonGraph('MIKEY', context_r='TTT')

    assert graph.right_context_node.transitions == {'TTT': graph.final_node}
    assert (graph.right_context_node, 'TTT') in graph.final_node.parents


def test_codon_nodes_by_pos_excludes_context_and_final_nodes():
    graph = CodonGraph('MIKEY')
    assert set(graph.codon_nodes_by_pos) == {1, 2, 3, 4, 5}


def test_codon_nodes_excludes_context_and_final_nodes():
    graph = CodonGraph('MIKEY')
    assert all(isinstance(node, CodonNode) for node in graph.codon_nodes)
    assert graph.left_context_node not in graph.codon_nodes
    assert graph.right_context_node not in graph.codon_nodes
    assert graph.final_node not in graph.codon_nodes


def test_only_two_context_nodes():
    graph = CodonGraph('MIKEY')
    context_nodes = [node for node in graph.nodes if isinstance(node, ContextNode)]
    assert len(context_nodes) == 2


def test_graph_has_one_end_node():
    graph = CodonGraph('MIKEY')
    end_nodes = [node for node in graph.nodes if isinstance(node, EndNode)]
    assert len(end_nodes) == 1


def test_can_instantiate_with_or_without_codon_weights_or_translation_table():
    tt = TranslationTable()
    graph = CodonGraph('MIKEY', translation_table=tt)
    _ = graph.view()
    assert not graph.tt.rna
    assert not graph.cw.rna

    tt = TranslationTable(rna=True)
    graph = CodonGraph('MIKEY', translation_table=tt)
    _ = graph.view()
    assert graph.tt.rna
    assert graph.cw.rna

    weights = CodonWeights.ecoli()
    graph = CodonGraph('MIKEY', weights=weights)
    _ = graph.view()
    assert not graph.tt.rna
    assert not graph.cw.rna

    weights = CodonWeights.ecoli(rna=True)
    graph = CodonGraph('MIKEY', weights=weights)
    _ = graph.view()
    assert graph.tt.rna
    assert graph.cw.rna

    tt = TranslationTable()
    weights = CodonWeights.ecoli()
    graph = CodonGraph('MIKEY', translation_table=tt, weights=weights)
    _ = graph.view()
    assert not graph.tt.rna
    assert not graph.cw.rna

    tt = TranslationTable(rna=True)
    weights = CodonWeights.ecoli(rna=True)
    graph = CodonGraph('MIKEY', translation_table=tt, weights=weights)
    _ = graph.view()
    assert graph.tt.rna
    assert graph.cw.rna

    tt = TranslationTable()
    weights = CodonWeights.ecoli(rna=True)
    with pytest.raises(ValueError):
        _ = CodonGraph('MIKEY', translation_table=tt, weights=weights)

    tt = TranslationTable(rna=True)
    weights = CodonWeights.ecoli()
    with pytest.raises(ValueError):
        _ = CodonGraph('MIKEY', translation_table=tt, weights=weights)


def test_codon_graph_pickle_preserves_enumeration():
    graph = CodonGraph('MIKEY')
    loaded = pickle.loads(pickle.dumps(graph))
    assert loaded.aa_seq == graph.aa_seq
    assert [*loaded.view().enumerate()] == [*graph.view().enumerate()]


def test_codon_graph_pickle_preserves_constraints():
    graph = CodonGraph('MIKEY', codon_restrictions={2: 'ATC'})
    loaded = pickle.loads(pickle.dumps(graph))
    assert loaded.codon_restrictions == graph.codon_restrictions
    assert loaded.view().n_valid_sequences == graph.view().n_valid_sequences
