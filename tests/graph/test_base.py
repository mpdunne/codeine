import pickle

import pytest

from codeine.graph.base import CodonGraph
from codeine.graph.nodes import ContextNode, CodonNode, EndNode
from codeine.translation.tables import TranslationTable


def test_empty_sequence_raises():
    with pytest.raises(ValueError):
        CodonGraph('')


def test_invalid_amino_acid_raises():
    with pytest.raises(ValueError):
        CodonGraph('MIXEY')


def test_lowercase_sequence_is_accepted():
    graph = CodonGraph('mikey', codon_restrictions={3: 'aaa'})
    assert graph.aa_seq == 'MIKEY'
    assert graph.codon_restrictions[3] == ['AAA']


def test_contexts_are_normalised_to_graph_molecule_type():
    graph = CodonGraph('MIKEY', context_l='aaa', context_r='uuu')
    assert graph.left_context_node.sequence == 'AAA'
    assert graph.right_context_node.sequence == 'TTT'

    tt = TranslationTable(rna=True)
    graph = CodonGraph('MIKEY', context_l='aaa', context_r='ttt', translation_table=tt)
    assert graph.left_context_node.sequence == 'AAA'
    assert graph.right_context_node.sequence == 'UUU'


def test_rna_codon_restrictions_are_normalised():
    tt = TranslationTable(rna=True)
    graph = CodonGraph('MIKEY', codon_restrictions={1: 'ATG'}, translation_table=tt)
    assert graph.codon_restrictions[1] == ['AUG']
    assert graph.codon_node_by_pos(1).codons == ('AUG',)


def test_duplicate_codon_restrictions_are_deduplicated():
    graph = CodonGraph('MIKEY', codon_restrictions={3: ['AAA', 'AAA', 'AAG']})
    assert set(graph.codon_restrictions[3]) == {'AAA', 'AAG'}


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
    assert set(graph.codon_restrictions[3]) == {'AAA', 'AAG'}


def test_single_codon_restriction_is_applied():
    graph = CodonGraph('MIKEY', codon_restrictions={3: 'AAA'})
    node = graph.codon_node_by_pos(3)
    assert node.codons == ('AAA',)


def test_multiple_codon_restriction_is_applied():
    graph = CodonGraph('MIKEY', codon_restrictions={3: ['AAA', 'AAG']})
    node = graph.codon_node_by_pos(3)
    assert set(node.codons) == {'AAA', 'AAG'}


def test_validate_codon_restrictions_respects_existing_restrictions():
    graph = CodonGraph('MIKEY', codon_restrictions={3: ['AAA']})

    with pytest.raises(ValueError):
        graph.validate_codon_restrictions({3: ['AAG']})


def test_graph_has_initial_and_final_nodes():
    graph = CodonGraph('MIKEY')
    assert isinstance(graph.initial_node, ContextNode)
    assert isinstance(graph.final_node, EndNode)


def test_context_nodes_store_contexts():
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
    first_node = graph.codon_node_by_pos(1)

    assert graph.left_context_node.transitions == {'AAA': first_node}
    assert (graph.left_context_node, 'AAA') in first_node.parents


def test_last_codon_node_points_to_right_context_node():
    aa_seq = 'MIKEY'
    graph = CodonGraph(aa_seq)
    last_node = graph.codon_node_by_pos(len(aa_seq))

    assert set(last_node.transitions) == set(last_node.codons)
    assert all(target is graph.right_context_node for target in last_node.transitions.values())
    assert len(graph.right_context_node.parents) == 2
    assert (last_node, 'TAC') in graph.right_context_node.parents
    assert (last_node, 'TAT') in graph.right_context_node.parents


def test_right_context_node_points_to_final_node():
    graph = CodonGraph('MIKEY', context_r='TTT')

    assert graph.right_context_node.transitions == {'TTT': graph.final_node}
    assert (graph.right_context_node, 'TTT') in graph.final_node.parents


def test_codon_nodes_excludes_context_and_final_nodes():
    graph = CodonGraph('MIKEY')
    codon_nodes = graph.codon_nodes
    assert all(isinstance(node, CodonNode) for node in codon_nodes)
    assert {node.pos for node in codon_nodes} == {1, 2, 3, 4, 5}
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


def test_graph_defaults_to_dna():
    graph = CodonGraph('MIKEY')

    assert graph.tt.rna is False


def test_graph_uses_provided_dna_translation_table():
    tt = TranslationTable()

    graph = CodonGraph('MIKEY', translation_table=tt)

    assert graph.tt is tt


def test_graph_uses_provided_rna_translation_table():
    tt = TranslationTable(rna=True)

    graph = CodonGraph('MIKEY', translation_table=tt)

    assert graph.tt is tt


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
