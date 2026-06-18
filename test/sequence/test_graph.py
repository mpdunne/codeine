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

def testfind_matching_subpaths_empty_sequence_raises():
    graph = CodonGraph('MIKEY')
    with pytest.raises(ValueError):
        graph.find_matching_subpaths('')


def testfind_matching_subpaths_no_match():
    graph = CodonGraph('MIKEY')

    assert not graph.find_matching_subpaths('AAAAAA')
    assert not graph.find_matching_subpaths('ATGATG')
    assert not graph.find_matching_subpaths('GGGG')


def testfind_matching_subpaths_no_match_long():
    graph = CodonGraph('MIKEY' * 1000)

    assert not graph.find_matching_subpaths('AAAAAA')
    assert not graph.find_matching_subpaths('ATGATG')
    assert not graph.find_matching_subpaths('GGGG')


def testfind_matching_subpaths_sequence_longer_than_possible_returns_empty():
    graph = CodonGraph('MIKEY')
    assert not graph.find_matching_subpaths('A' * 10_000)


def testfind_matching_subpaths_finds_lowercase_sequences():
    graph = CodonGraph('MIKEY')
    assert graph.find_matching_subpaths('atg')


def testfind_matching_subpaths_single_nt():
    graph = CodonGraph('MIKEY')

    matches = graph.find_matching_subpaths('T')
    assert matches
    for path, offset in matches:
        assert len(path) == 1
        node, codon = path[0]
        assert codon[offset] == 'T'


def testfind_matching_subpaths_matches_inside_coding_sequence():
    graph = CodonGraph('MIKEY')

    assert graph.find_matching_subpaths('TAAAAG')
    assert graph.find_matching_subpaths('TAAAAGAG')
    assert graph.find_matching_subpaths('ATGATA')


def testfind_matching_subpaths_matches_inside_coding_sequence_long():
    graph = CodonGraph('MIKEY' * 1000)

    assert graph.find_matching_subpaths('TAAAAG')
    assert graph.find_matching_subpaths('TAAAAGAG')
    assert graph.find_matching_subpaths('ATGATA')


def testfind_matching_subpaths_offset_correct():
    graph = CodonGraph('MIKEY')

    matches = graph.find_matching_subpaths('ATGATAAAGGAATAC')
    assert len(matches) == 1
    path, offset = matches[0]
    assert offset == 0

    matches = graph.find_matching_subpaths('TGATAAAGGAATAC')
    assert len(matches) == 1
    path, offset = matches[0]
    assert offset == 1

    matches = graph.find_matching_subpaths('GATAAAGGAATAC')
    assert len(matches) == 1
    path, offset = matches[0]
    assert offset == 2


def testfind_matching_subpaths_fully_in_contexts():
    graph = CodonGraph('MIKEY', context_l='AAGGAAGGAAGGAAGG')
    assert graph.find_matching_subpaths('GGAAGGAAGG')

    graph = CodonGraph('MIKEY', context_r='AATTAATTAATTAATT')
    assert graph.find_matching_subpaths('AATTAATTAA')


def testfind_matching_subpaths_overlapping_contexts():
    graph = CodonGraph('MIKEY', context_l='AAGGAAGGAAGGAAGG')
    matches = graph.find_matching_subpaths('GGAAGGAAGG' + 'ATG')
    assert matches
    assert len(matches) == 1
    path, offset = matches[0]
    assert isinstance(path[0][0], ContextNode)
    assert isinstance(path[1][0], CodonNode)

    graph = CodonGraph('MIKEY', context_r='AATTAATTAATTAATT')
    matches = graph.find_matching_subpaths('TAC' + 'AATTAATTAA')
    assert matches
    assert len(matches) == 1
    path, offset = matches[0]
    assert isinstance(path[0][0], CodonNode)
    assert isinstance(path[1][0], ContextNode)


def testfind_matching_subpaths_matches_entire_left_context_plus_cds():
    graph = CodonGraph('MIKEY', context_l='AAGG')

    matches = graph.find_matching_subpaths('AAGGATG')
    assert len(matches) == 1

    path, offset = matches[0]
    assert offset == 0
    assert isinstance(path[0][0], ContextNode)
    assert isinstance(path[1][0], CodonNode)


def testfind_matching_subpaths_matches_cds_plus_entire_right_context():
    graph = CodonGraph('MIKEY', context_r='AAGG')

    matches = graph.find_matching_subpaths('TACAAGG')
    assert len(matches) == 1

    path, offset = matches[0]
    assert offset == 0
    assert isinstance(path[-2][0], CodonNode)
    assert isinstance(path[-1][0], ContextNode)


def testfind_matching_subpaths_multiple_matches():
    graph = CodonGraph('MMMM')

    matches = graph.find_matching_subpaths('ATGATG')
    assert len(matches) == 3

    for path, offset in matches:
        assert offset == 0
        assert len(path) == 2
        assert ''.join(codon for node, codon in path) == 'ATGATG'


def testfind_matching_subpaths_single_codon_multiple_matches():
    graph = CodonGraph('KKK')
    matches = graph.find_matching_subpaths('AAA')
    assert len(matches) == 11


def testfind_matching_subpaths_respects_codon_restrictions():
    graph = CodonGraph('MIKEY', codon_restrictions={2: 'ATA'})

    assert graph.find_matching_subpaths('ATGATA')
    assert not graph.find_matching_subpaths('ATGATT')


def testfind_matching_subpaths_ends_inside_codon():
    graph = CodonGraph('MIKEY')

    matches = graph.find_matching_subpaths('ATTAAGG')
    for match in matches:
        path, offset = match
        assert offset == 0
        codons = [codon for node, codon in path]
        assert ''.join(codons).startswith('ATTAAGG')


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


SHORT_AA_SEQUENCES = (
    'M',
    'MIKEY',
    'MILDRED',
    'ELEPHANT',
    'REGINALD',
)

MEDIUM_AA_SEQUENCES = (
    'MIKEY' * 10,
    'MIKEY' * 20,
)

# TODO Improve behaviour for long sequences
LONG_AA_SEQUENCES = (
    'MIKEY' * 100,
#    'MIKEY' * 250,
#    'MIKEY' * 500,
#    'MIKEY' * 1000,
)

CONTEXTS_L = ('', 'aaggaaggaagg')
CONTEXTS_R = ('', 'ttccttccttcc')


@pytest.mark.parametrize(
    'aa_seq',
    SHORT_AA_SEQUENCES + MEDIUM_AA_SEQUENCES + LONG_AA_SEQUENCES,
)
@pytest.mark.parametrize('context_l', CONTEXTS_L)
@pytest.mark.parametrize('context_r', CONTEXTS_R)
def testfind_matching_subpaths_full_sequences(aa_seq, context_l, context_r):
    tt = TranslationTable()

    graph = CodonGraph(aa_seq, context_l=context_l, context_r=context_r, translation_table=tt)

    view = graph.view()
    seqs = [view[i] for i in range(min(view.n_valid_sequences, 10))]

    for seq in seqs:
        matches = graph.find_matching_subpaths(seq)

        # Only one path matches any given full sequence.
        assert len(matches) == 1
        path, offset = matches[0]

        # It should start at the beginning.
        assert offset == 0

        # All codon nodes please :)
        assert all(isinstance(node, CodonNode) for node, codon in path)

        # And the positions should be logical...
        positions = [node.pos for node, codon in path]
        assert positions == [*range(1, 1 + (len(seq) // 3))]

        codons = [codon.upper() for node, codon in path]
        assert seq.upper() == ''.join(codons)
