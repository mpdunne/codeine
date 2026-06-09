import pytest

from itertools import product

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
    assert graph.codon_nodes[2].codons == ['AAA']


def test_multiple_codon_restriction_is_applied():
    graph = CodonGraph('MIKEY', codon_restrictions={3: ['AAA', 'AAG']})
    assert graph.codon_nodes[2].codons == ['AAA', 'AAG']


def test_lowercase_sequence_is_accepted():
    graph = CodonGraph('mikey', codon_restrictions={3: 'aaa'})
    assert graph.aa_seq == 'MIKEY'
    assert graph.codon_restrictions[3] == ['AAA']


def test_view_can_pin_codons():
    view = CodonGraph('MIKEY').view()
    view.pin_codons({3: 'AAA'})
    assert view.pinned_codons[3] == ['AAA']


def test_view_can_unpin_codons():
    view = CodonGraph('MIKEY').view()
    view.pin_codons({3: 'AAA'})
    view.unpin_codons([3])
    assert 3 not in view.pinned_codons


def test_view_can_clear_pins():
    view = CodonGraph('MIKEY').view()
    view.pin_codons({3: 'AAA', 5: 'TAT'})
    view.clear_pins()
    assert view.pinned_codons == {}


def test_view_rejects_out_of_range_pin():
    view = CodonGraph('MIKEY').view()

    with pytest.raises(ValueError):
        view.pin_codons({0: 'ATG'})

    with pytest.raises(ValueError):
        view.pin_codons({6: 'ATG'})


def test_view_rejects_pin_outside_codon_restrictions():
    view = CodonGraph('MIKEY', codon_restrictions={3: ['AAA']}).view()
    with pytest.raises(ValueError):
        view.pin_codons({3: 'AAG'})


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
    first_node = graph.codon_nodes[0]

    assert graph.left_context_node.transitions == {'AAA': first_node}
    assert (graph.left_context_node, 'AAA') in first_node.parents


def test_last_codon_node_points_to_right_context_node():
    graph = CodonGraph('MIKEY')
    last_node = graph.codon_nodes[-1]

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


@pytest.fixture
def standard_codon_table():
    return {
        'F': ['TTT', 'TTC'],
        'L': ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'],
        'S': ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'],
        'Y': ['TAT', 'TAC'],
        'C': ['TGT', 'TGC'],
        'W': ['TGG'],
        'P': ['CCT', 'CCC', 'CCA', 'CCG'],
        'H': ['CAT', 'CAC'],
        'Q': ['CAA', 'CAG'],
        'R': ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
        'I': ['ATT', 'ATC', 'ATA'],
        'M': ['ATG'],
        'T': ['ACT', 'ACC', 'ACA', 'ACG'],
        'N': ['AAT', 'AAC'],
        'K': ['AAA', 'AAG'],
        'V': ['GTT', 'GTC', 'GTA', 'GTG'],
        'A': ['GCT', 'GCC', 'GCA', 'GCG'],
        'D': ['GAT', 'GAC'],
        'E': ['GAA', 'GAG'],
        'G': ['GGT', 'GGC', 'GGA', 'GGG']
    }


def helper_enumerate_sequences(aa_seq, aa_to_codons):
    codon_choices = [aa_to_codons[aa] for aa in aa_seq]
    seqs = [''.join(choices) for choices in product(*codon_choices)]
    return seqs


@pytest.mark.parametrize('aa_seq',
                         (
                                 'MIKEY',
                                 'MIKEY',
                                 'M' * 1000,
                                 'SSSSSS',
                                 'M',
                                 'MILDRED',
                                 'ELEPHANT',
                                 'REGINALD',
                         )
                         )
def test_n_valid_sequences_no_restrictions(aa_seq, standard_codon_table):
    view = CodonGraph(aa_seq).view()
    expected_n_all_seqs = len(helper_enumerate_sequences(aa_seq, standard_codon_table))
    assert view.n_valid_sequences == expected_n_all_seqs


def test_n_valid_sequences_fixed_codon(standard_codon_table):
    aa_seq = 'MIKEY'
    codon_restrictions = {2: 'ATC'}
    view = CodonGraph(aa_seq, codon_restrictions=codon_restrictions).view()

    sequences_all = helper_enumerate_sequences(aa_seq, standard_codon_table)
    sequences_restricted = [s for s in sequences_all if s[3:6] == 'ATC']

    assert len(sequences_restricted) != len(sequences_all)
    assert len(sequences_restricted) == view.n_valid_sequences


def test_n_valid_sequences_pinning_and_unpinning(standard_codon_table):
    aa_seq = 'MIKEY'
    view = CodonGraph(aa_seq).view()

    sequences_all = helper_enumerate_sequences(aa_seq, standard_codon_table)
    assert len(sequences_all) == view.n_valid_sequences

    codon_restrictions = {2: 'ATC'}
    sequences_restricted = [s for s in sequences_all if s[3:6] == 'ATC']
    assert len(sequences_restricted) != len(sequences_all)

    view.pin_codons(codon_restrictions)
    assert len(sequences_restricted) == view.n_valid_sequences

    view.clear_pins()
    assert len(sequences_all) == view.n_valid_sequences


@pytest.mark.parametrize('aa_seq',
                         (
                                 'M'
                                 'MIKEY',
                                 'MILDRED',
                                 'ELEPHANT',
                                 'REGINALD',
                         )
                         )
def test_contains_passes_on_valid_sequences(aa_seq, standard_codon_table):
    view = CodonGraph(aa_seq).view()
    expected_all_seqs = helper_enumerate_sequences(aa_seq, standard_codon_table)
    for seq in expected_all_seqs:
        assert view.contains(seq)


def test_contains_fails_on_wrong_length_sequences():
    view = CodonGraph('MIKEY').view()
    assert not view.contains('')
    assert not view.contains('ATG')
    assert not view.contains('ATG' * 10)
    assert not view.contains('ATGA')  # not multiple of 3


@pytest.mark.parametrize(
    'aa_seq, invalid_seq',
    (
        ('M', 'ATT'),
        ('MIKEY', 'ATGATCAAAGAGTAA'),
    ),
)
def test_contains_fails_on_invalid_sequences(aa_seq, invalid_seq):
    view = CodonGraph(aa_seq).view()
    assert not view.contains(invalid_seq)


def test_contains_respects_pinning():
    view = CodonGraph('MS').view()
    assert view.contains('ATGTCT')
    assert view.contains('ATGTCC')

    view.pin_codons({2: 'TCT'})
    assert view.contains('ATGTCT')
    assert not view.contains('ATGTCC')

    view.clear_pins()
    assert view.contains('ATGTCT')
    assert view.contains('ATGTCC')


def test_view_getitem():
    view = CodonGraph('MF').view()
    assert view[0] == 'ATGTTT'
    assert view[1] == 'ATGTTC'


def test_view_len():
    view = CodonGraph('MF').view()
    assert len(view) == 2


def test_view_iter():
    view = CodonGraph('MIKEY').view()
    seqs = [*view]
    assert len(seqs) == len(set(seqs)) == 24


@pytest.mark.parametrize('aa_seq',
                         (
                                 'MIKEY',
                                 'MIKEY',
                                 'M' * 1000,
                                 'SSSSSS',
                                 'M',
                                 'MILDRED',
                                 'ELEPHANT',
                                 'REGINALD',
                         )
                         )
def test_enumerate_sequences(aa_seq, standard_codon_table):
    view = CodonGraph(aa_seq).view()

    generated_all_seqs = [*view.enumerate()]
    expected_all_seqs = helper_enumerate_sequences(aa_seq, standard_codon_table)

    assert view.n_valid_sequences == len(view) == len(generated_all_seqs) == len(expected_all_seqs)
    assert len(generated_all_seqs) == len(expected_all_seqs) == len(set(expected_all_seqs))
    assert set(generated_all_seqs) == set(expected_all_seqs)


def test_enumerate_pinned_sequences(standard_codon_table):
    aa_seq = 'MIKEY'
    view = CodonGraph(aa_seq).view()

    generated_all_seqs = [*view.enumerate()]
    expected_all_seqs = helper_enumerate_sequences(aa_seq, standard_codon_table)

    assert 24 == view.n_valid_sequences == len(view) == len(generated_all_seqs) == len(expected_all_seqs)
    assert 24 == len(generated_all_seqs) == len(expected_all_seqs) == len(set(expected_all_seqs))
    assert set(generated_all_seqs) == set(expected_all_seqs)

    view.pin_codons({2: 'ATC'})
    generated_pinned_seqs = [*view.enumerate()]
    assert 8 == view.n_valid_sequences == len(view)
    assert 8 == len(generated_pinned_seqs)
    assert all(seq[3:6] == 'ATC' for seq in generated_pinned_seqs)

    view.clear_pins()
    generated_unpinned_seqs = [*view.enumerate()]

    assert 24 == view.n_valid_sequences == len(view) == len(generated_unpinned_seqs) == len(expected_all_seqs)
    assert 24 == len(generated_unpinned_seqs) == len(expected_all_seqs) == len(set(expected_all_seqs))
    assert set(generated_unpinned_seqs) == set(expected_all_seqs)


def test_get_works_for_very_large_sequences():
    view = CodonGraph('MIKEY' * 1000).view()
    _ = view[100]
    _ = view[1000000]
    _ = view[10**40]


def test_can_instantiate_with_or_without_codon_weights_or_translation_table():
    tt = TranslationTable()
    graph = CodonGraph('MIKEY', translation_table=tt)
    _ = graph.view()

    tt = TranslationTable(rna=True)
    graph = CodonGraph('MIKEY', translation_table=tt)
    _ = graph.view()

    weights = CodonWeights.ecoli()
    graph = CodonGraph('MIKEY', weights=weights)
    _ = graph.view()

    weights = CodonWeights.ecoli(rna=True)
    graph = CodonGraph('MIKEY', weights=weights)
    _ = graph.view()

    tt = TranslationTable()
    weights = CodonWeights.ecoli()
    graph = CodonGraph('MIKEY', translation_table=tt, weights=weights)
    _ = graph.view()

    tt = TranslationTable(rna=True)
    weights = CodonWeights.ecoli(rna=True)
    graph = CodonGraph('MIKEY', translation_table=tt, weights=weights)
    _ = graph.view()

    tt = TranslationTable()
    weights = CodonWeights.ecoli(rna=True)
    with pytest.raises(ValueError):
        _ = CodonGraph('MIKEY', translation_table=tt, weights=weights)

    tt = TranslationTable(rna=True)
    weights = CodonWeights.ecoli()
    with pytest.raises(ValueError):
        _ = CodonGraph('MIKEY', translation_table=tt, weights=weights)
