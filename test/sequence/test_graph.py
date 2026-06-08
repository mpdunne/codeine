import pytest

from itertools import product

from codeine.sequence.graph import CodonGraph, CodonNode, ContextNode, EndNode


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
        CodonGraph("MIKEY", codon_restrictions={1: "ATT"})

    with pytest.raises(ValueError):
        CodonGraph("MIKEY", codon_restrictions={2: ["TTT"]})


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
    graph = CodonGraph(aa_seq)
    expected_n_all_seqs = len(helper_enumerate_sequences(aa_seq, standard_codon_table))
    assert graph.n_valid_sequences == expected_n_all_seqs


def test_n_valid_sequences_fixed_codon(standard_codon_table):
    aa_seq = 'MIKEY'
    codon_restrictions = {2: 'ATC'}
    graph = CodonGraph(aa_seq, codon_restrictions=codon_restrictions)

    sequences_all = helper_enumerate_sequences(aa_seq, standard_codon_table)
    sequences_restricted = [s for s in sequences_all if s[3:6] == 'ATC']

    assert len(sequences_restricted) != len(sequences_all)
    assert len(sequences_restricted) == graph.n_valid_sequences


def test_n_valid_sequences_pinning_and_unpinning(standard_codon_table):
    aa_seq = 'MIKEY'
    graph = CodonGraph(aa_seq)

    sequences_all = helper_enumerate_sequences(aa_seq, standard_codon_table)
    assert len(sequences_all) == graph.n_valid_sequences

    codon_restrictions = {2: 'ATC'}
    sequences_restricted = [s for s in sequences_all if s[3:6] == 'ATC']
    assert len(sequences_restricted) != len(sequences_all)

    graph.pin_codons(codon_restrictions)
    assert len(sequences_restricted) == graph.n_valid_sequences

    graph.clear_pins()
    assert len(sequences_all) == graph.n_valid_sequences


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
    graph = CodonGraph(aa_seq)
    expected_all_seqs = helper_enumerate_sequences(aa_seq, standard_codon_table)
    for seq in expected_all_seqs:
        assert graph.contains(seq)


def test_contains_fails_on_wrong_length_sequences():
    graph = CodonGraph("MIKEY")
    assert not graph.contains("")
    assert not graph.contains("ATG")
    assert not graph.contains("ATG" * 10)
    assert not graph.contains("ATGA")  # not multiple of 3


@pytest.mark.parametrize(
    "aa_seq, invalid_seq",
    (
        ("M", "ATT"),
        ("MIKEY", "ATGATCAAAGAGTAA"),
    ),
)
def test_contains_fails_on_invalid_sequences(aa_seq, invalid_seq):
    graph = CodonGraph(aa_seq)
    assert not graph.contains(invalid_seq)


def test_contains_respects_pinning():
    graph = CodonGraph("MS")
    assert graph.contains("ATGTCT")
    assert graph.contains("ATGTCC")

    graph.pin_codons({2: "TCT"})
    assert graph.contains("ATGTCT")
    assert not graph.contains("ATGTCC")

    graph.clear_pins()
    assert graph.contains("ATGTCT")
    assert graph.contains("ATGTCC")
