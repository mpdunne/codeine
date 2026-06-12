import pytest
import random

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


def test_find_matching_subpaths_empty_sequence_raises():
    graph = CodonGraph('MIKEY')
    with pytest.raises(ValueError):
        graph._find_matching_subpaths('')


def test_find_matching_subpaths_no_match():
    graph = CodonGraph('MIKEY')

    assert not graph._find_matching_subpaths('AAAAAA')
    assert not graph._find_matching_subpaths('ATGATG')
    assert not graph._find_matching_subpaths('GGGG')


def test_find_matching_subpaths_no_match_long():
    graph = CodonGraph('MIKEY' * 1000)

    assert not graph._find_matching_subpaths('AAAAAA')
    assert not graph._find_matching_subpaths('ATGATG')
    assert not graph._find_matching_subpaths('GGGG')


def test_find_matching_subpaths_sequence_longer_than_possible_returns_empty():
    graph = CodonGraph('MIKEY')
    assert not graph._find_matching_subpaths('A' * 10_000)


def test_find_matching_subpaths_finds_lowercase_sequences():
    graph = CodonGraph('MIKEY')
    assert graph._find_matching_subpaths('atg')


def test_find_matching_subpaths_single_nt():
    graph = CodonGraph('MIKEY')

    matches = graph._find_matching_subpaths('T')
    assert matches
    for path, offset in matches:
        assert len(path) == 1
        node, codon = path[0]
        assert codon[offset] == 'T'


def test_find_matching_subpaths_matches_inside_coding_sequence():
    graph = CodonGraph('MIKEY')

    assert graph._find_matching_subpaths('TAAAAG')
    assert graph._find_matching_subpaths('TAAAAGAG')
    assert graph._find_matching_subpaths('ATGATA')


def test_find_matching_subpaths_matches_inside_coding_sequence_long():
    graph = CodonGraph('MIKEY' * 1000)

    assert graph._find_matching_subpaths('TAAAAG')
    assert graph._find_matching_subpaths('TAAAAGAG')
    assert graph._find_matching_subpaths('ATGATA')


def test_find_matching_subpaths_offset_correct():
    graph = CodonGraph('MIKEY')

    matches = graph._find_matching_subpaths('ATGATAAAGGAATAC')
    assert len(matches) == 1
    path, offset = matches[0]
    assert offset == 0

    matches = graph._find_matching_subpaths('TGATAAAGGAATAC')
    assert len(matches) == 1
    path, offset = matches[0]
    assert offset == 1

    matches = graph._find_matching_subpaths('GATAAAGGAATAC')
    assert len(matches) == 1
    path, offset = matches[0]
    assert offset == 2


def test_find_matching_subpaths_fully_in_contexts():
    graph = CodonGraph('MIKEY', context_l='AAGGAAGGAAGGAAGG')
    assert graph._find_matching_subpaths('GGAAGGAAGG')

    graph = CodonGraph('MIKEY', context_r='AATTAATTAATTAATT')
    assert graph._find_matching_subpaths('AATTAATTAA')


def test_find_matching_subpaths_overlapping_contexts():
    graph = CodonGraph('MIKEY', context_l='AAGGAAGGAAGGAAGG')
    matches = graph._find_matching_subpaths('GGAAGGAAGG' + 'ATG')
    assert matches
    assert len(matches) == 1
    path, offset = matches[0]
    assert isinstance(path[0][0], ContextNode)
    assert isinstance(path[1][0], CodonNode)

    graph = CodonGraph('MIKEY', context_r='AATTAATTAATTAATT')
    matches = graph._find_matching_subpaths('TAC' + 'AATTAATTAA')
    assert matches
    assert len(matches) == 1
    path, offset = matches[0]
    assert isinstance(path[0][0], CodonNode)
    assert isinstance(path[1][0], ContextNode)


def test_find_matching_subpaths_matches_entire_left_context_plus_cds():
    graph = CodonGraph('MIKEY', context_l='AAGG')

    matches = graph._find_matching_subpaths('AAGGATG')
    assert len(matches) == 1

    path, offset = matches[0]
    assert offset == 0
    assert isinstance(path[0][0], ContextNode)
    assert isinstance(path[1][0], CodonNode)


def test_find_matching_subpaths_matches_cds_plus_entire_right_context():
    graph = CodonGraph('MIKEY', context_r='AAGG')

    matches = graph._find_matching_subpaths('TACAAGG')
    assert len(matches) == 1

    path, offset = matches[0]
    assert offset == 0
    assert isinstance(path[-2][0], CodonNode)
    assert isinstance(path[-1][0], ContextNode)


def test_find_matching_subpaths_multiple_matches():
    graph = CodonGraph('MMMM')

    matches = graph._find_matching_subpaths('ATGATG')
    assert len(matches) == 3

    for path, offset in matches:
        assert offset == 0
        assert len(path) == 2
        assert ''.join(codon for node, codon in path) == 'ATGATG'


def test_find_matching_subpaths_single_codon_multiple_matches():
    graph = CodonGraph('KKK')
    matches = graph._find_matching_subpaths('AAA')
    assert len(matches) == 11


def test_find_matching_subpaths_respects_codon_restrictions():
    graph = CodonGraph('MIKEY', codon_restrictions={2: 'ATA'})

    assert graph._find_matching_subpaths('ATGATA')
    assert not graph._find_matching_subpaths('ATGATT')


def test_find_matching_subpaths_ends_inside_codon():
    graph = CodonGraph('MIKEY')

    matches = graph._find_matching_subpaths('ATTAAGG')
    for match in matches:
        path, offset = match
        assert offset == 0
        codons = [codon for node, codon in path]
        assert ''.join(codons).startswith('ATTAAGG')


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
def test_find_matching_subpaths_full_sequences(aa_seq, context_l, context_r):
    tt = TranslationTable()

    graph = CodonGraph(aa_seq, context_l=context_l, context_r=context_r, translation_table=tt)

    view = graph.view()
    seqs = [view[i] for i in range(min(view.n_valid_sequences, 10))]

    for seq in seqs:
        matches = graph._find_matching_subpaths(seq)

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


def test_banned_sequence_entirely_in_left_context_gives_empty_space():
    graph = CodonGraph('MIKEY', context_l='GAATTC', banned_sequences=['GAATTC'])
    view = graph.view()
    assert view.n_valid_sequences == 0

    graph = CodonGraph('MIKEY', context_l='GAATTC', banned_sequences=['AATTC'])
    view = graph.view()
    assert view.n_valid_sequences == 0

    graph = CodonGraph('MIKEY', context_l='GAATTC', banned_sequences=['GAATT'])
    view = graph.view()
    assert view.n_valid_sequences == 0

    graph = CodonGraph('MIKEY', context_l='GAATTC', banned_sequences=['AATT'])
    view = graph.view()
    assert view.n_valid_sequences == 0


def test_banned_sequence_entirely_in_right_context_gives_empty_space():
    graph = CodonGraph('MIKEY', context_r='GAATTC', banned_sequences=['GAATTC'])
    view = graph.view()
    assert view.n_valid_sequences == 0

    graph = CodonGraph('MIKEY', context_r='GAATTC', banned_sequences=['AATTC'])
    view = graph.view()
    assert view.n_valid_sequences == 0

    graph = CodonGraph('MIKEY', context_r='GAATTC', banned_sequences=['AATT'])
    view = graph.view()
    assert view.n_valid_sequences == 0


def helper_ban_sequences_and_check_comprehensive(
        aa_seq,
        banned_sequences,
        context_l='',
        context_r='',
):
    tt = TranslationTable()
    unconstrained_graph = CodonGraph(
        aa_seq,
        context_l=context_l,
        context_r=context_r,
        translation_table=tt,
    )
    unconstrained_view = unconstrained_graph.view()
    unconstrained_seqs = set(unconstrained_view)

    graph = CodonGraph(
        aa_seq,
        context_l=context_l,
        context_r=context_r,
        translation_table=tt,
        banned_sequences=banned_sequences,
    )
    view = graph.view()
    observed_seqs = set(view)

    for banned_sequence in banned_sequences:
        matches = graph._find_matching_subpaths(banned_sequence)
        assert not matches

    expected_seqs = {seq for seq in unconstrained_seqs
                     if all(banned_sequence.upper() not in seq.upper()
                            for banned_sequence in banned_sequences)}

    assert observed_seqs == expected_seqs
    assert view.n_valid_sequences == len(expected_seqs)


def helper_ban_sequences_and_check_probabilistic(
        aa_seq,
        banned_sequences,
        context_l='',
        context_r='',
        n_samples=100,
):
    tt = TranslationTable()
    unconstrained_graph = CodonGraph(
        aa_seq,
        context_l=context_l,
        context_r=context_r,
        translation_table=tt,
    )
    unconstrained_view = unconstrained_graph.view()

    graph = CodonGraph(
        aa_seq,
        context_l=context_l,
        context_r=context_r,
        translation_table=tt,
        banned_sequences=banned_sequences,
    )
    view = graph.view()

    assert view.n_valid_sequences >= 0

    for banned_sequence in banned_sequences:
        matches = graph._find_matching_subpaths(banned_sequence)
        assert not matches

    for _ in range(n_samples):
        seq = view.sample()
        assert seq in unconstrained_view
        for banned_sequence in banned_sequences:
            assert banned_sequence.upper() not in seq.upper()


def helper_arbitrary_coding_sequence(aa_seq, translation_table):
    return ''.join(translation_table.aa_to_codons[aa][0] for aa in aa_seq)


@pytest.mark.parametrize('aa_seq', SHORT_AA_SEQUENCES)
@pytest.mark.parametrize('context_l', CONTEXTS_L)
@pytest.mark.parametrize('context_r', CONTEXTS_R)
def test_regression_banned_sequences_short_aa_sequence(aa_seq, context_l, context_r):
    banned_seqs = (
        'ATG',
        'TAAAAG',
        'AAGGAA',
        'ATTAAGG',
        'GAATAC',
    )
    helper_ban_sequences_and_check_comprehensive(aa_seq, banned_seqs,
                                                 context_l=context_l, context_r=context_r)


@pytest.mark.parametrize('aa_seq', SHORT_AA_SEQUENCES)
@pytest.mark.parametrize('context_l', CONTEXTS_L)
@pytest.mark.parametrize('context_r', CONTEXTS_R)
def test_regression_banned_sequences_short_aa_sequence_overlapping_banned_sequences(aa_seq, context_l, context_r):
    banned_seqs = (
        'TAAAAG',
        'AAAGGA',
        'AAGGAA',
        'GGAATA',
    )
    helper_ban_sequences_and_check_comprehensive(aa_seq, banned_seqs,
                                                 context_l=context_l, context_r=context_r)


@pytest.mark.parametrize('aa_seq', SHORT_AA_SEQUENCES + MEDIUM_AA_SEQUENCES + LONG_AA_SEQUENCES)
@pytest.mark.parametrize('context_l', CONTEXTS_L)
@pytest.mark.parametrize('context_r', CONTEXTS_R)
def test_regression_banned_whole_sequences(aa_seq, context_l, context_r):
    tt = TranslationTable()
    cds = helper_arbitrary_coding_sequence(aa_seq, tt)

    unconstrained_graph = CodonGraph(
        aa_seq,
        context_l=context_l,
        context_r=context_r,
        translation_table=tt,
    )
    unconstrained_view = unconstrained_graph.view()

    assert cds in unconstrained_view
    unconstrained_n_sequences = unconstrained_view.n_valid_sequences
    assert unconstrained_n_sequences >= 0

    graph = CodonGraph(
        aa_seq,
        context_l=context_l,
        context_r=context_r,
        translation_table=tt,
        banned_sequences=[cds],
    )

    view = graph.view()

    assert cds not in view
    assert view.n_valid_sequences == unconstrained_n_sequences - 1

    # Same thing but with multiple sequences.
    seqs = [unconstrained_view[i] for i in range(min(unconstrained_n_sequences, 5))]

    graph = CodonGraph(
        aa_seq,
        context_l=context_l,
        context_r=context_r,
        translation_table=tt,
        banned_sequences=seqs,
    )

    view = graph.view()

    for seq in seqs:
        assert seq not in view

    assert view.n_valid_sequences == unconstrained_n_sequences - len(seqs)


@pytest.mark.parametrize('aa_seq', SHORT_AA_SEQUENCES + MEDIUM_AA_SEQUENCES + LONG_AA_SEQUENCES)
@pytest.mark.parametrize('context_l', CONTEXTS_L)
@pytest.mark.parametrize('context_r', CONTEXTS_R)
def test_regression_banned_sequences_that_arent_present_anyway(aa_seq, context_l, context_r):
    tt = TranslationTable()

    unconstrained_graph = CodonGraph(
        aa_seq,
        context_l=context_l,
        context_r=context_r,
        translation_table=tt,
    )

    unconstrained_view = unconstrained_graph.view()

    graph = CodonGraph(
        aa_seq,
        context_l=context_l,
        context_r=context_r,
        translation_table=tt,
        banned_sequences=[
            'CCCCCCCCCCCC',
            'GGGGGGGGGGGG',
            'TTTTTTTTTTTT',
        ],
    )

    view = graph.view()
    assert view.n_valid_sequences == unconstrained_view.n_valid_sequences


@pytest.mark.parametrize('aa_seq', MEDIUM_AA_SEQUENCES + LONG_AA_SEQUENCES)
@pytest.mark.parametrize('context_l', CONTEXTS_L)
@pytest.mark.parametrize('context_r', CONTEXTS_R)
def test_regression_banned_sequences_long_aa_sequence(aa_seq, context_l, context_r):
    banned_seqs = (
        'TAAAAG',
        'AAAGGA',
        'AAGGAA',
        'GGAATA',
    )

    helper_ban_sequences_and_check_probabilistic(aa_seq, banned_seqs,
                                                 context_l, context_r, n_samples=1000)


@pytest.mark.parametrize('aa_seq', MEDIUM_AA_SEQUENCES + LONG_AA_SEQUENCES)
@pytest.mark.parametrize('context_l', CONTEXTS_L)
@pytest.mark.parametrize('context_r', CONTEXTS_R)
def test_regression_banned_sequences_long_aa_sequence_overlapping_banned_sequences(aa_seq, context_l, context_r):
    tt = TranslationTable()

    seq = helper_arbitrary_coding_sequence(aa_seq, tt)

    banned_sequences = [
        seq[0:12],
        seq[3:15],
        seq[6:18],
        seq[9:21],
    ]

    helper_ban_sequences_and_check_probabilistic(aa_seq, banned_sequences,
                                                 context_l, context_r, n_samples=1000)


@pytest.mark.parametrize('aa_seq', MEDIUM_AA_SEQUENCES + LONG_AA_SEQUENCES)
@pytest.mark.parametrize('context_l', CONTEXTS_L)
@pytest.mark.parametrize('context_r', CONTEXTS_R)
def test_regression_banned_sequences_long_aa_sequence_nested_banned_sequences(aa_seq, context_l, context_r):
    tt = TranslationTable()

    seq = helper_arbitrary_coding_sequence(aa_seq, tt)

    banned_sequences = (
        seq[0:12],
        seq[3:12],
        seq[5:10],
    )

    helper_ban_sequences_and_check_probabilistic(aa_seq, banned_sequences,
                                                 context_l, context_r, n_samples=1000)


@pytest.mark.parametrize('aa_seq', MEDIUM_AA_SEQUENCES + LONG_AA_SEQUENCES)
@pytest.mark.parametrize('context_l', CONTEXTS_L)
@pytest.mark.parametrize('context_r', CONTEXTS_R)
def test_regression_banned_sequences_long_aa_sequence_long_banned_sequences(aa_seq, context_l, context_r):

    tt = TranslationTable()

    unconstrained_graph = CodonGraph(
        aa_seq,
        context_l=context_l,
        context_r=context_r,
        translation_table=tt,
    )

    unconstrained_view = unconstrained_graph.view()

    # Jenny?????
    rng = random.Random(8675309)

    # Grab 500 random 50nt-long sequences from here.
    # (I personally think that's "long", don't know about you!)
    banned_seqs = []
    for _ in range(500):
        ix = rng.randrange(unconstrained_view.n_valid_sequences)
        seq = unconstrained_view[ix]
        start = rng.randrange(len(seq) - 49)
        banned_seqs.append(seq[start:start + 50])

    helper_ban_sequences_and_check_probabilistic(
        aa_seq,
        banned_seqs,
        context_l=context_l,
        context_r=context_r,
    )


@pytest.mark.parametrize('aa_seq', SHORT_AA_SEQUENCES)
@pytest.mark.parametrize('context_l', CONTEXTS_L)
@pytest.mark.parametrize('context_r', CONTEXTS_R)
def test_regression_banned_sequences_short_aa_sequence_many_banned_sequences(aa_seq, context_l, context_r):
    banned_seqs = (
        'ATG',
        'TAA',
        'AAG',
        'GAA',
        'TAC',
        'ATT',
        'CTG',
        'GGT',
        'GGC',
        'AAA',
        'AAC',
        'GAG',
    )

    helper_ban_sequences_and_check_comprehensive(
        aa_seq,
        banned_seqs,
        context_l=context_l,
        context_r=context_r,
    )


@pytest.mark.parametrize('aa_seq', MEDIUM_AA_SEQUENCES + LONG_AA_SEQUENCES)
@pytest.mark.parametrize('context_l', CONTEXTS_L)
@pytest.mark.parametrize('context_r', CONTEXTS_R)
def test_regression_banned_sequences_long_aa_sequence_many_banned_sequences(aa_seq, context_l, context_r):
    tt = TranslationTable()

    unconstrained_graph = CodonGraph(
        aa_seq,
        context_l=context_l,
        context_r=context_r,
        translation_table=tt,
    )

    unconstrained_view = unconstrained_graph.view()

    rng = random.Random(8675309)

    banned_seqs = []
    for _ in range(500):
        ix = rng.randrange(unconstrained_view.n_valid_sequences)
        seq = unconstrained_view[ix]

        start = rng.randrange(len(seq) - 11)
        banned_seqs.append(seq[start:start + 12])

    helper_ban_sequences_and_check_probabilistic(
        aa_seq,
        banned_seqs,
        context_l=context_l,
        context_r=context_r,
        n_samples=1000,
    )
