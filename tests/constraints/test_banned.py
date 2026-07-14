import pytest

from itertools import product

from codeine.graph.base import CodonGraph
from codeine.constraints.base import DEAD_STATE, SAFE_STATE
from codeine.constraints.banned import BannedSequenceConstraint, _find_matching_subpaths


def test_empty_banned_sequence_raises():
    with pytest.raises(ValueError):
        BannedSequenceConstraint([''])


def test_banned_sequence_constraint_finds_ban_inside_left_context():
    graph = CodonGraph('REGINALD', context_l='aaggaaggaagg')
    banned_seqs = (
        'ATG',
        'TAAAAG',
        'AAGGAA',
        'ATTAAGG',
        'GAATAC',
    )
    bsc = BannedSequenceConstraint(banned_seqs)
    bsc.link(graph)
    assert bsc.paths


def helper_find_first_path_for(bsc, sequence):
    """
    Return the first tracked path corresponding to a banned sequence.
    """
    for path in bsc.paths:
        if path.sequence == sequence:
            return path

    raise AssertionError(f'No path found for {sequence!r}')


def helper_walk_path(bsc, path):
    """
    Walk a tracked path from the initial state until it either completes
    a ban or reaches the end of the path.
    """
    state_id = bsc.initial_state_id

    for step in path.steps:
        pos, choice = step
        result = bsc.advance(state_id, pos, choice)

        if result == DEAD_STATE:
            return result

        state_id = result

    return result


def test_banned_sequence_constraint_is_trivial_without_banned_sequences():
    graph = CodonGraph(aa_seq='MIKEY')
    bsc = BannedSequenceConstraint([])
    bsc.link(graph)

    assert bsc.is_trivial
    assert bsc.paths == ()
    assert bsc.starts == {}
    assert bsc.initial_state == 0
    assert bsc.initial_state_id == 0


def test_banned_sequence_constraint_finds_paths_for_possible_banned_sequence():
    graph = CodonGraph(aa_seq='MIKEY')
    bsc = BannedSequenceConstraint(['TCAAA'])
    bsc.link(graph)

    assert not bsc.is_trivial
    assert len(bsc.paths) > 0
    assert all(path.sequence == 'TCAAA' for path in bsc.paths)


def test_banned_sequence_constraint_has_no_paths_for_impossible_banned_sequence():
    graph = CodonGraph(aa_seq='MIKEY')
    bsc = BannedSequenceConstraint(['CCCCCC'])
    bsc.link(graph)

    assert bsc.is_trivial
    assert bsc.paths == ()
    assert bsc.starts == {}


def test_starts_are_built_from_first_path_step():
    graph = CodonGraph(aa_seq='MIKEY')
    bsc = BannedSequenceConstraint(['TCAAA'])
    bsc.link(graph)
    path = helper_find_first_path_for(bsc, 'TCAAA')
    first_step = path.steps[0]

    assert first_step in bsc.starts
    assert bsc.starts[first_step]


def test_safe_choice_returns_empty_state():
    graph = CodonGraph(aa_seq='MIKEY')
    bsc = BannedSequenceConstraint(['TCAAA'])
    bsc.link(graph)
    path = helper_find_first_path_for(bsc, 'TCAAA')

    pos, _choice = path.steps[0]
    result = bsc.advance(bsc.initial_state_id, pos, 'ATG')

    assert result == 0


def test_choice_can_start_watch():
    graph = CodonGraph(aa_seq='MIKEY')
    bsc = BannedSequenceConstraint(['TCAAA'])
    bsc.link(graph)
    path = helper_find_first_path_for(bsc, 'TCAAA')

    pos, choice = path.steps[0]
    result = bsc.advance(bsc.initial_state_id, pos, choice)

    assert result != DEAD_STATE
    assert result != 0


def test_choice_can_immediately_complete_banned_sequence():
    graph = CodonGraph(aa_seq='MIKEY')
    bsc = BannedSequenceConstraint(['ATG'])
    bsc.link(graph)
    path = helper_find_first_path_for(bsc, 'ATG')

    pos, choice = path.steps[0]
    result = bsc.advance(bsc.initial_state_id, pos, choice)

    assert result == DEAD_STATE


def test_existing_watch_can_complete_banned_sequence():
    graph = CodonGraph(aa_seq='MIKEY')
    bsc = BannedSequenceConstraint(['TCAAA'])
    bsc.link(graph)
    path = helper_find_first_path_for(bsc, 'TCAAA')

    pos_1, choice_1 = path.steps[0]
    result_1 = bsc.advance(bsc.initial_state_id, pos_1, choice_1)
    pos_2, choice_2 = path.steps[1]
    result_2 = bsc.advance(result_1, pos_2, choice_2)

    assert result_1 != DEAD_STATE
    assert result_2 == DEAD_STATE


def test_existing_watch_drops_if_choice_does_not_match():
    graph = CodonGraph(aa_seq='MIKEY')
    bsc = BannedSequenceConstraint(['TCAAA'])
    bsc.link(graph)
    path = helper_find_first_path_for(bsc, 'TCAAA')

    pos_2, _choice_2 = path.steps[1]

    pos_1, choice_1 = path.steps[0]
    result_1 = bsc.advance(bsc.initial_state_id, pos_1, choice_1)
    result_2 = bsc.advance(result_1, pos_2, 'GAG')

    assert result_1 != DEAD_STATE
    assert result_1 != 0
    assert result_2 == 0


def test_multiple_watches_can_be_active():
    graph = CodonGraph(aa_seq='MIKEY')
    bsc = BannedSequenceConstraint(['TCAAA', 'TCAAG'])
    bsc.link(graph)

    path = helper_find_first_path_for(bsc, 'TCAAA')

    pos, choice = path.steps[0]
    result = bsc.advance(bsc.initial_state_id, pos, choice)

    assert result != DEAD_STATE
    assert len(bsc.states[result]) >= 2


def test_one_of_multiple_watches_can_complete_ban():
    graph = CodonGraph(aa_seq='MIKEY')
    bsc = BannedSequenceConstraint(['TCAAA', 'TCAAG'])
    bsc.link(graph)

    path = helper_find_first_path_for(bsc, 'TCAAA')

    pos_1, choice_1 = path.steps[0]
    result_1 = bsc.advance(bsc.initial_state_id, pos_1, choice_1)
    pos_2, choice_2 = path.steps[1]
    result_2 = bsc.advance(result_1, pos_2, choice_2)

    assert len(bsc.states[result_1]) >= 2
    assert result_2 == DEAD_STATE


def test_state_is_a_frozenset():
    graph = CodonGraph(aa_seq='MIKEY')
    bsc = BannedSequenceConstraint(['TCAAA'])
    bsc.link(graph)
    path = helper_find_first_path_for(bsc, 'TCAAA')

    pos, choice = path.steps[0]
    result = bsc.advance(bsc.initial_state_id, pos, choice)

    assert isinstance(bsc.states[result], frozenset)


def test_banned_sequence_constraint_finds_banned_sequence_crossing_left_context():
    graph = CodonGraph(aa_seq='MIKEY', context_l='TCA')
    bsc = BannedSequenceConstraint(['TCAATG'])
    bsc.link(graph)
    assert not bsc.is_trivial


def test_left_context_can_immediately_complete_ban():
    graph = CodonGraph(aa_seq='MIKEY', context_l='TCA')
    bsc = BannedSequenceConstraint(['TCAATG'])
    bsc.link(graph)
    path = helper_find_first_path_for(bsc, 'TCAATG')
    result = helper_walk_path(bsc, path)

    assert result == DEAD_STATE


def test_banned_sequence_constraint_finds_banned_sequence_crossing_right_context():
    graph = CodonGraph(aa_seq='MIKEY', context_r='AAA')
    bsc = BannedSequenceConstraint(['ATACAAA'])
    bsc.link(graph)
    assert not bsc.is_trivial

    graph = CodonGraph(aa_seq='MIKEY', context_r='AAA')
    bsc = BannedSequenceConstraint(['GGGAAA'])
    bsc.link(graph)
    assert bsc.is_trivial


def test_existing_watch_can_complete_in_right_context():
    graph = CodonGraph(aa_seq='MIKEY', context_r='AAA')
    bsc = BannedSequenceConstraint(['ATACAAA'])
    bsc.link(graph)
    path = helper_find_first_path_for(bsc, 'ATACAAA')
    result = helper_walk_path(bsc, path)

    assert result == DEAD_STATE


@pytest.mark.parametrize(
    'context_l,context_r,sequence',
    [
        ('', '', 'A'),
        ('', '', 'AT'),
        ('', '', 'ATG'),
        ('', '', 'TG'),
        ('', '', 'ATTAAG'),
        ('', '', 'ATTAAGG'),
        ('', '', 'TTAAG'),
        ('', '', 'GATA'),
        ('AAGG', '', 'GGAT'),
        ('AAGGTT', '', 'GG'),
        ('AAGGTT', '', 'GGTTATG'),
        ('', 'AAGG', 'ACAAG'),
        ('', 'AAGGTT', 'TACAAG'),
    ],
)
def test_found_paths_are_walkable(context_l, context_r, sequence):
    graph = CodonGraph('MIKEY', context_l=context_l, context_r=context_r)
    bsc = BannedSequenceConstraint([sequence])
    bsc.link(graph)
    path = helper_find_first_path_for(bsc, sequence)

    assert helper_walk_path(bsc, path) == DEAD_STATE


@pytest.mark.parametrize(
    'sequence,expected_offset',
    [
        ('ATG', 0),
        ('TG', 1),
        ('G', 2),
    ],
)
def test_offsets_are_correct(sequence, expected_offset):
    graph = CodonGraph('MIKEY')
    bsc = BannedSequenceConstraint([sequence])
    bsc.link(graph)
    path = helper_find_first_path_for(bsc, sequence)

    assert path.offset == expected_offset


def test_watch_survives_until_partial_final_codon_match():
    graph = CodonGraph('MIKEY')
    bsc = BannedSequenceConstraint(['ATTAAG'])
    bsc.link(graph)
    path = helper_find_first_path_for(bsc, 'ATTAAG')

    state_id = bsc.initial_state_id

    for step in path.steps[:-1]:
        pos, choice = step
        result = bsc.advance(state_id, pos, choice)
        assert result != DEAD_STATE
        state_id = result

    pos, choice = path.steps[-1]
    result = bsc.advance(state_id, pos, choice)

    assert result == DEAD_STATE


def test_ban_longer_than_choice_keeps_watch_alive():
    graph = CodonGraph('MIKEY')
    bsc = BannedSequenceConstraint(['ATGATA'])
    bsc.link(graph)
    path = helper_find_first_path_for(bsc, 'ATGATA')

    pos, choice = path.steps[0]
    result = bsc.advance(bsc.initial_state_id, pos, choice)

    assert result != DEAD_STATE
    assert result

    state = bsc.states[result]
    path_ix, matched_length = next(iter(state))
    assert bsc.paths[path_ix] == path
    assert matched_length == 3


def test_duplicate_banned_sequences_do_not_break_tracking():
    graph = CodonGraph('MIKEY')
    bsc = BannedSequenceConstraint(['TCAAA', 'TCAAA'])
    bsc.link(graph)
    path = helper_find_first_path_for(bsc, 'TCAAA')

    result = helper_walk_path(bsc, path)
    assert result == DEAD_STATE


def test_different_bans_can_start_from_same_choice():
    graph = CodonGraph('MIKEY')
    bsc = BannedSequenceConstraint(['ATGA', 'ATGAT'])
    bsc.link(graph)
    path = helper_find_first_path_for(bsc, 'ATGA')

    pos, choice = path.steps[0]
    result = bsc.advance(bsc.initial_state_id, pos, choice)

    assert result != DEAD_STATE
    assert len(bsc.states[result]) >= 2


def test_shorter_ban_wins_when_multiple_bans_share_prefix():
    graph = CodonGraph('MIKEY')
    bsc = BannedSequenceConstraint(['ATGA', 'ATGATA'])
    bsc.link(graph)
    path = helper_find_first_path_for(bsc, 'ATGA')

    result = helper_walk_path(bsc, path)
    assert result == DEAD_STATE


def test_longer_ban_can_complete_after_shorter_related_ban_if_shorter_absent():
    graph = CodonGraph('MIKEY')
    bsc = BannedSequenceConstraint(['ATGATA'])
    bsc.link(graph)
    path = helper_find_first_path_for(bsc, 'ATGATA')

    result = helper_walk_path(bsc, path)
    assert result == DEAD_STATE


def test_unrelated_active_watch_does_not_prevent_new_watch_starting():
    graph = CodonGraph('MIKEY')
    bsc = BannedSequenceConstraint(['TCAAA', 'GAA'])
    bsc.link(graph)
    path = helper_find_first_path_for(bsc, 'TCAAA')

    pos_1, choice_1 = path.steps[0]
    result_1 = bsc.advance(bsc.initial_state_id, pos_1, choice_1)
    pos_2, choice_2 = path.steps[1]
    result_2 = bsc.advance(result_1, pos_2, choice_2)

    assert result_2 == DEAD_STATE


def test_path_steps_are_never_empty():
    graph = CodonGraph('MIKEY')
    bsc = BannedSequenceConstraint(['A', 'AT', 'ATG', 'TCAAA'])
    bsc.link(graph)

    assert bsc.paths
    assert all(path.steps for path in bsc.paths)


def test_starts_only_reference_real_paths():
    graph = CodonGraph('MIKEY')
    bsc = BannedSequenceConstraint(['A', 'ATG', 'TCAAA'])
    bsc.link(graph)

    for starts in bsc.starts.values():
        for watch in starts:
            if watch is None:
                continue

            path_ix, matched_length = watch

            assert 0 <= path_ix < len(bsc.paths)
            assert matched_length > 0


def test_all_start_keys_are_real_first_steps():
    graph = CodonGraph('MIKEY')
    bsc = BannedSequenceConstraint(['A', 'ATG', 'TCAAA'])
    bsc.link(graph)

    first_steps = {path.steps[0] for path in bsc.paths}
    assert set(bsc.starts) <= first_steps


def test_every_found_path_really_contains_banned_sequence():
    graph = CodonGraph('MIKEY', context_l='AAGG', context_r='TTCC')
    bsc = BannedSequenceConstraint(['GGATG', 'TACAAG', 'ATTAAG'])
    bsc.link(graph)

    for path in bsc.paths:
        emitted = ''.join(choice for pos, choice in path.steps)
        visible = emitted[path.offset:]
        assert visible.startswith(path.sequence)


def test_walking_every_found_path_completes_ban():
    graph = CodonGraph('MIKEY', context_l='AAGG', context_r='TTCC')
    bsc = BannedSequenceConstraint(['GGATG', 'TACAAG', 'ATTAAG'])
    bsc.link(graph)

    for path in bsc.paths:
        result = helper_walk_path(bsc, path)
        assert result == DEAD_STATE


def test_safe_walk_drops_all_active_watches():
    graph = CodonGraph('MIKEY')
    bsc = BannedSequenceConstraint(['TCAAA'])
    bsc.link(graph)
    path = helper_find_first_path_for(bsc, 'TCAAA')

    pos_2, _choice_2 = path.steps[1]

    pos_1, choice_1 = path.steps[0]
    result_1 = bsc.advance(bsc.initial_state_id, pos_1, choice_1)
    result_2 = bsc.advance(result_1, pos_2, 'GAG')

    assert result_1 != 0
    assert result_2 == 0


def test_every_bsc_path_is_walkable_from_initial_state():
    graph = CodonGraph('MIKEY', context_l='AAGGTT', context_r='CCAAGG')
    banned = ['A', 'AT', 'ATG', 'TGATA', 'GGTTATG', 'TACCCA']
    bsc = BannedSequenceConstraint(banned)
    bsc.link(graph)

    for path in bsc.paths:
        result = helper_walk_path(bsc, path)
        assert result == DEAD_STATE


def test_find_matching_subpaths_empty_sequence_raises():
    graph = CodonGraph('MIKEY')
    with pytest.raises(ValueError):
        _find_matching_subpaths(graph, '')


def test_find_matching_subpaths_no_match():
    graph = CodonGraph('MIKEY')

    assert not _find_matching_subpaths(graph, 'AAAAAA')
    assert not _find_matching_subpaths(graph, 'ATGATG')
    assert not _find_matching_subpaths(graph, 'GGGG')


def test_find_matching_subpaths_no_match_long():
    graph = CodonGraph('MIKEY' * 1000)

    assert not _find_matching_subpaths(graph, 'AAAAAA')
    assert not _find_matching_subpaths(graph, 'ATGATG')
    assert not _find_matching_subpaths(graph, 'GGGG')


def test_find_matching_subpaths_sequence_longer_than_possible_returns_empty():
    graph = CodonGraph('MIKEY')
    assert not _find_matching_subpaths(graph, 'A' * 10_000)


def test_find_matching_subpaths_finds_lowercase_sequences():
    graph = CodonGraph('MIKEY')
    assert _find_matching_subpaths(graph, 'atg')


def test_find_matching_subpaths_single_nt():
    graph = CodonGraph('MIKEY')

    matches = _find_matching_subpaths(graph, 'T')
    assert matches
    for match in matches:
        path = match.steps
        offset = match.offset
        assert len(path) == 1
        _pos, codon = path[0]
        assert codon[offset] == 'T'


def test_find_matching_subpaths_matches_inside_coding_sequence():
    graph = CodonGraph('MIKEY')

    assert _find_matching_subpaths(graph, 'TAAAAG')
    assert _find_matching_subpaths(graph, 'TAAAAGAG')
    assert _find_matching_subpaths(graph, 'ATGATA')


def test_find_matching_subpaths_matches_inside_coding_sequence_long():
    graph = CodonGraph('MIKEY' * 1000)

    assert _find_matching_subpaths(graph, 'TAAAAG')
    assert _find_matching_subpaths(graph, 'TAAAAGAG')
    assert _find_matching_subpaths(graph, 'ATGATA')


def test_find_matching_subpaths_offset_correct():
    graph = CodonGraph('MIKEY')

    matches = _find_matching_subpaths(graph, 'ATGATAAAGGAATAC')
    assert len(matches) == 1
    assert matches[0].offset == 0

    matches = _find_matching_subpaths(graph, 'TGATAAAGGAATAC')
    assert len(matches) == 1
    assert matches[0].offset == 1

    matches = _find_matching_subpaths(graph, 'GATAAAGGAATAC')
    assert len(matches) == 1
    assert matches[0].offset == 2


def test_find_matching_subpaths_fully_in_contexts():
    graph = CodonGraph('MIKEY', context_l='AAGGAAGGAAGGAAGG')
    assert _find_matching_subpaths(graph, 'GGAAGGAAGG')

    graph = CodonGraph('MIKEY', context_r='AATTAATTAATTAATT')
    assert _find_matching_subpaths(graph, 'AATTAATTAA')


def test_find_matching_subpaths_overlapping_contexts():
    graph = CodonGraph('MIKEY', context_l='AAGGAAGGAAGGAAGG')
    matches = _find_matching_subpaths(graph, 'GGAAGGAAGG' + 'ATG')
    assert matches
    assert len(matches) == 1
    match = matches[0]
    assert match.offset == 6
    assert ''.join(choice for _pos, choice in match.steps)[match.offset:] \
        .startswith(match.sequence)

    graph = CodonGraph('MIKEY', context_r='AATTAATTAATTAATT')
    matches = _find_matching_subpaths(graph, 'TAC' + 'AATTAATTAA')
    assert matches
    assert len(matches) == 1
    match = matches[0]
    assert match.offset == 0
    assert ''.join(choice for _pos, choice in match.steps)[match.offset:] \
        .startswith(match.sequence)


def test_find_matching_subpaths_matches_entire_left_context_plus_cds():
    graph = CodonGraph('MIKEY', context_l='AAGG')

    matches = _find_matching_subpaths(graph, 'AAGGATG')
    assert len(matches) == 1

    match = matches[0]
    assert match.offset == 0
    assert ''.join(choice for _pos, choice in match.steps).startswith('AAGGATG')


def test_find_matching_subpaths_matches_cds_plus_entire_right_context():
    graph = CodonGraph('MIKEY', context_r='AAGG')

    matches = _find_matching_subpaths(graph, 'TACAAGG')
    assert len(matches) == 1

    match = matches[0]
    assert match.offset == 0
    assert ''.join(choice for _pos, choice in match.steps).startswith('TACAAGG')


def test_find_matching_subpaths_multiple_matches():
    graph = CodonGraph('MMMM')

    matches = _find_matching_subpaths(graph, 'ATGATG')
    assert len(matches) == 3

    for match in matches:
        assert match.offset == 0
        assert len(match.steps) == 2
        assert ''.join(codon for _pos, codon in match.steps) == 'ATGATG'


def test_find_matching_subpaths_single_codon_multiple_matches():
    graph = CodonGraph('KKK')
    matches = _find_matching_subpaths(graph, 'AAA')
    assert len(matches) == 11


def test_find_matching_subpaths_respects_codon_restrictions():
    graph = CodonGraph('MIKEY', codon_restrictions={2: 'ATA'})

    assert _find_matching_subpaths(graph, 'ATGATA')
    assert not _find_matching_subpaths(graph, 'ATGATT')


def test_find_matching_subpaths_ends_inside_codon():
    graph = CodonGraph('MIKEY')

    matches = _find_matching_subpaths(graph, 'ATTAAGG')
    for match in matches:
        assert match.offset == 0
        codons = [codon for _pos, codon in match.steps]
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

LONG_AA_SEQUENCES = (
    'MIKEY' * 100,
    'MIKEY' * 250,
    'MIKEY' * 500,
)

CONTEXTS_L = ('', 'aaggaaggaagg')
CONTEXTS_R = ('', 'ttccttccttcc')


standard_codon_table = {
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


def helper_yield_sequences(aa_seq):
    codon_choices = [standard_codon_table[aa] for aa in aa_seq]
    for choices in product(*codon_choices):
         yield ''.join(choices)


@pytest.mark.parametrize(
    'aa_seq',
    SHORT_AA_SEQUENCES + MEDIUM_AA_SEQUENCES + LONG_AA_SEQUENCES,
)
@pytest.mark.parametrize('context_l', CONTEXTS_L)
@pytest.mark.parametrize('context_r', CONTEXTS_R)
def test_find_matching_subpaths_full_sequences(aa_seq, context_l, context_r):

    # Get 10 sequences the dumb way.
    gen = helper_yield_sequences(aa_seq)
    seqs = []
    for _ in range(10):
        try:
            seqs.append(next(gen))
        except StopIteration:
            break

    graph = CodonGraph(aa_seq, context_l=context_l, context_r=context_r )

    for seq in seqs:
        matches = _find_matching_subpaths(graph, seq)

        # Only one path matches any given full sequence.
        assert len(matches) == 1
        match = matches[0]

        # It should start at the beginning.
        assert match.offset == 0

        # And the positions should be logical...
        positions = [pos for pos, _codon in match.steps]
        assert positions == [*range(1, 1 + (len(seq) // 3))]

        codons = [codon.upper() for _pos, codon in match.steps]
        assert seq.upper() == ''.join(codons)


@pytest.mark.parametrize(
    'banned_sequence',
    (
        'GAATTC',
        'AATTC',
        'GAATT',
    ),
)
def test_banned_sequence_entirely_in_left_context_is_dead(banned_sequence):
    graph = CodonGraph('MIKEY', context_l='GAATTC')
    bsc = BannedSequenceConstraint([banned_sequence])
    bsc.link(graph)
    state = bsc.advance(bsc.initial_state, 0, graph.left_context_node.sequence)
    assert state == DEAD_STATE


@pytest.mark.parametrize(
    'banned_sequence',
    (
        'GAATTC',
        'AATTC',
        'GAATT',
    ),
)
def test_banned_sequence_entirely_in_right_context_gives_empty_space(banned_sequence):
    graph = CodonGraph('MIKEY', context_r='GAATTC')
    bsc = BannedSequenceConstraint([banned_sequence])
    bsc.link(graph)
    state = bsc.advance(bsc.initial_state, graph.right_context_node.pos, graph.right_context_node.sequence)
    assert state == DEAD_STATE


def test_banned_sequence_spanning_left_context_and_first_codon_is_dead():
    graph = CodonGraph('ELEPHANT', context_l='AAGGATGATG')
    constraint = BannedSequenceConstraint(['AAGGATGATGGAA'])
    constraint.link(graph)

    state = constraint.advance(constraint.initial_state, graph.left_context_node.pos, graph.left_context_node.sequence)
    state = constraint.advance(state, graph.codon_nodes[0].pos, 'GAA')

    assert state == DEAD_STATE


def test_banned_sequence_spanning_last_codon_and_right_context_is_dead():
    graph = CodonGraph('ELEPHANT', context_r='AAGGATGATG')
    constraint = BannedSequenceConstraint(['CGAAGGATGATG'])
    constraint.link(graph)

    state = constraint.advance(constraint.initial_state, graph.codon_nodes[-1].pos, 'ACG')
    state = constraint.advance(state, graph.right_context_node.pos, graph.right_context_node.sequence)

    assert state == DEAD_STATE


def test_banned_sequence_spanning_both_contexts_is_dead():
    coding_sequence = 'GAGCTTGAGCCGCATGCCAATACG'

    graph = CodonGraph('ELEPHANT', context_l='TTAA', context_r='AAGG')
    constraint = BannedSequenceConstraint(['AA' + coding_sequence + 'AA'])
    constraint.link(graph)

    state = constraint.advance(
        constraint.initial_state,
        graph.left_context_node.pos,
        graph.left_context_node.sequence
    )

    for node, i in zip(graph.codon_nodes, range(0, len(coding_sequence), 3)):
        state = constraint.advance(state, node.pos, coding_sequence[i:i + 3])

    state = constraint.advance(
        state,
        graph.right_context_node.pos,
        graph.right_context_node.sequence
    )

    assert state == DEAD_STATE


def test_banned_sequence_constraint_relink_resets_internal_state():
    constraint = BannedSequenceConstraint(['AAA'])

    graph1 = CodonGraph('KKK')
    graph2 = CodonGraph('MMMM')

    constraint.link(graph1)

    # Populate the registry/cache.
    constraint.advance(constraint.initial_state, 1, next(iter(graph1.initial_node.transitions)))

    assert len(constraint.states) > 1 or constraint.advance_cache

    constraint.link(graph2)

    assert constraint.states == [frozenset()]
    assert constraint.state_ids == {frozenset(): 0}
    assert constraint.advance_cache == {}


def test_banned_sequence_constraint_is_trivial_before_linking():
    constraint = BannedSequenceConstraint(['AAA'])
    assert constraint.is_trivial


def test_banned_sequence_constraint_preserves_terminal_states():
    constraint = BannedSequenceConstraint(['AAA'])

    assert constraint.advance(DEAD_STATE, 1, 'AAA') == DEAD_STATE
    assert constraint.advance(SAFE_STATE, 1, 'AAA') == SAFE_STATE