import pytest

from codeine.graph.base import CodonGraph
from codeine.graph.nodes import ContextNode, CodonNode
from codeine.translation.tables import TranslationTable
from codeine.constraints.banned import BannedSequenceTracker, AdvanceResult, _find_matching_subpaths


def helper_find_first_path_for(tracker, sequence):
    """
    Return the first tracked path corresponding to a banned sequence.
    """
    for path in tracker.paths:
        if path.sequence == sequence:
            return path

    raise AssertionError(f'No path found for {sequence!r}')


def helper_walk_path(tracker, path):
    """
    Walk a tracked path from the initial state until it either completes
    a ban or reaches the end of the path.
    """
    state = tracker.initial_state

    for step in path.steps:
        result = tracker.advance(step, state)

        if result.banned:
            return result

        state = result.state

    return result


def test_tracker_is_trivial_without_banned_sequences():
    tracker = BannedSequenceTracker(CodonGraph(aa_seq='MIKEY'), [])

    assert tracker.is_trivial
    assert tracker.paths == ()
    assert tracker.starts == {}
    assert tracker.initial_state == frozenset()


def test_tracker_finds_paths_for_possible_banned_sequence():
    tracker = BannedSequenceTracker(CodonGraph(aa_seq='MIKEY'), ['TCAAA'])

    assert not tracker.is_trivial
    assert len(tracker.paths) > 0
    assert all(path.sequence == 'TCAAA' for path in tracker.paths)


def test_tracker_has_no_paths_for_impossible_banned_sequence():
    tracker = BannedSequenceTracker(CodonGraph(aa_seq='MIKEY'), ['CCCCCC'])

    assert tracker.is_trivial
    assert tracker.paths == ()
    assert tracker.starts == {}


def test_starts_are_built_from_first_path_step():
    tracker = BannedSequenceTracker(CodonGraph(aa_seq='MIKEY'), ['TCAAA'])
    path = helper_find_first_path_for(tracker, 'TCAAA')
    first_step = path.steps[0]

    assert first_step in tracker.starts
    assert tracker.starts[first_step]


def test_safe_choice_returns_empty_state():
    graph = CodonGraph(aa_seq='MIKEY')
    tracker = BannedSequenceTracker(graph, ['TCAAA'])
    path = helper_find_first_path_for(tracker, 'TCAAA')

    pos, _choice = path.steps[0]
    result = tracker.advance((pos, 'ATG'), tracker.initial_state)

    assert result == AdvanceResult(banned=False, state=frozenset())


def test_choice_can_start_watch():
    graph = CodonGraph(aa_seq='MIKEY')
    tracker = BannedSequenceTracker(graph, ['TCAAA'])
    path = helper_find_first_path_for(tracker, 'TCAAA')

    result = tracker.advance(path.steps[0], tracker.initial_state)

    assert result.banned is False
    assert result.state


def test_choice_can_immediately_complete_banned_sequence():
    graph = CodonGraph(aa_seq='MIKEY')
    tracker = BannedSequenceTracker(graph, ['ATG'])
    path = helper_find_first_path_for(tracker, 'ATG')

    result = tracker.advance(path.steps[0], tracker.initial_state)

    assert result == AdvanceResult(banned=True)


def test_existing_watch_can_complete_banned_sequence():
    graph = CodonGraph(aa_seq='MIKEY')
    tracker = BannedSequenceTracker(graph, ['TCAAA'])
    path = helper_find_first_path_for(tracker, 'TCAAA')

    result_1 = tracker.advance(path.steps[0], tracker.initial_state)
    result_2 = tracker.advance(path.steps[1], result_1.state)

    assert result_1.banned is False
    assert result_2 == AdvanceResult(banned=True)


def test_existing_watch_drops_if_choice_does_not_match():
    graph = CodonGraph(aa_seq='MIKEY')
    tracker = BannedSequenceTracker(graph, ['TCAAA'])
    path = helper_find_first_path_for(tracker, 'TCAAA')

    pos_2, _choice_2 = path.steps[1]

    result_1 = tracker.advance(path.steps[0], tracker.initial_state)
    result_2 = tracker.advance((pos_2, 'GAG'), result_1.state)

    assert result_1.banned is False
    assert result_1.state
    assert result_2 == AdvanceResult(banned=False, state=frozenset())


def test_multiple_watches_can_be_active():
    graph = CodonGraph(aa_seq='MIKEY')
    tracker = BannedSequenceTracker(graph, ['TCAAA', 'TCAAG'])

    path = helper_find_first_path_for(tracker, 'TCAAA')

    result = tracker.advance(path.steps[0], tracker.initial_state)

    assert result.banned is False
    assert len(result.state) >= 2


def test_one_of_multiple_watches_can_complete_ban():
    graph = CodonGraph(aa_seq='MIKEY')
    tracker = BannedSequenceTracker(graph, ['TCAAA', 'TCAAG'])

    path = helper_find_first_path_for(tracker, 'TCAAA')

    result_1 = tracker.advance(path.steps[0], tracker.initial_state)
    result_2 = tracker.advance(path.steps[1], result_1.state)

    assert len(result_1.state) >= 2
    assert result_2 == AdvanceResult(banned=True)


def test_state_is_a_frozenset():
    graph = CodonGraph(aa_seq='MIKEY')
    tracker = BannedSequenceTracker(graph, ['TCAAA'])
    path = helper_find_first_path_for(tracker, 'TCAAA')

    result = tracker.advance(path.steps[0], tracker.initial_state)

    assert isinstance(result.state, frozenset)


def test_tracker_finds_banned_sequence_crossing_left_context():
    graph = CodonGraph(aa_seq='MIKEY', context_l='TCA')
    tracker = BannedSequenceTracker(graph, ['TCAATG'])
    assert not tracker.is_trivial


def test_left_context_can_immediately_complete_ban():
    graph = CodonGraph(aa_seq='MIKEY', context_l='TCA')
    tracker = BannedSequenceTracker(graph, ['TCAATG'])
    path = helper_find_first_path_for(tracker, 'TCAATG')
    result = helper_walk_path(tracker, path)

    assert result == AdvanceResult(banned=True)


def test_tracker_finds_banned_sequence_crossing_right_context():
    graph = CodonGraph(aa_seq='MIKEY', context_r='AAA')
    tracker = BannedSequenceTracker(graph, ['ATACAAA'])
    assert not tracker.is_trivial

    graph = CodonGraph(aa_seq='MIKEY', context_r='AAA')
    tracker = BannedSequenceTracker(graph, ['GGGAAA'])
    assert tracker.is_trivial


def test_existing_watch_can_complete_in_right_context():
    graph = CodonGraph(aa_seq='MIKEY', context_r='AAA')
    tracker = BannedSequenceTracker(graph, ['ATACAAA'])
    path = helper_find_first_path_for(tracker, 'ATACAAA')
    result = helper_walk_path(tracker, path)

    assert result == AdvanceResult(banned=True)


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
    tracker = BannedSequenceTracker(graph, [sequence])
    path = helper_find_first_path_for(tracker, sequence)

    assert helper_walk_path(tracker, path) == AdvanceResult(banned=True)


@pytest.mark.parametrize(
    'sequence,expected_offset',
    [
        ('ATG', 0),
        ('TG', 1),
        ('G', 2),
    ],
)
def test_offsets_are_correct(sequence, expected_offset):
    tracker = BannedSequenceTracker(CodonGraph('MIKEY'), [sequence])
    path = helper_find_first_path_for(tracker, sequence)

    assert path.offset == expected_offset


def test_watch_survives_until_partial_final_codon_match():
    graph = CodonGraph('MIKEY')
    tracker = BannedSequenceTracker(graph, ['ATTAAG'])
    path = helper_find_first_path_for(tracker, 'ATTAAG')

    state = tracker.initial_state

    for step in path.steps[:-1]:
        result = tracker.advance(step, state)
        assert result.banned is False
        state = result.state

    result = tracker.advance(path.steps[-1], state)

    assert result == AdvanceResult(banned=True)


def test_ban_longer_than_choice_keeps_watch_alive():
    graph = CodonGraph('MIKEY')
    tracker = BannedSequenceTracker(graph, ['ATGATA'])
    path = helper_find_first_path_for(tracker, 'ATGATA')

    result = tracker.advance(path.steps[0], tracker.initial_state)

    assert result.banned is False
    assert result.state

    path_ix, matched_length = next(iter(result.state))
    assert tracker.paths[path_ix] == path
    assert matched_length == 3


def test_duplicate_banned_sequences_do_not_break_tracking():
    tracker = BannedSequenceTracker(CodonGraph('MIKEY'), ['TCAAA', 'TCAAA'])
    path = helper_find_first_path_for(tracker, 'TCAAA')

    result = helper_walk_path(tracker, path)
    assert result == AdvanceResult(banned=True)


def test_different_bans_can_start_from_same_choice():
    graph = CodonGraph('MIKEY')
    tracker = BannedSequenceTracker(graph, ['ATGA', 'ATGAT'])
    path = helper_find_first_path_for(tracker, 'ATGA')

    result = tracker.advance(path.steps[0], tracker.initial_state)

    assert result.banned is False
    assert len(result.state) >= 2


def test_shorter_ban_wins_when_multiple_bans_share_prefix():
    tracker = BannedSequenceTracker(CodonGraph('MIKEY'), ['ATGA', 'ATGATA'])
    path = helper_find_first_path_for(tracker, 'ATGA')

    result = helper_walk_path(tracker, path)
    assert result == AdvanceResult(banned=True)


def test_longer_ban_can_complete_after_shorter_related_ban_if_shorter_absent():
    tracker = BannedSequenceTracker(CodonGraph('MIKEY'), ['ATGATA'])
    path = helper_find_first_path_for(tracker, 'ATGATA')

    result = helper_walk_path(tracker, path)
    assert result == AdvanceResult(banned=True)


def test_unrelated_active_watch_does_not_prevent_new_watch_starting():
    graph = CodonGraph('MIKEY')
    tracker = BannedSequenceTracker(graph, ['TCAAA', 'GAA'])
    path = helper_find_first_path_for(tracker, 'TCAAA')

    result_1 = tracker.advance(path.steps[0], tracker.initial_state)
    result_2 = tracker.advance(path.steps[1], result_1.state)

    assert result_2 == AdvanceResult(banned=True)


def test_path_steps_are_never_empty():
    tracker = BannedSequenceTracker(CodonGraph('MIKEY'), ['A', 'AT', 'ATG', 'TCAAA'])

    assert tracker.paths
    assert all(path.steps for path in tracker.paths)


def test_starts_only_reference_real_paths():
    tracker = BannedSequenceTracker(CodonGraph('MIKEY'), ['A', 'ATG', 'TCAAA'])

    for results in tracker.starts.values():
        for result in results:
            if result.banned:
                continue

            for path_ix, matched_length in result.state:
                assert 0 <= path_ix < len(tracker.paths)
                assert 0 < matched_length <= len(tracker.paths[path_ix].sequence)


def test_all_start_keys_are_real_first_steps():
    tracker = BannedSequenceTracker(CodonGraph('MIKEY'), ['A', 'ATG', 'TCAAA'])

    first_steps = {path.steps[0] for path in tracker.paths}
    assert set(tracker.starts) <= first_steps


def test_every_found_path_really_contains_banned_sequence():
    graph = CodonGraph('MIKEY', context_l='AAGG', context_r='TTCC')
    tracker = BannedSequenceTracker(graph, ['GGATG', 'TACAAG', 'ATTAAG'])

    for path in tracker.paths:
        emitted = ''.join(choice for pos, choice in path.steps)
        visible = emitted[path.offset:]
        assert visible.startswith(path.sequence)


def test_walking_every_found_path_completes_ban():
    graph = CodonGraph('MIKEY', context_l='AAGG', context_r='TTCC')
    tracker = BannedSequenceTracker(graph, ['GGATG', 'TACAAG', 'ATTAAG'])

    for path in tracker.paths:
        result = helper_walk_path(tracker, path)
        assert result == AdvanceResult(banned=True)


def test_safe_walk_drops_all_active_watches():
    tracker = BannedSequenceTracker(CodonGraph('MIKEY'), ['TCAAA'])
    path = helper_find_first_path_for(tracker, 'TCAAA')

    pos_2, _choice_2 = path.steps[1]

    result_1 = tracker.advance(path.steps[0], tracker.initial_state)
    result_2 = tracker.advance((pos_2, 'GAG'), result_1.state)

    assert result_1.state
    assert result_2 == AdvanceResult(banned=False, state=frozenset())


def test_every_tracker_path_is_walkable_from_initial_state():
    graph = CodonGraph('MIKEY', context_l='AAGGTT', context_r='CCAAGG')
    banned = ['A', 'AT', 'ATG', 'TGATA', 'GGTTATG', 'TACCCA']
    tracker = BannedSequenceTracker(graph, banned)

    for path in tracker.paths:
        result = helper_walk_path(tracker, path)
        assert result == AdvanceResult(banned=True)


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
    for path, offset in matches:
        assert len(path) == 1
        node, codon = path[0]
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
    path, offset = matches[0]
    assert offset == 0

    matches = _find_matching_subpaths(graph, 'TGATAAAGGAATAC')
    assert len(matches) == 1
    path, offset = matches[0]
    assert offset == 1

    matches = _find_matching_subpaths(graph, 'GATAAAGGAATAC')
    assert len(matches) == 1
    path, offset = matches[0]
    assert offset == 2


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
    path, offset = matches[0]
    assert isinstance(path[0][0], ContextNode)
    assert isinstance(path[1][0], CodonNode)

    graph = CodonGraph('MIKEY', context_r='AATTAATTAATTAATT')
    matches = _find_matching_subpaths(graph, 'TAC' + 'AATTAATTAA')
    assert matches
    assert len(matches) == 1
    path, offset = matches[0]
    assert isinstance(path[0][0], CodonNode)
    assert isinstance(path[1][0], ContextNode)


def test_find_matching_subpaths_matches_entire_left_context_plus_cds():
    graph = CodonGraph('MIKEY', context_l='AAGG')

    matches = _find_matching_subpaths(graph, 'AAGGATG')
    assert len(matches) == 1

    path, offset = matches[0]
    assert offset == 0
    assert isinstance(path[0][0], ContextNode)
    assert isinstance(path[1][0], CodonNode)


def test_find_matching_subpaths_matches_cds_plus_entire_right_context():
    graph = CodonGraph('MIKEY', context_r='AAGG')

    matches = _find_matching_subpaths(graph, 'TACAAGG')
    assert len(matches) == 1

    path, offset = matches[0]
    assert offset == 0
    assert isinstance(path[-2][0], CodonNode)
    assert isinstance(path[-1][0], ContextNode)


def test_find_matching_subpaths_multiple_matches():
    graph = CodonGraph('MMMM')

    matches = _find_matching_subpaths(graph, 'ATGATG')
    assert len(matches) == 3

    for path, offset in matches:
        assert offset == 0
        assert len(path) == 2
        assert ''.join(codon for node, codon in path) == 'ATGATG'


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

LONG_AA_SEQUENCES = (
    'MIKEY' * 100,
    'MIKEY' * 250,
    'MIKEY' * 500,
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

    graph = CodonGraph(
        aa_seq,
        context_l=context_l,
        context_r=context_r,
        translation_table=tt,
    )

    view = graph.view()
    seqs = [view[i] for i in range(min(view.n_valid_sequences, 10))]

    for seq in seqs:
        matches = _find_matching_subpaths(graph, seq)

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
