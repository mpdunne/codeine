from graph.tracking import BannedSequenceTracker, AdvanceResult
from graph.codon import CodonGraph


def helper_find_path_for(tracker, banned_sequence):
    for path in tracker.paths:
        if path.sequence == banned_sequence:
            return path

    raise AssertionError(f'No path found for {banned_sequence!r}')


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


def test_starts_are_built_from_first_path_part():
    tracker = BannedSequenceTracker(CodonGraph(aa_seq='MIKEY'), ['TCAAA'])
    path = helper_find_path_for(tracker, 'TCAAA')
    first_node, first_choice = path.parts[0]

    assert (first_node, first_choice) in tracker.starts
    assert tracker.starts[(first_node, first_choice)]


def test_safe_choice_returns_empty_state():
    tracker = BannedSequenceTracker(CodonGraph(aa_seq='MIKEY'), ['TCAAA'])
    path = helper_find_path_for(tracker, 'TCAAA')
    first_node, _first_choice = path.parts[0]

    result = tracker.advance(node=first_node, state=tracker.initial_state, choice='ATG')

    assert result == AdvanceResult(banned=False, state=frozenset())


def test_choice_can_start_watch():
    tracker = BannedSequenceTracker(CodonGraph(aa_seq='MIKEY'), ['TCAAA'])
    path = helper_find_path_for(tracker, 'TCAAA')
    node, choice = path.parts[0]

    result = tracker.advance(node=node, state=tracker.initial_state, choice=choice)

    assert result.banned is False
    assert result.state


def test_choice_can_immediately_complete_banned_sequence():
    tracker = BannedSequenceTracker(CodonGraph(aa_seq='MIKEY'), ['ATG'])
    path = helper_find_path_for(tracker, 'ATG')
    node, choice = path.parts[0]

    result = tracker.advance(node=node, state=tracker.initial_state, choice=choice)

    assert result == AdvanceResult(banned=True)


def test_existing_watch_can_complete_banned_sequence():
    tracker = BannedSequenceTracker(CodonGraph(aa_seq='MIKEY'), ['TCAAA'])
    path = helper_find_path_for(tracker, 'TCAAA')

    node_1, choice_1 = path.parts[0]
    node_2, choice_2 = path.parts[1]

    result_1 = tracker.advance(node_1, tracker.initial_state, choice_1)
    result_2 = tracker.advance(node_2, result_1.state, choice_2)

    assert result_1.banned is False
    assert result_2 == AdvanceResult(banned=True)


def test_existing_watch_drops_if_choice_does_not_match():
    tracker = BannedSequenceTracker(CodonGraph(aa_seq='MIKEY'), ['TCAAA'])
    path = helper_find_path_for(tracker, 'TCAAA')

    node_1, choice_1 = path.parts[0]
    node_2, _choice_2 = path.parts[1]

    result_1 = tracker.advance(node_1, tracker.initial_state, choice_1)
    result_2 = tracker.advance(node_2, result_1.state, 'GAG')

    assert result_1.banned is False
    assert result_1.state
    assert result_2 == AdvanceResult(banned=False, state=frozenset())


def test_multiple_watches_can_be_active():
    graph = CodonGraph(aa_seq='MIKEY')
    tracker = BannedSequenceTracker(graph, ['TCAAA', 'TCAAG'])

    path = helper_find_path_for(tracker, 'TCAAA')
    node, choice = path.parts[0]

    result = tracker.advance(node, tracker.initial_state, choice)

    assert result.banned is False
    assert len(result.state) >= 2


def test_one_of_multiple_watches_can_complete_ban():
    graph = CodonGraph(aa_seq='MIKEY')
    tracker = BannedSequenceTracker(graph, ['TCAAA', 'TCAAG'])

    path = helper_find_path_for(tracker, 'TCAAA')

    node_1, choice_1 = path.parts[0]
    node_2, choice_2 = path.parts[1]

    result_1 = tracker.advance(node_1, tracker.initial_state, choice_1)
    result_2 = tracker.advance(node_2, result_1.state, choice_2)

    assert len(result_1.state) >= 2
    assert result_2 == AdvanceResult(banned=True)


def test_state_is_a_frozenset():
    tracker = BannedSequenceTracker(CodonGraph(aa_seq='MIKEY'), ['TCAAA'])
    path = helper_find_path_for(tracker, 'TCAAA')
    node, choice = path.parts[0]

    result = tracker.advance(node, tracker.initial_state, choice)

    assert isinstance(result.state, frozenset)


def test_tracker_finds_banned_sequence_crossing_left_context():
    graph = CodonGraph(aa_seq='MIKEY', context_l='TCA')
    tracker = BannedSequenceTracker(graph, ['TCAATG'])

    assert not tracker.is_trivial


def helper_walk_path(tracker, path):
    state = tracker.initial_state

    for node, choice in path.parts:
        result = tracker.advance(node, state, choice)

        if result.banned:
            return result

        state = result.state

    return result


def test_left_context_can_immediately_complete_ban():
    graph = CodonGraph(aa_seq='MIKEY', context_l='TCA')
    tracker = BannedSequenceTracker(graph, ['TCAATG'])

    path = helper_find_path_for(tracker, 'TCAATG')
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

    path = helper_find_path_for(tracker, 'ATACAAA')
    result = helper_walk_path(tracker, path)

    assert result == AdvanceResult(banned=True)
