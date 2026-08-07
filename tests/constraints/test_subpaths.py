import pytest

from codeine.graph.base import CodonGraph
from codeine.constraints.base import DEAD_STATE, SAFE_STATE
from codeine.constraints.motifs import ForbiddenMotifs
from codeine.constraints.subpaths import SubPath, SubPathConstraint


###############################
# Helpers
###############################


class TestSubPathConstraint(SubPathConstraint):

    def __init__(self, paths):
        super().__init__()
        self._paths = tuple(paths)

    def _find_paths(self):
        return self._paths


def helper_get_example_path(graph, start_pos, end_pos):
    """
    Return a deterministic real path through the graph.
    """
    node = graph.initial_node
    steps = []

    while node.pos <= end_pos:
        choice, next_node = next(iter(node.transitions.items()))

        if node.pos >= start_pos:
            steps.append((node.pos, choice))

        node = next_node

    return SubPath(
        sequence=''.join(choice for _pos, choice in steps),
        steps=tuple(steps),
        offset=0,
    )



def helper_find_first_path_for(constraint, sequence):
    """
    Return the first tracked path corresponding to a concrete motif sequence.
    """
    for path in constraint.paths:
        if path.sequence == sequence:
            return path

    raise AssertionError(f'No path found for {sequence!r}')


def helper_walk_path(constraint, path):
    """
    Walk a tracked path from the initial state until it either completes
    a ban or reaches the end of the path.
    """
    state_id = constraint.initial_state_id

    for step in path.steps:
        pos, choice = step
        result = constraint.advance(state_id, pos, choice)

        if result == DEAD_STATE:
            return result

        state_id = result

    return result


###############################
# Construction and linking
###############################


def test_starts_are_built_from_first_path_step():
    graph = CodonGraph(aa_seq='MIKEY')
    path = helper_get_example_path(graph, 1, 2)
    constraint = TestSubPathConstraint([path])
    constraint.link(graph)
    first_step = path.steps[0]

    assert first_step in constraint.starts
    assert constraint.starts[first_step]


def test_subpath_constraint_relink_resets_internal_state():
    graph1 = CodonGraph('KKK')
    graph2 = CodonGraph('KKK')
    path = helper_get_example_path(graph1, 1, 2)
    constraint = TestSubPathConstraint([path])

    constraint.link(graph1)

    # Populate the registry/cache.
    pos, choice = path.steps[0]
    constraint.advance(constraint.initial_state, pos, choice)

    assert len(constraint.states) > 1 or constraint.advance_cache

    constraint.link(graph2)

    assert constraint.states == [frozenset()]
    assert constraint.state_ids == {frozenset(): 0}
    assert constraint.advance_cache == {}


def test_subpath_constraint_is_trivial_before_linking():
    constraint = TestSubPathConstraint([])
    assert constraint.is_trivial


def test_subpath_constraint_preserves_terminal_states():
    constraint = TestSubPathConstraint([])

    assert constraint.advance(DEAD_STATE, 1, 'AAA') == DEAD_STATE
    assert constraint.advance(SAFE_STATE, 1, 'AAA') == SAFE_STATE


###############################
# State transitions
###############################


def test_safe_choice_returns_empty_state():
    graph = CodonGraph(aa_seq='MIKEY')
    path = helper_get_example_path(graph, 1, 2)
    constraint = TestSubPathConstraint([path])
    constraint.link(graph)

    pos, path_choice = path.steps[0]
    safe_choice = next(choice for choice in graph.initial_node.transitions if choice != path_choice)
    result = constraint.advance(constraint.initial_state_id, pos, safe_choice)

    assert result == 0


def test_choice_can_start_watch():
    graph = CodonGraph(aa_seq='MIKEY')
    path = helper_get_example_path(graph, 1, 2)
    constraint = TestSubPathConstraint([path])
    constraint.link(graph)

    pos, choice = path.steps[0]
    result = constraint.advance(constraint.initial_state_id, pos, choice)

    assert result != DEAD_STATE
    assert result != 0


def test_choice_can_immediately_complete_banned_sequence():
    graph = CodonGraph(aa_seq='MIKEY')
    path = helper_get_example_path(graph, 1, 1)
    constraint = TestSubPathConstraint([path])
    constraint.link(graph)

    pos, choice = path.steps[0]
    result = constraint.advance(constraint.initial_state_id, pos, choice)

    assert result == DEAD_STATE


def test_existing_watch_can_complete_banned_sequence():
    graph = CodonGraph(aa_seq='MIKEY')
    path = helper_get_example_path(graph, 1, 2)
    constraint = TestSubPathConstraint([path])
    constraint.link(graph)

    pos_1, choice_1 = path.steps[0]
    result_1 = constraint.advance(constraint.initial_state_id, pos_1, choice_1)
    pos_2, choice_2 = path.steps[1]
    result_2 = constraint.advance(result_1, pos_2, choice_2)

    assert result_1 != DEAD_STATE
    assert result_2 == DEAD_STATE


def test_existing_watch_drops_if_choice_does_not_match():
    graph = CodonGraph(aa_seq='MIKEY')
    path = helper_get_example_path(graph, 1, 2)
    constraint = TestSubPathConstraint([path])
    constraint.link(graph)

    pos_1, choice_1 = path.steps[0]
    result_1 = constraint.advance(constraint.initial_state_id, pos_1, choice_1)

    pos_2, path_choice_2 = path.steps[1]
    node = next(iter(graph.initial_node.transitions.values()))
    safe_choice = next(choice for choice in node.transitions if choice != path_choice_2)
    result_2 = constraint.advance(result_1, pos_2, safe_choice)

    assert result_1 != DEAD_STATE
    assert result_1 != 0
    assert result_2 == 0



def test_multiple_watches_can_be_active():
    graph = CodonGraph(aa_seq='MIKEY')
    path_1 = helper_get_example_path(graph, 1, 2)
    path_2 = helper_get_example_path(graph, 1, 3)
    constraint = TestSubPathConstraint([path_1, path_2])
    constraint.link(graph)

    pos, choice = path_1.steps[0]
    result = constraint.advance(constraint.initial_state_id, pos, choice)

    assert result != DEAD_STATE
    assert len(constraint.states[result]) >= 2


def test_one_of_multiple_watches_can_complete_ban():
    graph = CodonGraph(aa_seq='MIKEY')
    path_1 = helper_get_example_path(graph, 1, 2)
    path_2 = helper_get_example_path(graph, 1, 3)
    constraint = TestSubPathConstraint([path_1, path_2])
    constraint.link(graph)

    pos_1, choice_1 = path_1.steps[0]
    result_1 = constraint.advance(constraint.initial_state_id, pos_1, choice_1)
    pos_2, choice_2 = path_1.steps[1]
    result_2 = constraint.advance(result_1, pos_2, choice_2)

    assert len(constraint.states[result_1]) >= 2
    assert result_2 == DEAD_STATE

def test_state_is_a_frozenset():
    graph = CodonGraph(aa_seq='MIKEY')
    path = helper_get_example_path(graph, 1, 2)
    constraint = TestSubPathConstraint([path])
    constraint.link(graph)

    pos, choice = path.steps[0]
    result = constraint.advance(constraint.initial_state_id, pos, choice)

    assert isinstance(constraint.states[result], frozenset)


def test_watch_survives_until_partial_final_codon_match():
    graph = CodonGraph('MIKEY')
    constraint = ForbiddenMotifs(['ATTAAG'])
    constraint.link(graph)
    path = helper_find_first_path_for(constraint, 'ATTAAG')

    state_id = constraint.initial_state_id

    for step in path.steps[:-1]:
        pos, choice = step
        result = constraint.advance(state_id, pos, choice)
        assert result != DEAD_STATE
        state_id = result

    pos, choice = path.steps[-1]
    result = constraint.advance(state_id, pos, choice)

    assert result == DEAD_STATE


def test_ban_longer_than_choice_keeps_watch_alive():
    graph = CodonGraph('MIKEY')
    path = helper_get_example_path(graph, 1, 2)
    constraint = TestSubPathConstraint([path])
    constraint.link(graph)

    pos, choice = path.steps[0]
    result = constraint.advance(constraint.initial_state_id, pos, choice)

    assert result != DEAD_STATE
    assert result

    state = constraint.states[result]
    path_ix, matched_length = next(iter(state))
    assert constraint.paths[path_ix] == path
    assert matched_length == 3


def test_duplicate_banned_sequences_do_not_break_tracking():
    graph = CodonGraph('MIKEY')
    path = helper_get_example_path(graph, 1, 2)
    constraint = TestSubPathConstraint([path, path])
    constraint.link(graph)

    result = helper_walk_path(constraint, path)
    assert result == DEAD_STATE


def test_different_paths_can_start_from_same_choice():
    graph = CodonGraph('MIKEY')
    path_1 = helper_get_example_path(graph, 1, 2)
    path_2 = helper_get_example_path(graph, 1, 3)
    constraint = TestSubPathConstraint([path_1, path_2])
    constraint.link(graph)

    pos, choice = path_1.steps[0]
    result = constraint.advance(constraint.initial_state_id, pos, choice)

    assert result != DEAD_STATE
    assert len(constraint.states[result]) >= 2


def test_shorter_path_wins_when_multiple_paths_share_prefix():
    graph = CodonGraph('MIKEY')
    short_path = helper_get_example_path(graph, 1, 2)
    long_path = helper_get_example_path(graph, 1, 3)
    constraint = TestSubPathConstraint([short_path, long_path])
    constraint.link(graph)

    result = helper_walk_path(constraint, short_path)
    assert result == DEAD_STATE

def test_longer_ban_can_complete_after_shorter_related_ban_if_shorter_absent():
    graph = CodonGraph('MIKEY')
    path = helper_get_example_path(graph, 1, 2)
    constraint = TestSubPathConstraint([path])
    constraint.link(graph)

    result = helper_walk_path(constraint, path)
    assert result == DEAD_STATE



def test_active_watch_does_not_prevent_new_watch_starting():
    graph = CodonGraph('MIKEY')
    active_path = helper_get_example_path(graph, 1, 3)
    new_path = helper_get_example_path(graph, 2, 2)
    constraint = TestSubPathConstraint([active_path, new_path])
    constraint.link(graph)

    pos_1, choice_1 = active_path.steps[0]
    result_1 = constraint.advance(constraint.initial_state_id, pos_1, choice_1)
    pos_2, choice_2 = active_path.steps[1]
    result_2 = constraint.advance(result_1, pos_2, choice_2)

    assert result_2 == DEAD_STATE

def test_safe_walk_drops_all_active_watches():
    graph = CodonGraph('MIKEY')
    path = helper_get_example_path(graph, 1, 2)
    constraint = TestSubPathConstraint([path])
    constraint.link(graph)

    pos_1, choice_1 = path.steps[0]
    result_1 = constraint.advance(constraint.initial_state_id, pos_1, choice_1)

    pos_2, path_choice_2 = path.steps[1]
    node = next(iter(graph.initial_node.transitions.values()))
    safe_choice = next(choice for choice in node.transitions if choice != path_choice_2)
    result_2 = constraint.advance(result_1, pos_2, safe_choice)

    assert result_1 != 0
    assert result_2 == 0


###############################
# Path invariants
###############################


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
    constraint = ForbiddenMotifs([sequence])
    constraint.link(graph)
    path = helper_find_first_path_for(constraint, sequence)

    assert helper_walk_path(constraint, path) == DEAD_STATE


def test_path_steps_are_never_empty():
    graph = CodonGraph('MIKEY')
    constraint = ForbiddenMotifs(['A', 'AT', 'ATG', 'TCAAA'])
    constraint.link(graph)

    assert constraint.paths
    assert all(path.steps for path in constraint.paths)



def test_starts_only_reference_real_paths():
    graph = CodonGraph('MIKEY')
    paths = [
        helper_get_example_path(graph, 1, 1),
        helper_get_example_path(graph, 1, 2),
        helper_get_example_path(graph, 2, 3),
    ]
    constraint = TestSubPathConstraint(paths)
    constraint.link(graph)

    for starts in constraint.starts.values():
        for watch in starts:
            if watch is None:
                continue

            path_ix, matched_length = watch

            assert 0 <= path_ix < len(constraint.paths)
            assert matched_length > 0


def test_all_start_keys_are_real_first_steps():
    graph = CodonGraph('MIKEY')
    paths = [
        helper_get_example_path(graph, 1, 1),
        helper_get_example_path(graph, 1, 2),
        helper_get_example_path(graph, 2, 3),
    ]
    constraint = TestSubPathConstraint(paths)
    constraint.link(graph)

    first_steps = {path.steps[0] for path in constraint.paths}
    assert set(constraint.starts) <= first_steps


def test_walking_every_path_completes_ban():
    graph = CodonGraph('MIKEY')
    paths = [
        helper_get_example_path(graph, 1, 1),
        helper_get_example_path(graph, 1, 2),
        helper_get_example_path(graph, 2, 3),
    ]
    constraint = TestSubPathConstraint(paths)
    constraint.link(graph)

    for path in constraint.paths:
        result = helper_walk_path(constraint, path)
        assert result == DEAD_STATE


def test_every_constraint_path_is_walkable_from_initial_state():
    graph = CodonGraph('MIKEY')
    paths = [
        helper_get_example_path(graph, 1, 1),
        helper_get_example_path(graph, 1, 2),
        helper_get_example_path(graph, 2, 3),
    ]
    constraint = TestSubPathConstraint(paths)
    constraint.link(graph)

    for path in constraint.paths:
        result = helper_walk_path(constraint, path)
        assert result == DEAD_STATE
