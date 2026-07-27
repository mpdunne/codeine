import pytest

from codeine.graph.base import CodonGraph
from codeine.constraints.base import DEAD_STATE, SAFE_STATE
from codeine.constraints.motifs import ForbiddenMotifConstraint


###############################
# Helpers
###############################


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
    constraint = ForbiddenMotifConstraint(['TCAAA'])
    constraint.link(graph)
    path = helper_find_first_path_for(constraint, 'TCAAA')
    first_step = path.steps[0]

    assert first_step in constraint.starts
    assert constraint.starts[first_step]


def test_forbidden_motif_constraint_relink_resets_internal_state():
    constraint = ForbiddenMotifConstraint(['AAA'])

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


def test_forbidden_motif_constraint_is_trivial_before_linking():
    constraint = ForbiddenMotifConstraint(['AAA'])
    assert constraint.is_trivial


def test_forbidden_motif_constraint_preserves_terminal_states():
    constraint = ForbiddenMotifConstraint(['AAA'])

    assert constraint.advance(DEAD_STATE, 1, 'AAA') == DEAD_STATE
    assert constraint.advance(SAFE_STATE, 1, 'AAA') == SAFE_STATE


###############################
# State transitions
###############################


def test_safe_choice_returns_empty_state():
    graph = CodonGraph(aa_seq='MIKEY')
    constraint = ForbiddenMotifConstraint(['TCAAA'])
    constraint.link(graph)
    path = helper_find_first_path_for(constraint, 'TCAAA')

    pos, _choice = path.steps[0]
    result = constraint.advance(constraint.initial_state_id, pos, 'ATG')

    assert result == 0


def test_choice_can_start_watch():
    graph = CodonGraph(aa_seq='MIKEY')
    constraint = ForbiddenMotifConstraint(['TCAAA'])
    constraint.link(graph)
    path = helper_find_first_path_for(constraint, 'TCAAA')

    pos, choice = path.steps[0]
    result = constraint.advance(constraint.initial_state_id, pos, choice)

    assert result != DEAD_STATE
    assert result != 0


def test_choice_can_immediately_complete_banned_sequence():
    graph = CodonGraph(aa_seq='MIKEY')
    constraint = ForbiddenMotifConstraint(['ATG'])
    constraint.link(graph)
    path = helper_find_first_path_for(constraint, 'ATG')

    pos, choice = path.steps[0]
    result = constraint.advance(constraint.initial_state_id, pos, choice)

    assert result == DEAD_STATE


def test_existing_watch_can_complete_banned_sequence():
    graph = CodonGraph(aa_seq='MIKEY')
    constraint = ForbiddenMotifConstraint(['TCAAA'])
    constraint.link(graph)
    path = helper_find_first_path_for(constraint, 'TCAAA')

    pos_1, choice_1 = path.steps[0]
    result_1 = constraint.advance(constraint.initial_state_id, pos_1, choice_1)
    pos_2, choice_2 = path.steps[1]
    result_2 = constraint.advance(result_1, pos_2, choice_2)

    assert result_1 != DEAD_STATE
    assert result_2 == DEAD_STATE


def test_existing_watch_drops_if_choice_does_not_match():
    graph = CodonGraph(aa_seq='MIKEY')
    constraint = ForbiddenMotifConstraint(['TCAAA'])
    constraint.link(graph)
    path = helper_find_first_path_for(constraint, 'TCAAA')

    pos_2, _choice_2 = path.steps[1]

    pos_1, choice_1 = path.steps[0]
    result_1 = constraint.advance(constraint.initial_state_id, pos_1, choice_1)
    result_2 = constraint.advance(result_1, pos_2, 'GAG')

    assert result_1 != DEAD_STATE
    assert result_1 != 0
    assert result_2 == 0


def test_multiple_watches_can_be_active():
    graph = CodonGraph(aa_seq='MIKEY')
    constraint = ForbiddenMotifConstraint(['TCAAA', 'TCAAG'])
    constraint.link(graph)

    path = helper_find_first_path_for(constraint, 'TCAAA')

    pos, choice = path.steps[0]
    result = constraint.advance(constraint.initial_state_id, pos, choice)

    assert result != DEAD_STATE
    assert len(constraint.states[result]) >= 2


def test_one_of_multiple_watches_can_complete_ban():
    graph = CodonGraph(aa_seq='MIKEY')
    constraint = ForbiddenMotifConstraint(['TCAAA', 'TCAAG'])
    constraint.link(graph)

    path = helper_find_first_path_for(constraint, 'TCAAA')

    pos_1, choice_1 = path.steps[0]
    result_1 = constraint.advance(constraint.initial_state_id, pos_1, choice_1)
    pos_2, choice_2 = path.steps[1]
    result_2 = constraint.advance(result_1, pos_2, choice_2)

    assert len(constraint.states[result_1]) >= 2
    assert result_2 == DEAD_STATE


def test_state_is_a_frozenset():
    graph = CodonGraph(aa_seq='MIKEY')
    constraint = ForbiddenMotifConstraint(['TCAAA'])
    constraint.link(graph)
    path = helper_find_first_path_for(constraint, 'TCAAA')

    pos, choice = path.steps[0]
    result = constraint.advance(constraint.initial_state_id, pos, choice)

    assert isinstance(constraint.states[result], frozenset)


def test_watch_survives_until_partial_final_codon_match():
    graph = CodonGraph('MIKEY')
    constraint = ForbiddenMotifConstraint(['ATTAAG'])
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
    constraint = ForbiddenMotifConstraint(['ATGATA'])
    constraint.link(graph)
    path = helper_find_first_path_for(constraint, 'ATGATA')

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
    constraint = ForbiddenMotifConstraint(['TCAAA', 'TCAAA'])
    constraint.link(graph)
    path = helper_find_first_path_for(constraint, 'TCAAA')

    result = helper_walk_path(constraint, path)
    assert result == DEAD_STATE


def test_different_bans_can_start_from_same_choice():
    graph = CodonGraph('MIKEY')
    constraint = ForbiddenMotifConstraint(['ATGA', 'ATGAT'])
    constraint.link(graph)
    path = helper_find_first_path_for(constraint, 'ATGA')

    pos, choice = path.steps[0]
    result = constraint.advance(constraint.initial_state_id, pos, choice)

    assert result != DEAD_STATE
    assert len(constraint.states[result]) >= 2


def test_shorter_ban_wins_when_multiple_bans_share_prefix():
    graph = CodonGraph('MIKEY')
    constraint = ForbiddenMotifConstraint(['ATGA', 'ATGATA'])
    constraint.link(graph)
    path = helper_find_first_path_for(constraint, 'ATGA')

    result = helper_walk_path(constraint, path)
    assert result == DEAD_STATE


def test_longer_ban_can_complete_after_shorter_related_ban_if_shorter_absent():
    graph = CodonGraph('MIKEY')
    constraint = ForbiddenMotifConstraint(['ATGATA'])
    constraint.link(graph)
    path = helper_find_first_path_for(constraint, 'ATGATA')

    result = helper_walk_path(constraint, path)
    assert result == DEAD_STATE


def test_unrelated_active_watch_does_not_prevent_new_watch_starting():
    graph = CodonGraph('MIKEY')
    constraint = ForbiddenMotifConstraint(['TCAAA', 'GAA'])
    constraint.link(graph)
    path = helper_find_first_path_for(constraint, 'TCAAA')

    pos_1, choice_1 = path.steps[0]
    result_1 = constraint.advance(constraint.initial_state_id, pos_1, choice_1)
    pos_2, choice_2 = path.steps[1]
    result_2 = constraint.advance(result_1, pos_2, choice_2)

    assert result_2 == DEAD_STATE


def test_safe_walk_drops_all_active_watches():
    graph = CodonGraph('MIKEY')
    constraint = ForbiddenMotifConstraint(['TCAAA'])
    constraint.link(graph)
    path = helper_find_first_path_for(constraint, 'TCAAA')

    pos_2, _choice_2 = path.steps[1]

    pos_1, choice_1 = path.steps[0]
    result_1 = constraint.advance(constraint.initial_state_id, pos_1, choice_1)
    result_2 = constraint.advance(result_1, pos_2, 'GAG')

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
    constraint = ForbiddenMotifConstraint([sequence])
    constraint.link(graph)
    path = helper_find_first_path_for(constraint, sequence)

    assert helper_walk_path(constraint, path) == DEAD_STATE


def test_path_steps_are_never_empty():
    graph = CodonGraph('MIKEY')
    constraint = ForbiddenMotifConstraint(['A', 'AT', 'ATG', 'TCAAA'])
    constraint.link(graph)

    assert constraint.paths
    assert all(path.steps for path in constraint.paths)


def test_starts_only_reference_real_paths():
    graph = CodonGraph('MIKEY')
    constraint = ForbiddenMotifConstraint(['A', 'ATG', 'TCAAA'])
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
    constraint = ForbiddenMotifConstraint(['A', 'ATG', 'TCAAA'])
    constraint.link(graph)

    first_steps = {path.steps[0] for path in constraint.paths}
    assert set(constraint.starts) <= first_steps


def test_walking_every_found_path_completes_ban():
    graph = CodonGraph('MIKEY', context_l='AAGG', context_r='TTCC')
    constraint = ForbiddenMotifConstraint(['GGATG', 'TACAAG', 'ATTAAG'])
    constraint.link(graph)

    for path in constraint.paths:
        result = helper_walk_path(constraint, path)
        assert result == DEAD_STATE


def test_every_constraint_path_is_walkable_from_initial_state():
    graph = CodonGraph('MIKEY', context_l='AAGGTT', context_r='CCAAGG')
    banned = ['A', 'AT', 'ATG', 'TGATA', 'GGTTATG', 'TACCCA']
    constraint = ForbiddenMotifConstraint(banned)
    constraint.link(graph)

    for path in constraint.paths:
        result = helper_walk_path(constraint, path)
        assert result == DEAD_STATE
