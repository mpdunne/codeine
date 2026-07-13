from codeine.constraints.mutations import MutationDistanceConstraint
from codeine.constraints.base import DEAD_STATE, SAFE_STATE


def test_mutation_distance_ignores_non_codon_positions():
    constraint = MutationDistanceConstraint('ATG', min_nts=1, min_codons=1)
    state = constraint.initial_state

    assert constraint.advance(state, 0, '') == state
    assert constraint.advance(state, 2, '') == state


def test_mutation_distance_initial_state_tracks_only_requested_distances():
    assert MutationDistanceConstraint('ATG').initial_state == (None, None)
    assert MutationDistanceConstraint('ATG', min_nts=1).initial_state == (0, None)
    assert MutationDistanceConstraint('ATG', max_nts=1).initial_state == (0, None)
    assert MutationDistanceConstraint('ATG', min_codons=1).initial_state == (None, 0)
    assert MutationDistanceConstraint('ATG', max_codons=1).initial_state == (None, 0)
    assert MutationDistanceConstraint('ATG', min_nts=1, min_codons=1).initial_state == (0, 0)


def test_mutation_distance_counts_nt_and_codon_differences():
    constraint = MutationDistanceConstraint('AAATTAAAT', min_nts=1, min_codons=1, max_nts=3, max_codons=2)
    state = constraint.advance(constraint.initial_state, 1, 'TTT')

    assert state == (3, 1)


def test_mutation_distance_accumulates_across_codons():
    constraint = MutationDistanceConstraint('AAAATGTTAAAT', min_nts=1, min_codons=1, max_nts=3, max_codons=2)
    state = constraint.initial_state

    state = constraint.advance(state, 1, 'AAT',)
    state = constraint.advance(state, 2, 'ATA')

    assert state == (2, 2)


def test_mutation_distance_dead_when_max_nts_exceeded():
    constraint = MutationDistanceConstraint('ATG', max_nts=0)
    state = constraint.advance(constraint.initial_state, 1, 'ATA',)

    assert state is DEAD_STATE


def test_mutation_distance_dead_when_max_codons_exceeded():
    constraint = MutationDistanceConstraint('ATG', max_codons=0)
    state = constraint.advance(constraint.initial_state, 1, 'ATA')

    assert state is DEAD_STATE


def test_mutation_distance_dead_when_min_nts_unreachable():
    constraint = MutationDistanceConstraint('ATG', min_nts=1)
    state = constraint.advance(constraint.initial_state, 1, 'ATG')

    assert state == DEAD_STATE


def test_mutation_distance_dead_when_min_codons_unreachable():
    constraint = MutationDistanceConstraint('ATG', min_codons=1)
    state = constraint.advance(constraint.initial_state, 1, 'ATG')

    assert state == DEAD_STATE


def test_mutation_distance_always_safe_when_no_constraints_provided():
    constraint = MutationDistanceConstraint('ATGAAACCCAAA')

    state = constraint.advance(constraint.initial_state, 1, 'ATA')
    assert state == SAFE_STATE

    state = constraint.advance(constraint.initial_state, 1, 'ATG')
    assert state == SAFE_STATE


def test_mutation_distance_safe_when_minimums_are_reached():
    constraint = MutationDistanceConstraint('ATGAAA', min_nts=1, min_codons=1)
    state = constraint.advance(constraint.initial_state, 1, 'ATA')

    assert state == SAFE_STATE


def test_mutation_distance_not_safe_when_maximum_can_still_be_exceeded():
    constraint = MutationDistanceConstraint('ATGAAA', min_nts=1, max_nts=2)
    state = constraint.advance(constraint.initial_state, 1, 'ATA')

    assert state == (1, None)


def test_mutation_distance_safe_when_maximum_is_unreachable():
    constraint = MutationDistanceConstraint('ATG', min_nts=1, max_nts=2)
    state = constraint.advance(constraint.initial_state, 1, 'ATA')

    assert state == SAFE_STATE


def test_mutation_distance_safe_at_valid_final_state():
    constraint = MutationDistanceConstraint('ATG', min_nts=1, min_codons=1)
    state = constraint.advance(constraint.initial_state, 1, 'ATA')

    assert state == SAFE_STATE

def test_mutation_distance_remains_live_until_minimum_becomes_unreachable():
    constraint = MutationDistanceConstraint('ATGAAATTAGGC', min_nts=4)
    state = constraint.initial_state

    state = constraint.advance(state, 1, 'ATA')  # 1 nt diff, 3 codons remain
    assert state == (1, None)

    state = constraint.advance(state, 2, 'AAA')  # 1 nt diff, 2 codons remain
    assert state == (1, None)

    state = constraint.advance(state, 3, 'TTA')  # 1 nt diff, 1 codon remains
    assert state == (1, None)

    state = constraint.advance(state, 4, 'GGC')  # Still only 1 nt diff
    assert state is DEAD_STATE


def test_mutation_distance_becomes_safe_partway_through_sequence():
    constraint = MutationDistanceConstraint('ATGAAATTAGGC', min_nts=3, min_codons=2)
    state = constraint.initial_state

    state = constraint.advance(state, 1, 'ATA')  # 1 nt, 1 codon
    assert state == (1, 1)

    state = constraint.advance(state, 2, 'TTT')  # +3 nt, +1 codon
    assert state is SAFE_STATE


def test_mutation_distance_tracks_nt_and_codon_limits_independently():
    constraint = MutationDistanceConstraint('AAAAAATTT', min_nts=3, min_codons=2, max_nts=4, max_codons=2)
    state = constraint.initial_state

    state = constraint.advance(state, 1, 'AAT')  # 1 nt, 1 codon
    assert state == (1, 1)

    state = constraint.advance(state, 2, 'ATA')  # +1 nt, +1 codon
    assert state == (2, 2)

    state = constraint.advance(state, 3, 'TTA')  # +1 nt, +1 codon
    assert state is DEAD_STATE


def test_mutation_distance_becomes_safe_when_maxima_cannot_be_exceeded():
    constraint = MutationDistanceConstraint('ATGAAATTAGGC', max_nts=4, max_codons=2)
    state = constraint.initial_state

    state = constraint.advance(state, 1, 'ATA')  # 1 nt, 1 codon
    assert state == (1, 1)

    state = constraint.advance(state, 2, 'AAA')  # unchanged
    assert state == (1, 1)

    state = constraint.advance(state, 3, 'TTA')  # unchanged
    assert state is SAFE_STATE
