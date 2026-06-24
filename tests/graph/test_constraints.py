from codeine.constraints.base import PathConstraint
from codeine.constraints.mutations import MutationDistanceConstraint


def test_base_path_constraint_does_nothing():
    constraint = PathConstraint()

    state = constraint.initial_state
    assert state == ()

    assert constraint.advance(state, 1, 'ATA') == state
    assert constraint.is_satisfied(state)


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
    constraint = MutationDistanceConstraint('AAA', min_nts=1, min_codons=1)

    state = constraint.advance(constraint.initial_state, 1, 'TTT')

    assert state == (3, 1)


def test_mutation_distance_accumulates_across_codons():
    constraint = MutationDistanceConstraint('AAAATG', min_nts=1, min_codons=1)
    state = constraint.initial_state

    state = constraint.advance(state, 1, 'AAT',)
    state = constraint.advance(state, 2, 'ATA')

    assert state == (2, 2)


def test_mutation_distance_rejects_when_max_nts_exceeded():
    constraint = MutationDistanceConstraint('ATG', max_nts=0)
    state = constraint.advance(constraint.initial_state, 1, 'ATA',)

    assert state is None


def test_mutation_distance_rejects_when_max_codons_exceeded():
    constraint = MutationDistanceConstraint('ATG', max_codons=0)
    state = constraint.advance(constraint.initial_state, 1, 'ATA')

    assert state is None


def test_mutation_distance_accepts_final_enforces_min_nts():
    constraint = MutationDistanceConstraint('ATG', min_nts=1)

    assert not constraint.is_satisfied((0, None))
    assert constraint.is_satisfied((1, None))


def test_mutation_distance_accepts_final_enforces_min_codons():
    constraint = MutationDistanceConstraint('ATG', min_codons=1)

    assert not constraint.is_satisfied((None, 0))
    assert constraint.is_satisfied((None, 1))


def test_mutation_distance_accepts_final_enforces_both_minimums():
    constraint = MutationDistanceConstraint('ATG', min_nts=2, min_codons=1)

    assert not constraint.is_satisfied((1, 1))
    assert not constraint.is_satisfied((2, 0))
    assert constraint.is_satisfied((2, 1))
