import pytest

from codeine.graph.constraints import MutationDistanceConstraint, PathConstraint
from codeine.graph.nodes import CodonNode, EndNode


def test_base_path_constraint_does_nothing():
    constraint = PathConstraint()

    state = constraint.initial_state
    assert state == ()

    node = CodonNode(aa='M', pos=1, codons=['ATG']),
    assert constraint.advance(state, node, 'ATA') == state
    assert constraint.accepts_final(state)


def test_mutation_distance_initial_state_tracks_only_requested_distances():
    assert MutationDistanceConstraint('ATG').initial_state == (None, None)
    assert MutationDistanceConstraint('ATG', min_nts=1).initial_state == (0, None)
    assert MutationDistanceConstraint('ATG', max_nts=1).initial_state == (0, None)
    assert MutationDistanceConstraint('ATG', min_codons=1).initial_state == (None, 0)
    assert MutationDistanceConstraint('ATG', max_codons=1).initial_state == (None, 0)
    assert MutationDistanceConstraint('ATG', min_nts=1, min_codons=1).initial_state == (0, 0)


def test_mutation_distance_ignores_non_codon_nodes():
    constraint = MutationDistanceConstraint('ATG', min_nts=1, min_codons=1)
    state = constraint.initial_state

    assert constraint.advance(state, EndNode(), '') == state


def test_mutation_distance_counts_nt_and_codon_differences():
    constraint = MutationDistanceConstraint('AAA', min_nts=1, min_codons=1)

    node = CodonNode(aa='K', pos=1, codons=['AAA'])
    state = constraint.advance(constraint.initial_state, node, 'TTT')

    assert state == (3, 1)


def test_mutation_distance_accumulates_across_codons():
    constraint = MutationDistanceConstraint('AAAATG', min_nts=1, min_codons=1)
    state = constraint.initial_state

    node = CodonNode(aa='K', pos=1, codons=['AAA'])
    state = constraint.advance(state, node, 'AAT',)

    node = CodonNode(aa='M', pos=2, codons=['ATG'])
    state = constraint.advance(state, node, 'ATA')

    assert state == (2, 2)


def test_mutation_distance_rejects_when_max_nts_exceeded():
    constraint = MutationDistanceConstraint('ATG', max_nts=0)

    node = CodonNode(aa='M', pos=1, codons=['ATG'])
    state = constraint.advance(constraint.initial_state, node, 'ATA',)

    assert state is None


def test_mutation_distance_rejects_when_max_codons_exceeded():
    constraint = MutationDistanceConstraint('ATG', max_codons=0)

    node = CodonNode(aa='M', pos=1, codons=['ATG'])
    state = constraint.advance(constraint.initial_state, node, 'ATA')

    assert state is None


def test_mutation_distance_accepts_final_enforces_min_nts():
    constraint = MutationDistanceConstraint('ATG', min_nts=1)

    assert not constraint.accepts_final((0, None))
    assert constraint.accepts_final((1, None))


def test_mutation_distance_accepts_final_enforces_min_codons():
    constraint = MutationDistanceConstraint('ATG', min_codons=1)

    assert not constraint.accepts_final((None, 0))
    assert constraint.accepts_final((None, 1))


def test_mutation_distance_accepts_final_enforces_both_minimums():
    constraint = MutationDistanceConstraint('ATG', min_nts=2, min_codons=1)

    assert not constraint.accepts_final((1, 1))
    assert not constraint.accepts_final((2, 0))
    assert constraint.accepts_final((2, 1))