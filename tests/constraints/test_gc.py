import pytest

from codeine.graph.base import CodonGraph
from codeine.constraints._gc import _GC3Constraint, _GCConstraint

from tests.data import NORMAL_PROTEINS, LARGE_PROTEINS


def gc_count(sequence):
    return sum(nt in 'GC' for nt in sequence)


def gc3_count(sequence):
    return sum(sequence[pos] in 'GC' for pos in range(2, len(sequence), 3))


@pytest.mark.parametrize(
    'constraint',
    [
        _GCConstraint(),
        _GCConstraint(min_count=3),
        _GCConstraint(max_count=4),
        _GCConstraint(min_perc=30, max_perc=70),
        _GCConstraint(min_frac=0.3, max_count=5),
    ],
)
@pytest.mark.parametrize('aa_seq', (
    'M',
    'MIKEY',
    'SASSY',
))
def test_gc_constraint_matches_enumerate_and_reject(aa_seq, constraint):
    graph = CodonGraph(aa_seq)

    min_count = constraint._resolve_min_count(len(aa_seq) * 3)
    max_count = constraint._resolve_max_count(len(aa_seq) * 3)

    unconstrained_view = graph.view()
    expected = {sequence for sequence in unconstrained_view.enumerate()
                if min_count <= gc_count(sequence) <= max_count}

    graph = CodonGraph(aa_seq)

    constrained_view = graph.view()
    constrained_view.set_constraints([constraint])
    observed = set(constrained_view.enumerate())

    assert observed == expected


@pytest.mark.parametrize(
    'constraint',
    [
        _GC3Constraint(),
        _GC3Constraint(min_count=2),
        _GC3Constraint(max_count=2),
        _GC3Constraint(min_perc=25, max_perc=75),
        _GC3Constraint(min_frac=0.25, max_count=3),
    ],
)
@pytest.mark.parametrize('aa_seq', (
    'M',
    'MIKEY',
    'SASSY',
))
def test_gc3_constraint_matches_enumerate_and_reject(aa_seq, constraint):

    graph = CodonGraph(aa_seq)

    min_count = constraint._resolve_min_count(len(aa_seq))
    max_count = constraint._resolve_max_count(len(aa_seq))

    unconstrained_view = graph.view()
    expected = {sequence for sequence in unconstrained_view.enumerate()
                if min_count <= gc3_count(sequence) <= max_count}

    graph = CodonGraph(aa_seq)

    constrained_view = graph.view()
    constrained_view.set_constraints([constraint])
    observed = set(constrained_view.enumerate())

    assert observed == expected


@pytest.mark.parametrize(
    'constraint,count_fn,total_fn',
    [
        (_GCConstraint(min_perc=40, max_perc=60), gc_count, lambda sequence: len(sequence)),
        (_GC3Constraint(min_perc=40, max_perc=60), gc3_count, lambda sequence: len(sequence) // 3)
    ],
)
@pytest.mark.parametrize('aa_seq', (
    'MIKEY' * 100,
    *NORMAL_PROTEINS.values(),
    *LARGE_PROTEINS.values(),
))
def test_sampled_sequences_satisfy_constraints(aa_seq, constraint, count_fn, total_fn):
    graph = CodonGraph(aa_seq)
    view = graph.view()

    view.set_constraints([constraint])

    for sequence in view.sample(n=100):
        count = count_fn(sequence)
        total = total_fn(sequence)

        assert count >= constraint._resolve_min_count(total)
        assert count <= constraint._resolve_max_count(total)


def test_gc_constraint_counts_single_final_codon():
    view = CodonGraph('M').view()
    view.set_constraints([_GCConstraint(min_count=2)])

    assert set(view.enumerate()) == set()


def test_combined_gc_constraints_match_enumerate_and_reject():
    aa_seq = 'SASSAFRAS'

    gc_constraint = _GCConstraint(min_perc=40, max_perc=60)
    gc3_constraint = _GC3Constraint(min_perc=20, max_perc=80)

    graph = CodonGraph(aa_seq)

    gc_min = gc_constraint._resolve_min_count(len(aa_seq) * 3)
    gc_max = gc_constraint._resolve_max_count(len(aa_seq) * 3)
    gc3_min = gc3_constraint._resolve_min_count(len(aa_seq))
    gc3_max = gc3_constraint._resolve_max_count(len(aa_seq))

    expected = {
        sequence
        for sequence in graph.view().enumerate()
        if gc_min <= gc_count(sequence) <= gc_max
        and gc3_min <= gc3_count(sequence) <= gc3_max
    }

    view = CodonGraph(aa_seq).view()
    view.set_constraints([gc_constraint, gc3_constraint])
    view.compile()

    assert set(view.enumerate()) == expected