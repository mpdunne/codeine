import random

import pytest

from dataclasses import replace
from unittest.mock import MagicMock

from codeine.constraints.base import Constraint, DEAD_STATE, SAFE_STATE
from codeine.graph.base import CodonGraph
from codeine.graph.compile import ViewCompiler
from codeine.graph.compiled import CompiledView
from codeine.translation.tables import TranslationTable
from codeine.translation.weights import CodonWeights


class RejectChoicesConstraint(Constraint):
    """
    Test constraint that rejects specified graph choices.
    """

    def __init__(self, rejected_choices):
        self.rejected_choices = frozenset(rejected_choices)

    @property
    def initial_state(self):
        return ()

    @property
    def is_trivial(self):
        return not self.rejected_choices

    def link(self, graph):
        pass

    def advance(self, state, pos, choice):
        if state == SAFE_STATE:
            return SAFE_STATE

        if choice in self.rejected_choices:
            return DEAD_STATE

        return state


def test_existing_state_does_not_interrupt_dense_state_ids():
    compiler = ViewCompiler(CodonGraph('MIKEY').view())

    assert compiler._get_or_register_state_id(compiler.initial_pos, ()) == (0, True)
    assert compiler._get_or_register_state_id(1, ()) == (1, True)
    assert compiler._get_or_register_state_id(compiler.initial_pos, ()) == (0, False)
    assert compiler._get_or_register_state_id(2, ()) == (2, True)

    assert compiler.states == [
        (compiler.initial_pos, ()),
        (1, ()),
        (2, ()),
    ]


def test_extend_matches_full_compile():
    graph = CodonGraph('MIKEY')
    constraint = RejectChoicesConstraint({'ATA'})

    base_compiled = ViewCompiler(graph.view()).compile()
    extended_compiled = ViewCompiler(graph.view()).extend(base_compiled, [constraint])
    full_compiled = ViewCompiler(graph.view(constraints=(constraint,))).compile()

    assert extended_compiled == full_compiled


def test_extend_with_multiple_constraints_matches_full_compile():
    graph = CodonGraph('MIKEY')
    base_constraint = RejectChoicesConstraint({'ATG'})
    new_constraints = (
        RejectChoicesConstraint({'AAG'}),
        RejectChoicesConstraint({'TAT'}),
    )

    base_compiled = ViewCompiler(
        graph.view(constraints=(base_constraint,))
    ).compile()

    extended_view = graph.view(
        constraints=(base_constraint, *new_constraints),
    )
    extended_compiled = ViewCompiler(extended_view).extend(
        base_compiled,
        new_constraints,
    )
    full_compiled = ViewCompiler(extended_view).compile()

    assert extended_compiled == full_compiled


def test_extend_from_constrained_compile_matches_full_compile():
    graph = CodonGraph('MIKEY')
    base_constraint = RejectChoicesConstraint({'ATA'})
    new_constraint = RejectChoicesConstraint({'AAG'})

    base_compiled = ViewCompiler(graph.view(constraints=[base_constraint])).compile()

    extended_compiled = ViewCompiler(graph.view()).extend(base_compiled, [new_constraint])
    full_compiled = ViewCompiler(graph.view(constraints=[base_constraint, new_constraint])).compile()

    assert extended_compiled == full_compiled


def test_extend_with_no_constraints_matches_shallow_compile():
    graph = CodonGraph('MIKEY')

    constraint = RejectChoicesConstraint({'ATA'})
    constrained_view = graph.view(constraints=[constraint])

    base_compiled = ViewCompiler(constrained_view).compile()

    extended_compiled = ViewCompiler(graph.view()).extend(base_compiled, [])
    shallow_compiled = ViewCompiler(graph.view()).compile_shallow(base_compiled)

    assert extended_compiled == shallow_compiled


def test_extend_ignores_trivial_constraints():
    view = CodonGraph('MIKEY').view()

    base_compiled = ViewCompiler(view).compile()
    extended_compiled = ViewCompiler(view).extend(base_compiled, [RejectChoicesConstraint(set())])
    shallow_compiled = ViewCompiler(view).compile_shallow(base_compiled)

    assert extended_compiled == shallow_compiled


def test_extend_does_not_add_sequences():
    graph = CodonGraph('MIKEY')

    base_compiled = ViewCompiler(graph.view()).compile()
    extended_compiled = ViewCompiler(graph.view()).extend(base_compiled, [RejectChoicesConstraint({'AAG'})])

    assert extended_compiled.n_valid_sequences <= base_compiled.n_valid_sequences


def test_extend_can_remove_all_sequences():
    graph = CodonGraph('MIKEY')
    unconstrained_view = graph.view()
    constraint = RejectChoicesConstraint({'ATG'})

    base_compiled = ViewCompiler(unconstrained_view).compile()
    extended_compiled = ViewCompiler(unconstrained_view).extend(base_compiled, [constraint])

    constrained_view = graph.view(constraints=[constraint])
    full_compiled = ViewCompiler(constrained_view).compile()

    assert extended_compiled == full_compiled
    assert extended_compiled.n_valid_sequences == 0


def test_compiled_state_ids_are_dense():
    view = CodonGraph('MIKEY').view(constraints=[RejectChoicesConstraint({'ATA'})])

    compiled = ViewCompiler(view).compile()

    referenced_state_ids = {
        child_id
        for child_results in compiled.child_results_by_state_id
        for _choice, child_id in child_results
    }

    assert compiled.initial_state_id == 0
    assert referenced_state_ids == set(range(1, len(compiled.states)))


def test_extended_state_ids_are_dense():
    view = CodonGraph('MIKEY').view()

    base_compiled = ViewCompiler(view).compile()
    extended_compiled = ViewCompiler(view).extend(base_compiled, [RejectChoicesConstraint({'ATA'})])

    referenced_state_ids = {
        child_id
        for child_results in extended_compiled.child_results_by_state_id
        for _choice, child_id in child_results
    }

    assert extended_compiled.initial_state_id == 0
    assert referenced_state_ids == set(range(1, len(extended_compiled.states)))


@pytest.mark.parametrize('rna', [False, True])
def test_compiled_result_supports_queries_without_the_view(rna):
    graph = CodonGraph(
        'MK',
        translation_table=TranslationTable(rna=rna),
        context_l='CCC',
        context_r='GGG',
    )
    compiled: CompiledView = ViewCompiler(graph.view()).compile()
    expected = ['AUGAAA', 'AUGAAG'] if rna else ['ATGAAA', 'ATGAAG']

    assert compiled.n_valid_sequences == 2
    assert list(compiled.enumerate()) == expected
    assert compiled.sequence_at(-1) == expected[-1]
    assert compiled.sequences_at(slice(None, None, -1)) == expected[::-1]
    assert list(compiled.enumerate_range(1, 2)) == expected[1:]
    assert compiled.contains('atg aaa')
    assert not compiled.contains('ATGCCC')
    assert not compiled.contains('ATG')


def test_compiled_result_is_a_snapshot_of_view_configuration():
    view = CodonGraph('MK').view()
    view.compile()
    compiled = view._compiled

    view.pin_codons({2: 'AAA'})
    view.compile()

    assert view._compiled is not compiled
    assert view.n_valid_sequences == 1
    assert compiled.n_valid_sequences == 2
    assert list(compiled.enumerate()) == ['ATGAAA', 'ATGAAG']
    assert compiled.contains('ATGAAG')


@pytest.mark.parametrize('aaa_weight, expected', [
    (1, ['ATGAAA', 'ATGAAA', 'ATGAAG']),
    (0.25, ['ATGAAA', 'ATGAAG', 'ATGAAG']),
])
def test_compiled_sampling_uses_the_supplied_random_generator(aaa_weight, expected):
    table = TranslationTable()
    data = {aa: {codon: 1 for codon in codons} for aa, codons in table.aa_to_codons.items()}
    data['K'] = {'AAA': aaa_weight, 'AAG': 1}

    compiled = ViewCompiler(CodonGraph('MK').view(weights=CodonWeights(data))).compile()
    rng = MagicMock(spec=random.Random)
    rng.random.side_effect = [0.1, 0.3, 0.9]

    assert compiled.sample(rng, n=3) == expected
    assert rng.random.call_count == 3


def test_compiled_singleton_sampling_does_not_advance_random_generator():
    compiled = ViewCompiler(CodonGraph('M').view()).compile()
    rng = random.Random(123)
    state = rng.getstate()

    assert compiled.sample(rng, n=3) == ['ATG'] * 3
    assert rng.getstate() == state


def test_sequence_at_raises_on_unexpected_dead_end():
    compiled = ViewCompiler(CodonGraph('M').view()).compile()
    choice_results = list(compiled.choice_results_by_state_id)
    choice_results[compiled.initial_state_id] = ()

    compiled = replace(compiled, choice_results_by_state_id=tuple(choice_results))

    with pytest.raises(RuntimeError):
        compiled.sequence_at(0)
