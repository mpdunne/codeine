import random
import pickle

from copy import deepcopy

import pytest

from dataclasses import replace
from unittest.mock import MagicMock

from codeine import CodingSpace
from codeine.constraints import DirectRepeats, ForbiddenMotifs, Hairpins, MaxHomopolymer, TandemRepeats
from codeine.constraints._gc import _GCConstraint, _GC3Constraint
from codeine.constraints.mutations import MutationDistanceConstraint
from codeine.constraints.base import Constraint, DEAD_STATE, SAFE_STATE
from codeine.graph.base import CodonGraph
from codeine.graph.compiler import ViewCompiler
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


@pytest.mark.parametrize('rna', [False, True])
@pytest.mark.parametrize('contexts', [('', ''), ('AT', 'GC')])
def test_factorised_compiled_view_matches_flat_operations(rna, contexts):

    graph = CodonGraph('LRFK', context_l=contexts[0], context_r=contexts[1],
                       translation_table=TranslationTable(table_id=1, rna=rna))
    flat = graph.view(constraints=[DirectRepeats(3)])
    factorised = graph.view(constraints=[DirectRepeats(3)], compiler='factorised', seed=42)
    expected = list(flat.enumerate())

    assert factorised.n_valid_sequences == len(expected)
    assert list(factorised.enumerate()) == expected
    assert all(seq in factorised for seq in expected)
    assert 'AAA' not in factorised

    with pytest.raises(ValueError):
        factorised.contains('invalid')

    for index, seq in enumerate(expected):
        assert factorised[index] == seq
        assert factorised[index - len(expected)] == seq

    for index_slice in [slice(None), slice(None, None, -1), slice(2, 8), slice(8, 2), slice(1, None, 3)]:
        assert factorised[index_slice] == expected[index_slice]

    assert all(seq in expected for seq in factorised.sample(30))


def test_factorised_lifecycle_and_pickle():

    view = CodonGraph('KKK').view(compiler='factorised', seed=42)
    view.compile()
    original = view.copy()
    snapshot = view._compiled
    weights = CodonWeights({
        aa: {codon: float(codon != 'AAA') for codon in codons}
        for aa, codons in view.graph.tt.aa_to_codons.items()
    })
    view.set_weights(weights)
    view.compile()

    assert view._compiled.counter is snapshot.counter
    assert view.sample(5) == ['AAGAAGAAG'] * 5
    assert original.n_valid_sequences == 8

    view.pin_codons({1: 'AAA'})

    assert view.n_valid_sequences == 4
    assert original.n_valid_sequences == 8

    view.unpin_codons(1)
    view.add_constraints([DirectRepeats(3)])
    expected = CodonGraph('KKK').view(constraints=[DirectRepeats(3)])

    assert list(view.enumerate()) == list(expected.enumerate())
    assert snapshot.n_valid_sequences == 8

    restored = pickle.loads(pickle.dumps(original))

    assert original.sample(20) == restored.sample(20)


def test_factorised_rejects_unsupported_constraints():
    view = CodonGraph('MIKEY').view(compiler='factorised', constraints=[RejectChoicesConstraint({'ATA'})])

    with pytest.raises(ValueError, match='does not support'):
        view.compile()


@pytest.mark.parametrize('rna', [False, True])
@pytest.mark.parametrize('case', range(10))
def test_factorised_constraint_families_and_combinations(case, rna):

    reference = 'CTTCGTTTTAAA'.replace('T', 'U') if rna else 'CTTCGTTTTAAA'
    cases = [
        [ForbiddenMotifs(['AATT', 'CGT'])],
        [MaxHomopolymer(3)],
        [TandemRepeats(2, 2)],
        [Hairpins(2, 1, 5)],
        [_GCConstraint(min_count=3, max_count=6)],
        [_GC3Constraint(min_count=1, max_count=2)],
        [MutationDistanceConstraint(reference, min_nts=1, max_nts=4, min_codons=1, max_codons=2)],
        [_GCConstraint(min_count=99)],
        [DirectRepeats(3), ForbiddenMotifs(['AATT']), MaxHomopolymer(4), _GCConstraint(max_count=6)],
        [Hairpins(3, 1, 5), _GC3Constraint(min_count=1), MutationDistanceConstraint(reference, max_codons=2)],
    ]
    graph = CodonGraph('LRFK', context_l='AT', context_r='GC',
                       translation_table=TranslationTable(table_id=1, rna=rna))
    flat = graph.view(constraints=deepcopy(cases[case]))
    factorised = graph.view(constraints=deepcopy(cases[case]), compiler='factorised', seed=1)
    expected = list(flat.enumerate())

    assert factorised.n_valid_sequences == len(expected)
    assert list(factorised.enumerate()) == expected
    assert factorised[::-1] == expected[::-1]

    if expected:
        assert all(seq in expected for seq in factorised.sample(20))
    else:
        with pytest.raises(ValueError, match='empty'):
            factorised.sample()


def test_factorised_mutation_space_updates_and_copy():

    base = CodingSpace('LRFK', compiler='factorised', constraints=[MaxHomopolymer(4)], seed=1)
    reference = base.sample()
    mutants = base.mutants(reference, free_positions=[1, 2, 3], max_codons=2)
    flat = CodingSpace('LRFK', constraints=[MaxHomopolymer(4)]).mutants(
        reference, free_positions=[1, 2, 3], max_codons=2,
    )

    assert mutants.compiler == 'factorised'
    assert list(mutants.enumerate()) == list(flat.enumerate())
    assert mutants.copy().sample(10) == mutants.sample(10)


def test_factorised_long_count_constraint():

    view = CodonGraph('K' * 1100).view(compiler='factorised', constraints=[_GC3Constraint(max_count=1)], seed=1)

    assert view.n_valid_sequences == 1101
    assert view[0] == 'AAA' * 1100
    assert view[-1] == 'AAG' + 'AAA' * 1099
    assert view.sample() in view

    restored = pickle.loads(pickle.dumps(view))

    assert restored.sample() == view.sample()


@pytest.mark.parametrize('protein,repeat_length', [
    ('MIKEYMIKEY', 9),
    ('MIKEYAAAAAMIKEY', 9),
    ('MIKEYAAAAAMIKEY', 15),
    ('RNYKQT', 4),
    ('IHERQW', 4),
])
def test_factorised_original_prototype_cases(protein, repeat_length):
    graph = CodonGraph(protein)
    flat = graph.view(constraints=[DirectRepeats(repeat_length)])
    factorised = graph.view(compiler='factorised', constraints=[DirectRepeats(repeat_length)], seed=42)

    assert factorised.n_valid_sequences == flat.n_valid_sequences
    assert factorised[:20] == flat[:20]
    assert factorised[-1] == flat[-1]
    assert all(sequence in flat for sequence in factorised.sample(100))
