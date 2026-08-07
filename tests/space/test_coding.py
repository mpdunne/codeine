import pickle

import pytest

from Bio.Seq import Seq

from codeine.translation.weights import CodonWeights
from codeine.translation.tables import TranslationTable
from codeine.space.coding import CodingSpace
from codeine.motifs.restriction import RestrictionSite
from codeine.constraints.base import Constraint, SAFE_STATE
from codeine.constraints.motifs import ForbiddenMotifs
from codeine.constraints.homopolymers import MaxHomopolymer


@pytest.mark.parametrize('aa_seq', ('MIKEY', 'MILDRED', 'STEVEN', 'WILLIAM'))
def test_coding_space_sequences_translate_correctly(aa_seq):
    space = CodingSpace(aa_seq=aa_seq)
    for _ in range(1000):
        cds = space.sample()
        translated = Seq(cds).translate()
        assert translated == aa_seq


def test_coding_space_fixed_codons_are_fixed():
    space = CodingSpace('MIKEY', fixed_codons={2: 'ATA'})

    for _ in range(1000):
        cds = space.sample()
        codons = [cds[i:i + 3] for i in range(0, len(cds), 3)]

        assert codons[1] == 'ATA'
        assert str(Seq(cds).translate()) == 'MIKEY'


def test_coding_space_generates_different_sequences():
    space = CodingSpace('MIKEY')
    generated = [space.sample() for _ in range(1000)]
    assert len(generated) > 1


def test_coding_space_generates_all_cdss_for_small_aa_seq():
    space = CodingSpace('MF')
    generated = {space.sample() for _ in range(1000)}
    expected = {'ATGTTT', 'ATGTTC'}
    assert generated == expected


def test_coding_space_pinned_codons_are_fixed():
    space = CodingSpace('MIKEY')
    space.pin_codons({3: 'AAA'})
    for _ in range(100):
        cds = space.sample()
        codons = [cds[i:i + 3] for i in range(0, len(cds), 3)]
        assert codons[2] == 'AAA'


def test_coding_space_unpin_codons_restores_sampling():
    space = CodingSpace('MIKEY')
    space.pin_codons({3: 'AAA'})
    space.unpin_codons(3)
    sampled = set()
    for _ in range(100):
        cds = space.sample()
        codons = [cds[i:i + 3] for i in range(0, len(cds), 3)]
        sampled.add(codons[2])
    assert sampled == {'AAA', 'AAG'}


def test_coding_space_clear_pins_restores_sampling():
    space = CodingSpace('MIKEY')

    space.pin_codons({3: 'AAA'})
    space.clear_pins()

    sampled = set()
    for _ in range(100):
        cds = space.sample()
        codons = [cds[i:i + 3] for i in range(0, len(cds), 3)]
        sampled.add(codons[2])

    assert sampled == {'AAA', 'AAG'}


def test_coding_space_rejects_invalid_pin():
    space = CodingSpace('MIKEY')
    with pytest.raises(ValueError):
        space.pin_codons({3: 'GCT'})


def test_coding_space_rejects_out_of_range_pin():
    space = CodingSpace('MIKEY')

    with pytest.raises(ValueError):
        space.pin_codons({0: 'ATG'})

    with pytest.raises(ValueError):
        space.pin_codons({6: 'ATG'})


def test_coding_space_sample_excludes_context_by_default():
    space = CodingSpace(
        aa_seq='MF',
        context_l='AAAA',
        context_r='CCCC',
    )

    cds = space.sample()

    assert cds in {'ATGTTT', 'ATGTTC'}
    assert not cds.startswith('AAAA')
    assert not cds.endswith('CCCC')


def test_coding_space_mutants_raises_if_seq_is_invalid():
    space = CodingSpace('MIKEY')

    with pytest.raises(ValueError):
        _ = space.mutants('', [1, 2])

    with pytest.raises(ValueError):
        _ = space.mutants('ATG', [1, 2])

    with pytest.raises(ValueError):
        _ = space.mutants('ATTATTAAAGAATAT', [1, 2])

    with pytest.raises(ValueError):
        _ = space.mutants('ATGATTAAAGAATATATG', [1, 2])


def test_mutants_normalises_dna_reference_for_rna_space():
    space = CodingSpace('M', rna=True)

    mutants = space.mutants('ATG')

    assert mutants.cds == 'AUG'
    assert mutants.contains('AUG')


def test_coding_space_mutants_returns_mutation_space():
    space = CodingSpace('MIKEY')
    muts = space.mutants('ATGATTAAAGAATAT', free_positions=[2])

    assert muts.cds == 'ATGATTAAAGAATAT'
    assert muts.free_positions == {2}


def test_coding_space_mutants_does_not_change_base_space():
    aa_seq = 'MIKEY'
    ref_cds = 'ATGATTAAAGAATAT'

    space = CodingSpace(aa_seq)

    sampled_seqs = [space.sample() for _ in range(1000)]
    assert not all(s[6:] == ref_cds[6:] for s in sampled_seqs)
    assert not all(s[3:6] == ref_cds[3:6] for s in sampled_seqs)

    mut = space.mutants(ref_cds, [2])
    sampled_seqs = [mut.sample() for _ in range(1000)]
    assert all(s[6:] == ref_cds[6:] for s in sampled_seqs)
    assert not all(s[3:6] == ref_cds[3:6] for s in sampled_seqs)

    sampled_seqs = [space.sample() for _ in range(1000)]
    assert not all(s[6:] == ref_cds[6:] for s in sampled_seqs)
    assert not all(s[3:6] == ref_cds[3:6] for s in sampled_seqs)


def test_coding_space_mutants_passes_constraints():
    space = CodingSpace('MIKEY')
    muts = space.mutants(
        'ATGATTAAAGAATAT',
        free_positions=[2, 3],
        min_nts=1,
        max_nts=3,
        min_codons=1,
        max_codons=2,
    )

    assert set(muts.free_positions) == {2, 3}
    assert muts.min_nts == 1
    assert muts.max_nts == 3
    assert muts.min_codons == 1
    assert muts.max_codons == 2


def test_coding_space_getitem():
    space = CodingSpace('MM')
    assert space[0] == 'ATGATG'


def test_coding_space_iter():
    space = CodingSpace('MIKEY')
    seqs = [*space]
    assert len(seqs) == len(set(seqs)) == 24


def test_coding_space_enumerate():
    space = CodingSpace('F')
    seqs = [*space.enumerate()]
    assert len(seqs) == 2
    assert set(seqs) == {'TTT', 'TTC'}


def test_coding_space_contains():
    space = CodingSpace('F')
    assert space.contains('TTT')
    assert space.contains('TTC')
    assert not space.contains('ATG')


def test_coding_space_sample_many():
    space = CodingSpace('MIKEY', seed=88)
    seqs = space.sample(n=10)
    assert len(seqs) == 10
    assert all(seq in space for seq in seqs)


def test_coding_space_count():
    space = CodingSpace('MIKEY')
    assert space.count() == space.n_valid_sequences == 24


def helper_get_banned_sequences_from_constraints(space):
    banned_sequence_constraints = [c for c in space.view.constraints if isinstance(c, ForbiddenMotifs)]
    if len(banned_sequence_constraints) == 0:
        return set()
    else:
        seqs = []
        for c in banned_sequence_constraints:
            seqs += list(c.forbidden_sequences)
        return set(seqs)


def test_forbidden_motif_constraint_is_stored():
    space = CodingSpace('MIKEY', constraints=[ForbiddenMotifs([RestrictionSite.EcoRI, 'AAAA'])])
    assert helper_get_banned_sequences_from_constraints(space) == {'AAAA', 'GAATTC'}


def test_homopolymer_constraint_is_stored():
    constraint = MaxHomopolymer(4)
    space = CodingSpace('MIKEY', constraints=[constraint])
    assert space.constraints == (constraint,)


def test_homopolymer_constraint_is_expanded_correctly():
    space = CodingSpace('MIKEY')
    banned_sequences = helper_get_banned_sequences_from_constraints(space)
    assert not banned_sequences

    space = CodingSpace('MIKEY', constraints=[MaxHomopolymer(4)])
    banned_sequences = helper_get_banned_sequences_from_constraints(space)
    assert all(nt * 5 in banned_sequences for nt in 'ACGT')


def test_mixed_constraints():
    space = CodingSpace('MIKEY', constraints=[
        ForbiddenMotifs([RestrictionSite.BsaI, 'GGTTCC']),
        MaxHomopolymer(4),
    ])
    banned_sequences = helper_get_banned_sequences_from_constraints(space)
    assert banned_sequences == {'GAGACC', 'GGTCTC', 'GGTTCC', 'AAAAA', 'CCCCC', 'GGGGG', 'TTTTT'}


def test_forbidden_motif_constraint_repr():
    space = CodingSpace('MIKEY', constraints=[ForbiddenMotifs([RestrictionSite.EcoRI, 'AAAA'])])

    text = repr(space)
    assert 'Forbidden motifs:' in text
    assert 'EcoRI' in text
    assert 'GAATTC' in text
    assert 'AAAA' in text


def test_homopolymer_constraint_repr():
    space = CodingSpace('MIKEY', constraints=[MaxHomopolymer(4)])

    text = repr(space)
    assert 'Maximum homopolymer' in text
    assert '4' in text


def test_coding_space_seed_consistent():
    samples_by_rep = []

    for _ in range(10):
        space = CodingSpace('MIKEY', seed=8675309)
        samples = [space.sample() for _ in range(100)]
        samples_by_rep.append(samples)

    assert all(samples == samples_by_rep[0] for samples in samples_by_rep)


def test_coding_space_no_seed_not_consistent():
    samples_by_rep = []

    for _ in range(10):
        space = CodingSpace('MIKEY')
        samples = [space.sample() for _ in range(100)]
        samples_by_rep.append(samples)

    assert not all(samples == samples_by_rep[0] for samples in samples_by_rep)


def test_coding_space_pickle_preserves_random_state():
    space = CodingSpace('MIKEY', seed=8675309)

    _ = [space.sample() for _ in range(100)]

    loaded = pickle.loads(pickle.dumps(space))

    assert [loaded.sample() for _ in range(100)] == [space.sample() for _ in range(100)]


def test_coding_space_pickle_preserves_pins():
    space = CodingSpace('MIKEY', seed=8675309)
    space.pin_codons({2: 'ATC'})

    loaded = pickle.loads(pickle.dumps(space))

    assert loaded.view.pinned_codons == space.view.pinned_codons
    assert loaded.n_valid_sequences == space.n_valid_sequences
    assert [*loaded.enumerate()] == [*space.enumerate()]


def test_coding_space_pickle_preserves_constraints():
    space = CodingSpace(
        'MIKEY',
        fixed_codons={2: 'ATC'},
        constraints=[ForbiddenMotifs(['AAAA']), MaxHomopolymer(4)],
        seed=8675309,
    )

    loaded = pickle.loads(pickle.dumps(space))

    assert len(loaded.constraints) == 2
    assert all(type(a) is type(b) for a, b in zip(loaded.constraints, space.constraints))
    assert loaded.n_valid_sequences == space.n_valid_sequences


def test_coding_space_pickle_preserves_future_sampling():
    space = CodingSpace('MIKEY', seed=8675309)

    _ = [space.sample() for _ in range(100)]

    loaded = pickle.loads(pickle.dumps(space))

    expected = [space.sample() for _ in range(1000)]
    observed = [loaded.sample() for _ in range(1000)]

    assert observed == expected


def test_coding_space_pickle_preserves_future_sampling_with_pins():
    space = CodingSpace('MIKEY', seed=8675309)
    space.pin_codons({2: 'ATC'})

    _ = [space.sample() for _ in range(100)]

    loaded = pickle.loads(pickle.dumps(space))

    expected = [space.sample() for _ in range(1000)]
    observed = [loaded.sample() for _ in range(1000)]

    assert observed == expected


def test_coding_space_pickle_preserves_initial_sampling():
    space = CodingSpace('MIKEY', seed=8675309)

    loaded = pickle.loads(pickle.dumps(space))

    expected = [space.sample() for _ in range(1000)]
    observed = [loaded.sample() for _ in range(1000)]

    assert observed == expected


def test_coding_space_save_load(tmp_path):
    path = tmp_path / 'space.pkl'

    space = CodingSpace('MIKEY', seed=8675309)
    space.save(path)

    loaded = CodingSpace.load(path)

    assert loaded.n_valid_sequences == space.n_valid_sequences
    assert [*loaded.enumerate()] == [*space.enumerate()]


def test_coding_space_save_load_preserves_random_state(tmp_path):
    path = tmp_path / 'space.pkl'

    space = CodingSpace('MIKEY', seed=8675309)

    _ = [space.sample() for _ in range(100)]

    space.save(path)
    loaded = CodingSpace.load(path)

    expected = [space.sample() for _ in range(1000)]
    observed = [loaded.sample() for _ in range(1000)]

    assert observed == expected


def test_coding_space_save_load_preserves_pins(tmp_path):
    path = tmp_path / 'space.pkl'

    space = CodingSpace('MIKEY', seed=8675309)
    space.pin_codons({2: 'ATC'})

    space.save(path)
    loaded = CodingSpace.load(path)

    assert loaded.view.pinned_codons == space.view.pinned_codons
    assert loaded.n_valid_sequences == space.n_valid_sequences
    assert [*loaded.enumerate()] == [*space.enumerate()]


def test_coding_space_exposes_graph_properties():
    cw = CodonWeights.ecoli()
    space = CodingSpace(
        'MIKEY',
        fixed_codons={2: 'ATC'},
        codon_weights=cw,
        context_l='AAA',
        context_r='CCC',
    )

    assert space.aa_seq == 'MIKEY'
    assert space.translation_table is space.view.translation_table
    assert space.codon_weights is space.view.codon_weights
    assert space.codon_weights is cw
    assert space.fixed_codons == space.view.fixed_codons
    assert space.context_l == 'AAA'
    assert space.context_r == 'CCC'


def test_coding_space_exposes_pins():
    space = CodingSpace('MIKEY')

    space.pin_codons({2: 'ATC'})
    assert space.pinned_codons == {2: ['ATC']}

    space.clear_pins()
    assert space.pinned_codons == {}


def test_coding_space_can_set_forbidden_motif_constraint():
    space = CodingSpace('MIKEY')
    constraint = ForbiddenMotifs(['AAA'])

    space.set_constraints(constraint)

    assert space.constraints == (constraint,)
    assert helper_get_banned_sequences_from_constraints(space) == {'AAA'}

    space.clear_constraints()

    assert space.constraints == ()
    assert helper_get_banned_sequences_from_constraints(space) == set()


def test_coding_space_can_set_homopolymer_constraint():
    space = CodingSpace('MIKEY')
    constraint = MaxHomopolymer(2)

    space.set_constraints(constraint)

    banned_sequences = helper_get_banned_sequences_from_constraints(space)
    assert banned_sequences == {'AAA', 'CCC', 'GGG', 'TTT'}

    space.clear_constraints()

    assert space.constraints == ()
    assert helper_get_banned_sequences_from_constraints(space) == set()


def test_coding_space_constraints_combine():
    space = CodingSpace('MIKEY')
    forbidden = ForbiddenMotifs(['GAG'])
    homopolymer = MaxHomopolymer(2)

    space.set_constraints([forbidden, homopolymer])

    banned_sequences = helper_get_banned_sequences_from_constraints(space)
    assert banned_sequences == {'AAA', 'CCC', 'GGG', 'TTT', 'GAG'}


def test_coding_space_can_set_codon_weights():
    space = CodingSpace('MIKEY')
    weights = CodonWeights.ecoli()

    space.set_codon_weights(weights)

    assert space.codon_weights is weights
    assert space.view.codon_weights is weights


def test_resolve_tables_defaults_to_dna():
    tt, cw = CodingSpace._resolve_tables(None, None, None)

    assert tt.rna is False
    assert 'ATG' in cw.weights
    assert 'AUG' not in cw.weights


def test_resolve_tables_uses_rna_flag():
    tt, cw = CodingSpace._resolve_tables(None, None, True)

    assert tt.rna is True
    assert 'AUG' in cw.weights
    assert 'ATG' not in cw.weights


def test_resolve_tables_builds_weights_from_table():
    tt_in = TranslationTable(rna=True)

    tt, cw = CodingSpace._resolve_tables(tt_in, None, None)

    assert tt is tt_in
    assert 'AUG' in cw.weights
    assert 'ATG' not in cw.weights


def test_resolve_tables_defaults_to_dna_when_only_weights_are_given():
    rna_table = TranslationTable(rna=True)
    cw_in = CodonWeights.uniform(table=rna_table)

    tt, cw = CodingSpace._resolve_tables(None, cw_in, None)

    assert tt.rna is False
    assert 'ATG' in cw.weights
    assert 'AUG' not in cw.weights


def test_resolve_tables_resolves_dna_weights_against_rna_table():
    dna_weights = CodonWeights.uniform()
    rna_table = TranslationTable(rna=True)

    tt, cw = CodingSpace._resolve_tables(rna_table, dna_weights, None)

    assert tt is rna_table
    assert cw is not dna_weights
    assert 'AUG' in cw.weights
    assert 'ATG' not in cw.weights


def test_resolve_tables_resolves_rna_weights_against_dna_table():
    rna_table = TranslationTable(rna=True)
    rna_weights = CodonWeights.uniform(table=rna_table)
    dna_table = TranslationTable()

    tt, cw = CodingSpace._resolve_tables(dna_table, rna_weights, None)

    assert tt is dna_table
    assert cw is not rna_weights
    assert 'ATG' in cw.weights
    assert 'AUG' not in cw.weights


def test_resolve_tables_rejects_incompatible_translation_table():
    weights = CodonWeights.uniform()

    with pytest.raises(ValueError):
        alt_table = TranslationTable(table_id=2)
        CodingSpace._resolve_tables(alt_table, weights, None)


def test_resolve_tables_rejects_rna_flag_table_mismatch():
    tt = TranslationTable(rna=False)

    with pytest.raises(ValueError, match='translation table'):
        CodingSpace._resolve_tables(tt, None, True)


def test_resolve_tables_does_not_modify_input_weights():
    weights = CodonWeights.uniform()

    rna_table = TranslationTable(rna=True)
    _, resolved = CodingSpace._resolve_tables(rna_table, weights, None)

    assert 'ATG' in weights.weights
    assert 'AUG' not in weights.weights

    assert 'AUG' in resolved.weights
    assert 'ATG' not in resolved.weights


def test_coding_space_rna_flag_creates_rna_space():
    space = CodingSpace('MKT', rna=True)

    assert space.translation_table.rna is True
    assert 'AUG' in space.codon_weights.weights
    assert 'ATG' not in space.codon_weights.weights
    assert 'U' in space.sample()
    assert 'T' not in space.sample()


def test_coding_space_rna_flag_normalises_inputs_to_rna():
    space = CodingSpace('M', fixed_codons={1: 'ATG'}, context_l='TTT', context_r='TAA', rna=True)

    assert space.fixed_codons == {1: ['AUG']}
    assert space.context_l == 'UUU'
    assert space.context_r == 'UAA'


def test_forbidden_motif_constraint_normalises_to_space_molecule_type():
    space = CodingSpace('M', rna=True, constraints=[ForbiddenMotifs(['ATG'])])
    assert space.n_valid_sequences == 0


def test_contains_normalises_dna_input_for_rna_space():
    space = CodingSpace('M', rna=True)
    assert space.contains('ATG')
    assert 'ATG' in space


def test_coding_space_accepts_constraints():
    constraint = ForbiddenMotifs(['GAA'])
    space = CodingSpace('E')
    space.set_constraints([constraint])

    assert space.constraints == (constraint,)
    assert set(space.enumerate()) == {'GAG'}


def test_coding_space_constraints_are_view_constraints():
    constraint = ForbiddenMotifs(['GAA'])
    space = CodingSpace('E', constraints=constraint)

    assert space.constraints is space.view.constraints

    replacement = ForbiddenMotifs(['GAG'])
    space.set_constraints(replacement)

    assert space.constraints is space.view.constraints
    assert space.constraints == (replacement,)


def test_coding_space_can_set_and_clear_constraints():
    space = CodingSpace('E')
    constraint = ForbiddenMotifs(['GAA'])

    space.set_constraints([constraint])
    assert set(space.enumerate()) == {'GAG'}

    space.clear_constraints()
    assert set(space.enumerate()) == {'GAA', 'GAG'}


def test_coding_space_can_set_multiple_constraints():
    space = CodingSpace('E')
    constraints = [
        ForbiddenMotifs(['AAA']),
        ForbiddenMotifs(['GAA']),
    ]

    space.set_constraints(constraints)

    assert space.constraints == tuple(constraints)
    assert set(space.enumerate()) == {'GAG'}


def test_coding_space_can_add_constraints():
    space = CodingSpace('MIKEY')
    constraint = ForbiddenMotifs(['AAA'])

    space.add_constraints(constraint)

    assert constraint in space.constraints
    assert constraint in space.view.constraints


def test_coding_space_accepts_single_constraint_at_initialisation():
    constraint = ForbiddenMotifs(['GAA'])

    space = CodingSpace('E', constraints=constraint)

    assert space.constraints == (constraint,)
    assert set(space.enumerate()) == {'GAG'}


def test_coding_space_accepts_constraints_at_initialisation():
    constraint = ForbiddenMotifs(['GAA'])

    space = CodingSpace('E', constraints=constraint)

    assert space.constraints == (constraint,)
    assert set(space.enumerate()) == {'GAG'}


class SafeConstraint(Constraint):

    def initial_state(self):
        return SAFE_STATE

    def advance(self, state, pos, choice):
        return SAFE_STATE

    def link(self, graph):
        pass


def test_coding_space_pickle_preserves_constraints():
    constraint = SafeConstraint()

    space = CodingSpace(
        'MIKEY',
        constraints=[constraint, ForbiddenMotifs(['CCC']), MaxHomopolymer(4)],
        seed=8675309,
    )

    loaded = pickle.loads(pickle.dumps(space))

    assert len(loaded.constraints) == 3
    assert isinstance(loaded.constraints[0], SafeConstraint)
    assert isinstance(loaded.constraints[1], ForbiddenMotifs)
    assert isinstance(loaded.constraints[2], MaxHomopolymer)
    assert loaded.n_valid_sequences == space.n_valid_sequences
    assert set(loaded.enumerate()) == set(space.enumerate())


def test_zero_weight_codons_are_valid_but_not_sampled():
    table = TranslationTable()
    data = {aa: {c: 1 for c in codons} for aa, codons in table.aa_to_codons.items()}
    data['F'] = {
        'TTT': 0.05,
        'TTC': 0.95,
    }

    weights = CodonWeights(data).threshold(0.1)
    space = CodingSpace('F', codon_weights=weights)

    assert space.n_valid_sequences == 2
    assert set(space.enumerate()) == {'TTT', 'TTC'}
    assert 'TTT' in space
    assert 'TTC' in space

    assert {space.sample() for _ in range(100)} == {'TTC'}


def test_sampling_raises_if_all_valid_sequences_have_zero_weight():
    table = TranslationTable()
    data = {aa: {c: 1 for c in codons} for aa, codons in table.aa_to_codons.items()}
    data['F'] = {
        'TTT': 0,
        'TTC': 1,
    }

    weights = CodonWeights(data)
    space = CodingSpace('F', fixed_codons={1: 'TTT'}, codon_weights=weights)

    assert space.n_valid_sequences == 1
    assert set(space.enumerate()) == {'TTT'}

    with pytest.raises(ValueError):
        space.sample()


def test_restricted_weights_remove_zero_weight_codons_from_space():
    table = TranslationTable()
    data = {aa: {c: 1 for c in codons} for aa, codons in table.aa_to_codons.items()}
    data['F'] = {
        'TTT': 0.05,
        'TTC': 0.95,
    }

    weights = CodonWeights(data).threshold(0.1)
    table, weights = weights.restrict(table)

    space = CodingSpace('F', translation_table=table, codon_weights=weights)

    assert space.n_valid_sequences == 1
    assert set(space.enumerate()) == {'TTC'}
    assert 'TTC' in space
    assert 'TTT' not in space


def test_coding_space_codon_options():
    space = CodingSpace('MIF', fixed_codons={2: ['ATT', 'ATC']})
    space.pin_codons({3: 'TTC'})

    assert space.codon_options[1] == ('ATG',)
    assert set(space.codon_options[2]) == {'ATT', 'ATC'}
    assert space.codon_options[3] == ('TTC',)

def test_coding_space_copy():
    space = CodingSpace('MIKEY', seed=8675309)
    space.pin_codons({2: 'ATC'})

    copied = space.copy()

    assert copied is not space
    assert copied.view is not space.view
    assert copied.pinned_codons == space.pinned_codons
    assert copied.n_valid_sequences == space.n_valid_sequences

    copied.clear_pins()

    assert copied.pinned_codons == {}
    assert space.pinned_codons == {2: ['ATC']}
