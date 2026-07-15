import pickle

import pytest

from codeine.translation.tables import TranslationTable
from codeine.translation.weights import CodonWeights


def make_weights_data(table):
    return {
        aa: {codon: 1.0 for codon in codons}
        for aa, codons in table.aa_to_codons.items()
    }


def test_uniform_weights_sum_to_one_per_amino_acid():
    table = TranslationTable()
    weights = CodonWeights.uniform(table)

    for codons in table.aa_to_codons.values():
        total = sum(weights[codon] for codon in codons)
        assert total == pytest.approx(1.0)


def test_uniform_weights_contains_all_table_codons():
    table = TranslationTable()
    weights = CodonWeights.uniform(table)

    assert set(weights.weights.keys()) == set(table.codons_to_aa.keys())


def test_getitem_returns_codon_weight():
    weights = CodonWeights.uniform()

    assert weights['GCT'] == pytest.approx(0.25)


def test_by_aa_returns_all_codon_weights():
    table = TranslationTable()
    weights = CodonWeights.uniform(table)

    assert set(weights.by_aa('A').keys()) == set(table.aa_to_codons['A'])


def test_uniform_weights_are_equal_within_each_amino_acid():
    table = TranslationTable()
    weights = CodonWeights.uniform(table)

    for codons in table.aa_to_codons.values():
        expected = 1.0 / len(codons)

        for codon in codons:
            assert weights[codon] == pytest.approx(expected)


def test_weights_are_normalised_per_amino_acid():
    table = TranslationTable()

    data = {
        aa: {codon: i + 1 for i, codon in enumerate(codons)}
        for aa, codons in table.aa_to_codons.items()
    }

    weights = CodonWeights(data)

    for codons in table.aa_to_codons.values():
        total = sum(weights[codon] for codon in codons)
        assert total == pytest.approx(1.0)


def test_constructor_preserves_dna_codons():
    table = TranslationTable()
    data = make_weights_data(table)

    weights = CodonWeights(data)

    assert 'GCT' in weights.weights
    assert 'GCU' not in weights.weights


def test_constructor_preserves_rna_codons():
    table = TranslationTable(rna=True)
    data = make_weights_data(table)

    weights = CodonWeights(data)

    assert 'GCU' in weights.weights
    assert 'GCT' not in weights.weights


def test_resolve_converts_dna_codons_to_rna():
    dna_table = TranslationTable()
    rna_table = TranslationTable(rna=True)

    weights = CodonWeights(make_weights_data(dna_table))
    resolved = weights.for_table(rna_table)

    assert all('T' not in codon for codon in resolved.weights.keys())
    assert 'GCU' in resolved.weights
    assert 'GCT' not in resolved.weights


def test_resolve_converts_rna_codons_to_dna():
    dna_table = TranslationTable()
    rna_table = TranslationTable(rna=True)

    weights = CodonWeights(make_weights_data(rna_table))
    resolved = weights.for_table(dna_table)

    assert all('U' not in codon for codon in resolved.weights.keys())
    assert 'GCT' in resolved.weights
    assert 'GCU' not in resolved.weights


def test_resolve_returns_new_codon_weights_object():
    dna_table = TranslationTable()
    rna_table = TranslationTable(rna=True)

    weights = CodonWeights(make_weights_data(dna_table))
    resolved = weights.for_table(rna_table)

    assert resolved is not weights
    assert type(resolved) is CodonWeights
    assert 'GCT' in weights.weights
    assert 'GCU' in resolved.weights


def test_resolve_preserves_weights():
    dna_table = TranslationTable()
    rna_table = TranslationTable(rna=True)

    data = {
        aa: {codon: i + 1 for i, codon in enumerate(codons)}
        for aa, codons in dna_table.aa_to_codons.items()
    }

    weights = CodonWeights(data)
    resolved = weights.for_table(rna_table)

    for dna_codon, weight in weights.weights.items():
        rna_codon = dna_codon.replace('T', 'U')
        assert resolved[rna_codon] == pytest.approx(weight)


def test_resolve_rejects_incompatible_translation_table():
    weights = CodonWeights.uniform()

    with pytest.raises(ValueError):
        weights.for_table(TranslationTable(table_id=2))


def test_uniform_can_use_rna_table():
    table = TranslationTable(rna=True)
    weights = CodonWeights.uniform(table)

    assert 'GCU' in weights.weights
    assert 'GCT' not in weights.weights


def test_uniform_defaults_to_dna():
    weights = CodonWeights.uniform()

    assert 'ATG' in weights.weights
    assert 'AUG' not in weights.weights


def test_negative_weight_raises_value_error():
    table = TranslationTable()
    data = make_weights_data(table)
    data['A']['GCT'] = -1.0

    with pytest.raises(ValueError, match='cannot be negative'):
        CodonWeights(data)


def test_zero_total_weight_raises_value_error():
    table = TranslationTable()
    data = make_weights_data(table)
    data['A'] = {codon: 0.0 for codon in data['A']}

    with pytest.raises(ValueError, match='must sum to > 0'):
        CodonWeights(data)


def test_duplicate_amino_acid_after_normalisation_raises_value_error():
    data = {
        'A': {'GCT': 1.0},
        'a': {'GCC': 1.0},
    }

    with pytest.raises(ValueError, match='Duplicate weights for amino acid A'):
        CodonWeights(data)


def test_duplicate_codon_within_amino_acid_raises_value_error():
    data = {
        'A': {
            'GCT': 1.0,
            'gct': 1.0,
        },
    }

    with pytest.raises(ValueError, match='Duplicate weight for codon GCT'):
        CodonWeights(data)


def test_duplicate_codon_between_amino_acids_raises_value_error():
    data = {
        'A': {'GCT': 1.0},
        'R': {'GCT': 1.0},
    }

    with pytest.raises(ValueError, match='Duplicate weight for codon GCT'):
        CodonWeights(data)


def test_weights_is_immutable():
    weights = CodonWeights.uniform()

    with pytest.raises(AttributeError):
        weights.weights = {}

    with pytest.raises(AttributeError):
        weights.chicken = 'duck'

    with pytest.raises(TypeError):
        weights.weights['GCT'] = 0.5


def test_getitem_is_strict():
    weights = CodonWeights.uniform()

    with pytest.raises(KeyError):
        _ = weights['GCU']

    with pytest.raises(KeyError):
        _ = weights['GCt']

    with pytest.raises(KeyError):
        _ = weights['gct']

    weights_rna = weights.for_table(TranslationTable(rna=True))

    with pytest.raises(KeyError):
        _ = weights_rna['GCT']

    with pytest.raises(KeyError):
        _ = weights_rna['GCu']

    with pytest.raises(KeyError):
        _ = weights_rna['gcu']


constructors = [
    CodonWeights.uniform,
    CodonWeights.ecoli,
    CodonWeights.yeast,
    CodonWeights.human,
    CodonWeights.mouse,
    CodonWeights.arabidopsis,
    CodonWeights.drosophila,
]


@pytest.mark.parametrize('rna', (True, False))
@pytest.mark.parametrize('constructor', constructors)
def test_builtin_weights_are_valid(constructor, rna):
    table = TranslationTable(rna=rna)
    weights = constructor().for_table(table)

    assert set(weights.weights.keys()) == set(table.codons_to_aa.keys())

    for codons in table.aa_to_codons.values():
        total = sum(weights[codon] for codon in codons)
        assert total == pytest.approx(1.0)


@pytest.mark.parametrize('constructor', constructors)
def test_builtin_dna_and_rna_weights_match_after_codon_conversion(constructor):
    dna_weights = constructor().for_table(TranslationTable())
    rna_weights = constructor().for_table(TranslationTable(rna=True))

    for dna_codon, dna_weight in dna_weights.weights.items():
        rna_codon = dna_codon.replace('T', 'U')
        assert rna_weights[rna_codon] == pytest.approx(dna_weight)


@pytest.mark.parametrize('constructor', constructors)
def test_codon_weights_pickle(constructor):
    weights = constructor()
    loaded = pickle.loads(pickle.dumps(weights))

    assert type(loaded) is CodonWeights
    assert loaded.aa_to_codons == weights.aa_to_codons
    assert loaded.weights == weights.weights
    assert loaded['ATG'] == weights['ATG']
    assert loaded.by_aa('K') == weights.by_aa('K')


def test_resolved_rna_codon_weights_pickle():
    weights = CodonWeights.uniform().for_table(
        TranslationTable(rna=True)
    )

    loaded = pickle.loads(pickle.dumps(weights))

    assert type(loaded) is CodonWeights
    assert loaded.aa_to_codons == weights.aa_to_codons
    assert loaded.weights == weights.weights
    assert loaded['AUG'] == weights['AUG']


def test_codon_weights_pickle_preserves_immutability():
    weights = CodonWeights.uniform()
    loaded = pickle.loads(pickle.dumps(weights))

    with pytest.raises(AttributeError):
        loaded.chicken = 'duck'


def test_by_aa_is_strict():
    weights = CodonWeights.uniform()

    for x in ('atg', 'u', 'x', 'AMINO', 'MIKEY', '#'):
        with pytest.raises(KeyError):
            weights.by_aa(x)


def test_by_aa_returns_copy():
    weights = CodonWeights.uniform()
    aa_weights = weights.by_aa('A')

    aa_weights['GCT'] = 0.9

    assert weights['GCT'] == pytest.approx(0.25)


def test_missing_amino_acid_raises_value_error_when_resolved():
    table = TranslationTable()
    data = make_weights_data(table)
    del data['A']

    weights = CodonWeights(data)

    with pytest.raises(
        ValueError,
        match='Missing weights for amino acid',
    ):
        weights.for_table(table)


def test_extra_amino_acid_raises_value_error_when_resolved():
    table = TranslationTable()
    data = make_weights_data(table)
    data['?'] = {'NNN': 1.0}

    weights = CodonWeights(data)

    with pytest.raises(ValueError, match='Unknown amino acid'):
        weights.for_table(table)


def test_missing_codon_raises_value_error_when_resolved():
    table = TranslationTable()
    data = make_weights_data(table)
    del data['A']['GCT']

    weights = CodonWeights(data)

    with pytest.raises(ValueError, match='Missing weights for codon'):
        weights.for_table(table)


def test_zero_weights_are_allowed():
    weights = CodonWeights(
        {
            'F': {
                'TTT': 0.0,
                'TTC': 1.0,
            },
        }
    )

    assert weights['TTT'] == 0.0
    assert weights['TTC'] == 1.0


def test_threshold_sets_weights_below_minimum_to_zero():
    weights = CodonWeights(
        {
            'F': {
                'TTT': 9,
                'TTC': 91,
            },
        }
    )

    thresholded = weights.threshold(0.1)

    assert thresholded['TTT'] == 0.0
    assert thresholded['TTC'] == 1.0


def test_threshold_keeps_weights_equal_to_minimum():
    weights = CodonWeights(
        {
            'F': {
                'TTT': 0.1,
                'TTC': 0.9,
            },
        }
    )

    thresholded = weights.threshold(0.1)

    assert thresholded is weights
    assert thresholded['TTT'] == pytest.approx(0.1)
    assert thresholded['TTC'] == pytest.approx(0.9)


def test_threshold_renormalises_remaining_weights():
    weights = CodonWeights(
        {
            'I': {
                'ATT': 0.05,
                'ATC': 0.25,
                'ATA': 0.70,
            },
        }
    )

    thresholded = weights.threshold(0.1)

    assert thresholded['ATT'] == 0.0
    assert thresholded['ATC'] == pytest.approx(0.25 / 0.95)
    assert thresholded['ATA'] == pytest.approx(0.70 / 0.95)


def test_threshold_does_not_modify_original_weights():
    weights = CodonWeights(
        {
            'F': {
                'TTT': 0.05,
                'TTC': 0.95,
            },
        }
    )

    thresholded = weights.threshold(0.1)

    assert weights['TTT'] == pytest.approx(0.05)
    assert weights['TTC'] == pytest.approx(0.95)
    assert thresholded is not weights


@pytest.mark.parametrize('min_weight', (-0.1, 1.1))
def test_threshold_rejects_invalid_minimum_weight(min_weight):
    weights = CodonWeights.uniform()

    with pytest.raises(ValueError, match='between 0 and 1'):
        weights.threshold(min_weight)


def test_restrict_removes_zero_weight_codons_from_table_and_weights():
    table = TranslationTable()
    data = make_weights_data(table)
    data['F'] = {
        'TTT': 0.0,
        'TTC': 1.0,
    }

    weights = CodonWeights(data)
    restricted_table, restricted_weights = weights.restrict(table)

    assert 'TTT' not in restricted_table.codons_to_aa
    assert 'TTT' not in restricted_weights.weights

    assert restricted_table['TTC'] == 'F'
    assert restricted_weights['TTC'] == 1.0


def test_restrict_returns_compatible_table_and_weights():
    table = TranslationTable()
    data = make_weights_data(table)
    data['F'] = {
        'TTT': 0.0,
        'TTC': 1.0,
    }

    restricted_table, restricted_weights = CodonWeights(data).restrict(table)

    assert set(restricted_weights.weights) == set(restricted_table.codons_to_aa)
    assert restricted_weights.for_table(restricted_table) is restricted_weights


def test_restrict_preserves_relative_positive_weights():
    table = TranslationTable()
    data = make_weights_data(table)
    data['I'] = {
        'ATT': 0.0,
        'ATC': 1.0,
        'ATA': 3.0,
    }

    _, restricted_weights = CodonWeights(data).restrict(table)

    assert restricted_weights['ATC'] == pytest.approx(0.25)
    assert restricted_weights['ATA'] == pytest.approx(0.75)


def test_restrict_returns_original_objects_when_no_weights_are_zero():
    table = TranslationTable()
    weights = CodonWeights.uniform(table)

    restricted_table, restricted_weights = weights.restrict(table)

    assert restricted_table is table
    assert restricted_weights is weights


def test_restrict_resolves_weights_to_table_molecule_type():
    dna_table = TranslationTable()
    rna_table = TranslationTable(rna=True)

    data = make_weights_data(dna_table)
    data['F'] = {
        'TTT': 0.0,
        'TTC': 1.0,
    }

    restricted_table, restricted_weights = CodonWeights(data).restrict(rna_table)

    assert restricted_table.rna is True
    assert 'UUU' not in restricted_table.codons_to_aa
    assert 'UUU' not in restricted_weights.weights
    assert 'UUC' in restricted_table.codons_to_aa
    assert 'UUC' in restricted_weights.weights


def test_restricted_table_name_lists_omitted_codons():
    table = TranslationTable()
    data = make_weights_data(table)
    data['F'] = {
        'TTT': 0.0,
        'TTC': 1.0,
    }

    restricted_table, _ = CodonWeights(data).restrict(table)

    assert table.name in restricted_table.name
    assert 'TTT' in restricted_table.name


def test_restricted_table_name_summarises_many_omitted_codons():
    table = TranslationTable()
    data = make_weights_data(table)

    omitted_codons = []
    for aa, codons in table.aa_to_codons.items():
        for codon in codons[:-1]:
            data[aa][codon] = 0.0
            omitted_codons.append(codon)

    restricted_table, _ = CodonWeights(data).restrict(table)

    assert f'omitting {len(omitted_codons)} zero-weight codons' in restricted_table.name
