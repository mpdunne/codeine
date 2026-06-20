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
    assert set(weights.weights) == set(table.codons_to_aa)


def test_getitem_returns_codon_weight():
    weights = CodonWeights.uniform()
    assert weights['GCT'] == pytest.approx(0.25)


def test_by_aa_returns_all_codon_weights():
    table = TranslationTable()
    weights = CodonWeights.uniform(table)
    assert set(weights.by_aa('A')) == set(table.aa_to_codons['A'])


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


def test_rna_weights_convert_codons_to_rna():
    table = TranslationTable()
    data = make_weights_data(table)
    weights = CodonWeights(data, rna=True)
    assert weights.rna is True
    assert all('T' not in codon for codon in weights.weights)
    assert 'GCU' in weights.weights
    assert 'GCT' not in weights.weights


def test_dna_weights_convert_codons_to_dna():
    table = TranslationTable(rna=True)
    data = make_weights_data(table)
    weights = CodonWeights(data, rna=False)
    assert weights.rna is False
    assert all('U' not in codon for codon in weights.weights)
    assert 'GCT' in weights.weights
    assert 'GCU' not in weights.weights


def test_uniform_can_use_rna_flag_without_table():
    weights = CodonWeights.uniform(rna=True)
    assert weights.rna is True
    assert 'GCU' in weights.weights
    assert 'GCT' not in weights.weights


def test_uniform_uses_table_dna_if_rna_is_not_given():
    table = TranslationTable(rna=True)
    weights = CodonWeights.uniform(table)
    assert weights.rna
    assert 'GCU' in weights.weights
    assert 'GCT' not in weights.weights


def test_uniform_rna_flag_can_override_table_molecule_type():
    table = TranslationTable(rna=False)
    weights = CodonWeights.uniform(table, rna=True)
    assert weights.rna is True
    assert 'GCU' in weights.weights
    assert 'GCT' not in weights.weights


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

    weights_rna = CodonWeights.uniform(rna=True)

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
    weights = constructor(rna=rna)

    assert weights.rna is rna
    assert set(weights.weights) == set(table.codons_to_aa)

    for codons in table.aa_to_codons.values():
        total = sum(weights[codon] for codon in codons)
        assert total == pytest.approx(1.0)


@pytest.mark.parametrize('constructor', constructors)
def test_builtin_dna_and_rna_weights_match_after_codon_conversion(constructor):
    dna_weights = constructor()
    rna_weights = constructor(rna=True)

    for dna_codon, dna_weight in dna_weights.weights.items():
        rna_codon = dna_codon.replace('T', 'U')
        assert rna_weights[rna_codon] == pytest.approx(dna_weight)


@pytest.mark.parametrize('rna', (True, False))
@pytest.mark.parametrize('constructor', constructors)
def test_codon_weights_pickle(constructor, rna):
    weights = constructor(rna=rna)
    loaded = pickle.loads(pickle.dumps(weights))
    assert type(loaded) is CodonWeights
    assert loaded.rna == weights.rna
    assert loaded.aa_to_codons == weights.aa_to_codons
    assert loaded.weights == weights.weights

    test_codon = 'AUG' if rna else 'ATG'
    assert loaded[test_codon] == weights[test_codon]
    assert loaded.by_aa('K') == weights.by_aa('K')


def test_codon_weights_pickle_preserves_immutability():
    weights = CodonWeights.uniform()
    loaded = pickle.loads(pickle.dumps(weights))
    with pytest.raises(AttributeError):
        loaded.rna = True


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


def test_missing_amino_acid_raises_value_error():
    table = TranslationTable()
    data = make_weights_data(table)
    del data['A']

    with pytest.raises(ValueError):
        CodonWeights(data)


def test_missing_codon_raises_value_error():
    table = TranslationTable()
    data = make_weights_data(table)
    del data['A']['GCT']

    with pytest.raises(ValueError):
        CodonWeights(data)


def test_extra_codon_raises_value_error():
    table = TranslationTable()
    data = make_weights_data(table)
    data['A']['AAA'] = 1.0

    with pytest.raises(ValueError):
        CodonWeights(data)
