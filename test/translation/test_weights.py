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

    weights = CodonWeights(data, table=table)

    for codons in table.aa_to_codons.values():
        total = sum(weights[codon] for codon in codons)
        assert total == pytest.approx(1.0)


def test_rna_table_requires_rna_codons():
    table = TranslationTable(rna=True)
    data = make_weights_data(table)

    weights = CodonWeights(data, table=table)

    assert all('T' not in codon for codon in weights.weights)


def test_dna_table_requires_dna_codons():
    table = TranslationTable(rna=False)
    data = make_weights_data(table)

    weights = CodonWeights(data, table=table)

    assert all('U' not in codon for codon in weights.weights)


def test_rna_table_rejects_dna_codons():
    table = TranslationTable(rna=True)

    data = {
        aa: {codon.replace('U', 'T'): 1.0 for codon in codons}
        for aa, codons in table.aa_to_codons.items()
    }

    with pytest.raises(ValueError, match='Missing codon weights|Unexpected codons'):
        CodonWeights(data, table=table)


def test_dna_table_rejects_rna_codons():
    table = TranslationTable(rna=False)

    data = {
        aa: {codon.replace('T', 'U'): 1.0 for codon in codons}
        for aa, codons in table.aa_to_codons.items()
    }

    with pytest.raises(ValueError, match='Missing codon weights|Unexpected codons'):
        CodonWeights(data, table=table)


def test_missing_amino_acid_raises_value_error():
    table = TranslationTable()
    data = make_weights_data(table)
    data.pop('A')
    with pytest.raises(ValueError, match='Missing weights'):
        CodonWeights(data, table=table)


def test_extra_amino_acid_raises_value_error():
    table = TranslationTable()
    data = make_weights_data(table)
    data['B'] = {'BBB': 1.0}
    with pytest.raises(ValueError, match='Unexpected amino acids'):
        CodonWeights(data, table=table)


def test_missing_codon_raises_value_error():
    table = TranslationTable()
    data = make_weights_data(table)
    data['A'].pop('GCT')
    with pytest.raises(ValueError, match='Missing codon weights'):
        CodonWeights(data, table=table)


def test_extra_codon_raises_value_error():
    table = TranslationTable()
    data = make_weights_data(table)
    data['A']['TTT'] = 1.0
    with pytest.raises(ValueError, match='Unexpected codons'):
        CodonWeights(data, table=table)


def test_negative_weight_raises_value_error():
    table = TranslationTable()
    data = make_weights_data(table)
    data['A']['GCT'] = -1.0
    with pytest.raises(ValueError, match='cannot be negative'):
        CodonWeights(data, table=table)


def test_zero_total_weight_raises_value_error():
    table = TranslationTable()
    data = make_weights_data(table)
    data['A'] = {codon: 0.0 for codon in data['A']}
    with pytest.raises(ValueError, match='must sum to > 0'):
        CodonWeights(data, table=table)


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

    table_rna = TranslationTable(rna=True)
    weights_rna = CodonWeights.uniform(table_rna)

    with pytest.raises(KeyError):
        _ = weights_rna['GCT']

    with pytest.raises(KeyError):
        _ = weights_rna['GCu']

    with pytest.raises(KeyError):
        _ = weights_rna['gcu']

