import pytest

from codeine.translation.tables import TranslationTable
from codeine.translation.weights import CodonWeights


def test_uniform_weights_sum_to_one_per_amino_acid():
    table = TranslationTable()
    weights = CodonWeights.uniform(table)

    for aa, codons in table.aa_to_codons.items():
        total = sum(weights[codon] for codon in codons)
        assert total == pytest.approx(1.0)


def test_uniform_weights_contains_all_table_codons():
    table = TranslationTable()
    weights = CodonWeights.uniform(table)

    assert set(weights.weights) == set(table.codons_to_aa)


def test_getitem_returns_codon_weight():
    weights = CodonWeights.uniform()

    assert weights['GCT'] == pytest.approx(0.25)


def test_for_amino_acid_returns_all_codon_weights():
    table = TranslationTable()
    weights = CodonWeights.uniform(table)

    assert set(weights.for_amino_acid('A')) == set(table.aa_to_codons['A'])


def test_uniform_weights_are_equal_within_each_amino_acid():
    table = TranslationTable()
    weights = CodonWeights.uniform(table)

    for aa, codons in table.aa_to_codons.items():
        expected = 1.0 / len(codons)

        for codon in codons:
            assert weights[codon] == pytest.approx(expected)


def test_weights_are_normalised_per_amino_acid():
    table = TranslationTable()

    raw = {
        aa: {codon: i + 1 for i, codon in enumerate(codons)}
        for aa, codons in table.aa_to_codons.items()
    }

    weights = CodonWeights(raw, table=table)

    for aa, codons in table.aa_to_codons.items():
        total = sum(weights[codon] for codon in codons)
        assert total == pytest.approx(1.0)


def test_rna_table_normalises_codons_to_rna():
    table = TranslationTable(rna=True)

    raw = {
        aa: {codon.replace('U', 'T'): 1.0 for codon in codons}
        for aa, codons in table.aa_to_codons.items()
    }

    weights = CodonWeights(raw, table=table)

    for codon in weights.weights:
        assert 'T' not in codon


def test_dna_table_normalises_codons_to_dna():
    table = TranslationTable(rna=False)

    raw = {
        aa: {codon.replace('T', 'U'): 1.0 for codon in codons}
        for aa, codons in table.aa_to_codons.items()
    }

    weights = CodonWeights(raw, table=table)

    for codon in weights.weights:
        assert 'U' not in codon


def test_missing_amino_acid_raises_value_error():
    table = TranslationTable()

    raw = {
        aa: {codon: 1.0 for codon in codons}
        for aa, codons in table.aa_to_codons.items()
    }
    raw.pop('A')

    with pytest.raises(ValueError, match='Missing weights'):
        CodonWeights(raw, table=table)


def test_extra_amino_acid_raises_value_error():
    table = TranslationTable()

    raw = {
        aa: {codon: 1.0 for codon in codons}
        for aa, codons in table.aa_to_codons.items()
    }
    raw['B'] = {'BBB': 1.0}

    with pytest.raises(ValueError, match='Unexpected amino acids'):
        CodonWeights(raw, table=table)


def test_missing_codon_raises_value_error():
    table = TranslationTable()

    raw = {
        aa: {codon: 1.0 for codon in codons}
        for aa, codons in table.aa_to_codons.items()
    }

    aa = 'A'
    codon = next(iter(raw[aa]))
    raw[aa].pop(codon)

    with pytest.raises(ValueError, match='Missing codon weights'):
        CodonWeights(raw, table=table)


def test_extra_codon_raises_value_error():
    table = TranslationTable()

    raw = {
        aa: {codon: 1.0 for codon in codons}
        for aa, codons in table.aa_to_codons.items()
    }
    raw['A']['TTT'] = 1.0

    with pytest.raises(ValueError, match='Unexpected codons'):
        CodonWeights(raw, table=table)


def test_negative_weight_raises_value_error():
    table = TranslationTable()

    raw = {
        aa: {codon: 1.0 for codon in codons}
        for aa, codons in table.aa_to_codons.items()
    }

    raw['A'][next(iter(raw['A']))] = -1.0

    with pytest.raises(ValueError, match='cannot be negative'):
        CodonWeights(raw, table=table)


def test_zero_total_weight_raises_value_error():
    table = TranslationTable()

    raw = {
        aa: {codon: 1.0 for codon in codons}
        for aa, codons in table.aa_to_codons.items()
    }

    raw['A'] = {codon: 0.0 for codon in raw['A']}

    with pytest.raises(ValueError, match='must sum to > 0'):
        CodonWeights(raw, table=table)


def test_from_codon_weights_stores_flat_weights():
    table = TranslationTable()

    flat = {
        codon: 1.0
        for codon in table.codons_to_aa
    }

    weights = CodonWeights.from_codon_weights(flat, table=table)

    assert set(weights.weights) == set(table.codons_to_aa)


def test_from_codon_weights_normalises_per_amino_acid():
    table = TranslationTable()

    flat = {
        codon: i + 1
        for i, codon in enumerate(table.codons_to_aa)
    }

    weights = CodonWeights.from_codon_weights(flat, table=table)

    for aa, codons in table.aa_to_codons.items():
        total = sum(weights[codon] for codon in codons)
        assert total == pytest.approx(1.0)


def test_from_codon_weights_unknown_codon_raises_key_error():
    table = TranslationTable()

    flat = {
        codon: 1.0
        for codon in table.codons_to_aa
    }
    flat['XXX'] = 1.0

    with pytest.raises(KeyError):
        CodonWeights.from_codon_weights(flat, table=table)


def test_codon_weights_is_immutable():
    weights = CodonWeights.uniform()

    with pytest.raises(AttributeError):
        weights.weights = {}


def test_weights_mapping_is_immutable():
    weights = CodonWeights.uniform()

    with pytest.raises(TypeError):
        weights.weights['GCT'] = 0.5


def test_getitem_normalises_input_codon():
    table = TranslationTable(rna=False)
    weights = CodonWeights.uniform(table)

    assert weights['GCU'] == weights['GCT']


def test_repr():
    table = TranslationTable(table_id=1, rna=False)
    weights = CodonWeights.uniform(table)

    assert repr(weights) == 'CodonWeights(table=TranslationTable(table_id=1, rna=False))'
