import pytest

from codeine.translation.tables import CodonTable


def test_default_codon_table_works_as_expected():
    ct = CodonTable()

    assert ct.codons_to_aa["ATG"] == "M"
    assert ct.codons_to_aa["CCC"] == "P"
    assert ct.codons_to_aa["GAC"] == "D"

    assert ct.aa_to_codons["M"] == ["ATG"]
    assert set(ct.aa_to_codons["P"]) == {"CCT", "CCC", "CCA", "CCG"}
    assert set(ct.aa_to_codons["D"]) == {"GAT", "GAC"}

    # Should be no stop codons
    assert "TAA" not in ct.codons_to_aa
    assert "TAG" not in ct.codons_to_aa
    assert "TGA" not in ct.codons_to_aa


def test_all_forward_codons_appear_in_reverse_table():
    ct = CodonTable()
    reverse_codons = {codon for codons in ct.aa_to_codons.values() for codon in codons}
    assert reverse_codons == set(ct.codons_to_aa)


def test_codon_probabilities_sum_to_one_per_amino_acid():
    ct = CodonTable()
    for aa, codon_probs in ct.codon_probabilities.items():
        assert sum(codon_probs.values()) == pytest.approx(1.0)


def test_codon_probabilities_are_uniform_per_amino_acid():
    ct = CodonTable()
    for aa, codons in ct.aa_to_codons.items():
        expected = 1 / len(codons)
        for codon in codons:
            assert ct.codon_probabilities[aa][codon] == pytest.approx(expected)


def test_codon_probabilities_match_reverse_table():
    ct = CodonTable()
    for aa, codon_probs in ct.codon_probabilities.items():
        assert set(codon_probs) == set(ct.aa_to_codons[aa])