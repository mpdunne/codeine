import pytest

from codeine.translation.tables import CodonTable


def test_default_codon_table_works_as_expected():
    ct = CodonTable()

    assert ct.rna is False

    assert ct.codons_to_aa["ATG"] == "M"
    assert ct.codons_to_aa["CCC"] == "P"
    assert ct.codons_to_aa["GAC"] == "D"

    assert ct.aa_to_codons["M"] == ("ATG",)
    assert set(ct.aa_to_codons["P"]) == {"CCT", "CCC", "CCA", "CCG"}
    assert set(ct.aa_to_codons["D"]) == {"GAT", "GAC"}

    assert "TAA" not in ct.codons_to_aa
    assert "TAG" not in ct.codons_to_aa
    assert "TGA" not in ct.codons_to_aa


def test_rna_codon_table_works_as_expected():
    ct = CodonTable(rna=True)

    assert ct.rna is True

    assert ct.codons_to_aa["AUG"] == "M"
    assert ct.codons_to_aa["CCC"] == "P"
    assert ct.codons_to_aa["GAC"] == "D"

    assert "ATG" not in ct.codons_to_aa

    assert ct.aa_to_codons["M"] == ("AUG",)
    assert set(ct.aa_to_codons["P"]) == {"CCU", "CCC", "CCA", "CCG"}
    assert set(ct.aa_to_codons["D"]) == {"GAU", "GAC"}

    assert "UAA" not in ct.codons_to_aa
    assert "UAG" not in ct.codons_to_aa
    assert "UGA" not in ct.codons_to_aa


def test_all_forward_codons_appear_in_reverse_table():
    ct = CodonTable()
    reverse_codons = {
        codon
        for codons in ct.aa_to_codons.values()
        for codon in codons
    }
    assert reverse_codons == set(ct.codons_to_aa)


def test_all_forward_rna_codons_appear_in_reverse_table():
    ct = CodonTable(rna=True)
    reverse_codons = {
        codon
        for codons in ct.aa_to_codons.values()
        for codon in codons
    }
    assert reverse_codons == set(ct.codons_to_aa)


@pytest.mark.parametrize("rna", [False, True])
def test_codon_probabilities_sum_to_one_per_amino_acid(rna):
    ct = CodonTable(rna=rna)
    for codon_probs in ct.codon_probabilities.values():
        assert sum(codon_probs.values()) == pytest.approx(1.0)


@pytest.mark.parametrize("rna", [False, True])
def test_codon_probabilities_are_uniform_per_amino_acid(rna):
    ct = CodonTable(rna=rna)
    for aa, codons in ct.aa_to_codons.items():
        expected = 1 / len(codons)
        for codon in codons:
            assert ct.codon_probabilities[aa][codon] == pytest.approx(expected)


@pytest.mark.parametrize("rna", [False, True])
def test_codon_probabilities_match_reverse_table(rna):
    ct = CodonTable(rna=rna)
    for aa, codon_probs in ct.codon_probabilities.items():
        assert set(codon_probs) == set(ct.aa_to_codons[aa])


def test_dna_table_contains_no_rna_codons():
    ct = CodonTable()
    assert all("U" not in codon for codon in ct.codons_to_aa)
    assert all(
        "U" not in codon
        for codons in ct.aa_to_codons.values()
        for codon in codons
    )


def test_rna_table_contains_no_dna_thymine_codons():
    ct = CodonTable(rna=True)
    assert all("T" not in codon for codon in ct.codons_to_aa)
    assert all(
        "T" not in codon
        for codons in ct.aa_to_codons.values()
        for codon in codons
    )


def test_codon_tables_are_immutable():
    ct = CodonTable()

    with pytest.raises(TypeError):
        ct.codons_to_aa["AAA"] = "X"

    with pytest.raises(TypeError):
        ct.aa_to_codons["X"] = ("AAA",)

    with pytest.raises(TypeError):
        ct.codon_probabilities["X"] = {"AAA": 1.0}


def test_aa_to_codons_values_are_immutable():
    ct = CodonTable()
    assert isinstance(ct.aa_to_codons["M"], tuple)
    with pytest.raises(AttributeError):
        ct.aa_to_codons["M"].append("XXX")


def test_nested_codon_probability_tables_are_immutable():
    ct = CodonTable()
    with pytest.raises(TypeError):
        ct.codon_probabilities["M"]["ATG"] = 0.5


def test_rna_attribute_is_read_only():
    ct = CodonTable(rna=True)
    with pytest.raises(AttributeError):
        ct.rna = False
