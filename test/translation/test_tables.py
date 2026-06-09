import pytest

from codeine.translation.tables import TranslationTable


def test_default_codon_table_works_as_expected():
    tt = TranslationTable()

    assert tt.rna is False

    assert tt.codons_to_aa["ATG"] == "M"
    assert tt.codons_to_aa["CCC"] == "P"
    assert tt.codons_to_aa["GAC"] == "D"

    assert tt.aa_to_codons["M"] == ("ATG",)
    assert set(tt.aa_to_codons["P"]) == {"CCT", "CCC", "CCA", "CCG"}
    assert set(tt.aa_to_codons["D"]) == {"GAT", "GAC"}


def test_rna_codon_table_works_as_expected():
    tt = TranslationTable(rna=True)

    assert tt.rna is True

    assert tt.codons_to_aa["AUG"] == "M"
    assert tt.codons_to_aa["CCC"] == "P"
    assert tt.codons_to_aa["GAC"] == "D"

    assert "ATG" not in tt.codons_to_aa

    assert tt.aa_to_codons["M"] == ("AUG",)
    assert set(tt.aa_to_codons["P"]) == {"CCU", "CCC", "CCA", "CCG"}
    assert set(tt.aa_to_codons["D"]) == {"GAU", "GAC"}


def test_all_forward_codons_appear_in_reverse_table():
    tt = TranslationTable()
    reverse_codons = {
        codon
        for codons in tt.aa_to_codons.values()
        for codon in codons
    }
    assert reverse_codons == set(tt.codons_to_aa)


def test_all_forward_rna_codons_appear_in_reverse_table():
    tt = TranslationTable(rna=True)
    reverse_codons = {
        codon
        for codons in tt.aa_to_codons.values()
        for codon in codons
    }
    assert reverse_codons == set(tt.codons_to_aa)


def test_dna_table_contains_no_rna_codons():
    tt = TranslationTable()
    assert all("U" not in codon for codon in tt.codons_to_aa)
    assert all(
        "U" not in codon
        for codons in tt.aa_to_codons.values()
        for codon in codons
    )


def test_rna_table_contains_no_dna_thymine_codons():
    tt = TranslationTable(rna=True)
    assert all("T" not in codon for codon in tt.codons_to_aa)
    assert all(
        "T" not in codon
        for codons in tt.aa_to_codons.values()
        for codon in codons
    )


def test_codon_tables_are_immutable():
    tt = TranslationTable()

    with pytest.raises(TypeError):
        tt.codons_to_aa["AAA"] = "X"

    with pytest.raises(TypeError):
        tt.aa_to_codons["X"] = ("AAA",)


def test_aa_to_codons_values_are_immutable():
    tt = TranslationTable()
    assert isinstance(tt.aa_to_codons["M"], tuple)
    with pytest.raises(AttributeError):
        tt.aa_to_codons["M"].append("XXX")


def test_rna_attribute_is_read_only():
    tt = TranslationTable(rna=True)
    with pytest.raises(AttributeError):
        tt.rna = False
