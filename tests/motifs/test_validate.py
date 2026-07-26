import pytest

from codeine.motifs.restriction import RestrictionSite
from codeine.motifs.validate import expand_and_validate_forbidden_motifs


def test_expand_forbidden_motifs_single_string_dna():
    motifs = 'gaattc'
    validated = expand_and_validate_forbidden_motifs(motifs, rna=False)
    assert validated == ['GAATTC']


def test_expand_forbidden_motifs_single_string_rna():
    motifs = 'GAATTC'
    validated = expand_and_validate_forbidden_motifs(motifs, rna=True)
    assert validated == ['GAAUUC']


def test_expand_forbidden_motifs_sequence_of_strings_deduplicates_and_sorts():
    motifs = ['tttt', 'UUUU', 'aaaa']
    validated = expand_and_validate_forbidden_motifs(motifs, rna=False)
    assert validated == ['AAAA', 'TTTT']


def test_expand_forbidden_motifs_restriction_sites():
    validated = expand_and_validate_forbidden_motifs([RestrictionSite.EcoRI], rna=False)
    assert validated == ['GAATTC']

    validated = expand_and_validate_forbidden_motifs([RestrictionSite.BsaI], rna=False)
    assert validated == ['GAGACC', 'GGTCTC']

    validated = expand_and_validate_forbidden_motifs([RestrictionSite.BsaI, 'GGTTCC'], rna=False)
    assert validated == ['GAGACC', 'GGTCTC', 'GGTTCC']


def test_expand_forbidden_motifs_mixed():
    motifs = [RestrictionSite.EcoRI, 'AAAA']
    validated = expand_and_validate_forbidden_motifs(motifs, rna=False)
    assert validated == ['AAAA', 'GAATTC']


def test_forbidden_motif_empty_raises():
    with pytest.raises(ValueError, match='Forbidden motifs cannot be empty'):
        expand_and_validate_forbidden_motifs('', rna=False)


def test_forbidden_motif_invalid_type_raises():
    with pytest.raises(TypeError, match='Forbidden motifs must be strings or codeine.RestrictionSite.'):
        expand_and_validate_forbidden_motifs([420], rna=False)


def test_forbidden_motif_invalid_nucleotide_raises():
    with pytest.raises(ValueError, match='Forbidden motifs must be nucleotide sequences'):
        expand_and_validate_forbidden_motifs('MANCHEGO', rna=False)

    with pytest.raises(ValueError, match='Forbidden motifs must be nucleotide sequences'):
        expand_and_validate_forbidden_motifs('MANCHEGO', rna=True)


def test_restriction_site_converts_to_rna():
    result = expand_and_validate_forbidden_motifs(RestrictionSite.BsaI, rna=True)
    assert result == ['GAGACC', 'GGUCUC']
