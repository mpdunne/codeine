import pytest

from codeine.motifs.restriction import RestrictionSite
from codeine.motifs.constraints import expand_and_validate_forbidden_motifs, \
    expand_and_validate_max_homopolymer, expand_and_validate_sequence_constraints


def test_expand_forbidden_motifs_none():
    assert expand_and_validate_forbidden_motifs(None, rna=False) == []


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


def test_orbidden_motif_empty_raises():
    with pytest.raises(ValueError, match='Forbidden motifs cannot be empty'):
        expand_and_validate_forbidden_motifs('', rna=False)


def test_forbidden_motif_invalid_type_raises():
    with pytest.raises(TypeError, match='Forbidden motifs must be strings or codeine.RestrictionSite.'):
        expand_and_validate_forbidden_motifs([420], rna=False)


def test_validate_max_homopolymer_none():
    assert expand_and_validate_max_homopolymer(None) == []


def test_validate_max_homopolymer_int():
    assert expand_and_validate_max_homopolymer(4) == ['AAAAA', 'CCCCC', 'GGGGG', 'TTTTT']


def test_validate_max_homopolymer_rejects_non_int():
    with pytest.raises(TypeError, match='max_homopolymer must be an integer'):
        expand_and_validate_max_homopolymer(4.5)


def test_validate_max_homopolymer_rejects_less_than_one():
    with pytest.raises(ValueError, match='max_homopolymer must be at least 1'):
        expand_and_validate_max_homopolymer(0)


def test_validate_mixed_restrictions():
    max_homopolymer = 4
    motifs = [RestrictionSite.BsaI, 'GGTTCC']
    result = expand_and_validate_sequence_constraints(max_homopolymer=max_homopolymer, forbidden_motifs=motifs)
    assert set(result) == {'GAGACC', 'GGTCTC', 'GGTTCC', 'AAAAA', 'CCCCC', 'GGGGG', 'TTTTT'}
