from codeine.motifs.restriction import RestrictionSite, reverse_complement


def test_restriction_site_forward():
    assert RestrictionSite.EcoRI.forward == 'GAATTC'
    assert RestrictionSite.BsaI.forward == 'GGTCTC'


def test_palindromic_site_has_one_motif():
    assert RestrictionSite.EcoRI.motifs == ('GAATTC',)


def test_non_palindromic_site_has_forward_and_reverse():
    assert RestrictionSite.BsaI.motifs == ('GGTCTC', 'GAGACC')


def test_reverse_complement():
    assert reverse_complement('GGTCTC') == 'GAGACC'
