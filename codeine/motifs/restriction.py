"""
Common restriction site motifs, obtained from New England BioLabs here:

https://www.neb.com/en/tools-and-resources/selection-charts/alphabetized-list-of-recognition-specificities

Retrieved: 2026-06-09
"""

from enum import Enum
from typing import Tuple

_COMPLEMENT = str.maketrans('ACGTacgt', 'TGCAtgca')


def reverse_complement(seq: str) -> str:
    """
    Return the reverse complement of a DNA sequence.
    """
    return seq.translate(_COMPLEMENT)[::-1]


class RestrictionSite(Enum):
    """
    Common restriction enzyme recognition sequences.
    """

    # BioBricks
    EcoRI = 'GAATTC'
    XbaI = 'TCTAGA'
    SpeI = 'ACTAGT'
    PstI = 'CTGCAG'

    # Cloning
    BamHI = 'GGATCC'
    HindIII = 'AAGCTT'
    XhoI = 'CTCGAG'
    SalI = 'GTCGAC'
    KpnI = 'GGTACC'
    SacI = 'GAGCTC'
    NcoI = 'CCATGG'
    NdeI = 'CATATG'
    NotI = 'GCGGCCGC'
    MluI = 'ACGCGT'
    AgeI = 'ACCGGT'
    AvrII = 'CCTAGG'
    BglII = 'AGATCT'

    # Golden Gate
    BsaI = 'GGTCTC'
    BsmBI = 'CGTCTC'
    BbsI = 'GAAGAC'
    SapI = 'GCTCTTC'

    def __repr__(self) -> str:
        return f'RestrictionSite.{self.name} ({" / ".join(self.motifs)})'

    def __str__(self) -> str:
        return f'RestrictionSite.{self.name}'

    @property
    def forward(self) -> str:
        """
        Forward recognition sequence.
        """
        return self.value

    @property
    def reverse(self) -> str:
        """
        Reverse-complemented recognition sequence.
        """
        return reverse_complement(self.value)

    @property
    def motifs(self) -> Tuple[str, ...]:
        """
        All motifs corresponding to this restriction site (forward and reverse).
        Palindromic sites return a single motif.
        """
        if self.forward == self.reverse:
            return self.forward,

        return self.forward, self.reverse
