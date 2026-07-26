
from typing import List, Optional, Sequence, Union

from codeine.motifs.restriction import RestrictionSite

ForbiddenMotif = Union[str, RestrictionSite]
ForbiddenMotifs = Union[ForbiddenMotif, Sequence[ForbiddenMotif]]


def expand_and_validate_forbidden_motifs(
        forbidden_motifs: ForbiddenMotifs,
        rna: bool
) -> List[str]:
    """
    Convert a set of forbidden motifs into a set of sequences to ban.

    Parameters
    ----------
    forbidden_motifs
        A sequence of either dna/rna strings or RestrictionSite objects.
    rna
        Whether to use RNA.

    Returns
    -------
    A list of forbidden nucleotide sequences.
    """
    all_sequences = []

    if forbidden_motifs is None:
        return []

    if isinstance(forbidden_motifs, (str, RestrictionSite)):
        forbidden_motifs = [forbidden_motifs]

    for motif in forbidden_motifs:
        if isinstance(motif, RestrictionSite):
            sequences = [*motif.motifs]

        elif isinstance(motif, str):
            if len(motif) == 0:
                raise ValueError('Forbidden motifs cannot be empty.')

            sequences = [motif]

        else:
            raise TypeError('Forbidden motifs must be strings or codeine.RestrictionSite.')

        sequences = [seq.upper() for seq in sequences]
        sequences = [seq.replace('T', 'U') if rna else seq.replace('U', 'T') for seq in sequences]

        allowed = set('ACGU' if rna else 'ACGT')
        for seq in sequences:
            if not set(seq) <= allowed:
                raise ValueError('Forbidden motifs must be nucleotide sequences.')

        all_sequences += sequences

    return sorted(set(all_sequences))
