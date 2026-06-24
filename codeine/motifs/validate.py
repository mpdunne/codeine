
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


def expand_and_validate_max_homopolymer(
        max_homopolymer: int,
        rna: bool = False
) -> List[str]:
    """
    Convert a max homopolymer constraint into a set of banned sequences.

    Parameters
    ----------
    max_homopolymer
        The max length of homopolymer
    rna
        Whether to use RNA

    Returns
    -------
    A list of forbidden nucleotide sequences.
    """

    if not isinstance(max_homopolymer, int):
        raise TypeError('max_homopolymer must be an integer.')

    if max_homopolymer < 1:
        raise ValueError('max_homopolymer must be at least 1.')

    nts = 'ACGU' if rna else 'ACGT'
    return [nt * (max_homopolymer + 1) for nt in nts]


def expand_and_validate_sequence_constraints(
        forbidden_motifs: Optional[ForbiddenMotifs] = None,
        max_homopolymer: Optional[int] = None,
        rna: bool = False,
):
    """
    Convert forbidden sequences and/or max homopolymer constraints into sets of banned sequences.

    Parameters
    ----------
    forbidden_motifs
        A sequence of either dna/rna strings or RestrictionSite objects.
    max_homopolymer
        The max allowed homopolymer.
    rna
        Whether to use RNA.

    Returns
    -------
    A list of forbidden nucleotide sequences.
    """
    forbidden_sequences = []

    if forbidden_motifs is not None:
        forbidden_sequences += expand_and_validate_forbidden_motifs(forbidden_motifs, rna=rna)

    if max_homopolymer is not None:
        forbidden_sequences += expand_and_validate_max_homopolymer(max_homopolymer, rna=rna)

    return sorted(set(forbidden_sequences))
