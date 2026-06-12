from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from codeine.sequence.space import ForbiddenMotif

from typing import Dict, List, Sequence, Union

from codeine.motifs.restriction import RestrictionSite


CodonRestriction = Union[str, Sequence[str]]


def format_count(n: int) -> str:
    """
    Format large counts in scientific notation.
    """
    if n < 10**9:
        return f'{n:,}'

    exponent = len(str(n)) - 1
    mantissa = n / (10 ** exponent)
    return f'{mantissa:.3g} × 10^{exponent}'


def format_restrictions(
        restrictions: Dict[int, CodonRestriction],
        label: str = 'positions',
        max_lines: int = 8,
) -> List[str]:
    """
    Format codon restrictions or pins for repr output.
    """
    if not restrictions:
        return ['    None']

    items = sorted(restrictions.items())

    if len(items) <= max_lines:
        lines = []

        for pos, codons in items:
            if isinstance(codons, str):
                codons = [codons]
            lines.append(f'    {pos}: {" ".join(codons)}')

        return lines

    lines = [f'    {len(items)} {label}']

    for pos, codons in items[:max_lines]:
        if isinstance(codons, str):
            codons = [codons]
        lines.append(f'        {pos}: {" ".join(codons)}')

    lines.append(f'        ... {len(items) - max_lines} more')
    return lines


def format_banned_sequences(
        sequences: Sequence[str],
        max_lines: int = 8,
) -> List[str]:
    """
    Format banned sequences for repr output.
    """
    if not sequences:
        return ['    None']

    if len(sequences) <= max_lines:
        return [f'    {sequence}' for sequence in sequences]

    lines = [f'    {len(sequences)} sequences']
    for sequence in sequences[:max_lines]:
        lines.append(f'        {sequence}')

    lines.append(f'        ... {len(sequences) - max_lines} more')
    return lines


def format_forbidden_motif(motif: 'ForbiddenMotif', rna: bool) -> str:

    if isinstance(motif, RestrictionSite):
        sequences = [
            seq.upper().replace('T', 'U') if rna else seq.upper().replace('U', 'T')
            for seq in motif.motifs
        ]
        return f'{motif.name} ({", ".join(sequences)})'

    return motif.upper().replace('T', 'U') if rna else motif.upper().replace('U', 'T')
