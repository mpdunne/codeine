import textwrap

from collections import Counter
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from codeine.constraints.base import Constraint
    from codeine.constraints.motifs import Motif

from typing import Dict, List, Sequence, Union

from codeine.motifs.restriction import RestrictionSite


CodonRestriction = Union[str, Sequence[str]]


def format_count(n: int) -> str:
    """
    Format large counts in scientific notation.

    Parameters
    ----------
    n
        The count to format.

    Returns
    -------
    str
        The formatted count.
    """
    if n < 10**9:
        return f'{n:,}'

    exponent = len(str(n)) - 1
    mantissa = n / (10 ** exponent)
    return f'{mantissa:.3g} × 10^{exponent}'


import textwrap

def format_sequence(
        sequence: str,
        max_lines: int = 4,
        line_length: int = 80,
) -> List[str]:
    """
    Format and, if necessary, truncate a sequence for repr output.
    """
    max_length = max_lines * line_length

    if len(sequence) > max_length:
        sequence = sequence[:max_length - 3] + '...'

    return textwrap.wrap(sequence, width=line_length)


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


def format_forbidden_motifs(
        motifs: Sequence['Motif'],
        rna: bool,
        max_motifs: int = 8,
) -> str:
    """
    Format forbidden motifs for repr output.
    """
    formatted = sorted({
        format_forbidden_motif(motif, rna)
        for motif in motifs
    })

    if len(formatted) <= max_motifs:
        return ', '.join(formatted)

    return f'{len(formatted)} motifs'


def format_constraints(
        constraints: Sequence['Constraint'],
        rna: bool,
        max_motifs: int = 8,
) -> List[str]:
    """
    Summarise constraints for repr output.
    """
    from codeine.constraints.hairpins import HairpinConstraint
    from codeine.constraints.homopolymers import HomopolymerConstraint
    from codeine.constraints.motifs import ForbiddenMotifConstraint
    from codeine.constraints.repeats import DirectRepeatConstraint, InvertedRepeatConstraint
    from codeine.constraints.tandem import TandemRepeatConstraint

    forbidden_motifs = []
    homopolymer_lengths = []
    constraint_counts = Counter()
    other_count = 0

    for constraint in constraints:
        if isinstance(constraint, HomopolymerConstraint):
            homopolymer_lengths.append(constraint.max_length)

        elif isinstance(constraint, ForbiddenMotifConstraint):
            forbidden_motifs.extend(constraint.motifs)

        elif isinstance(constraint, HairpinConstraint):
            constraint_counts['Hairpins'] += 1

        elif isinstance(constraint, DirectRepeatConstraint):
            constraint_counts['Direct repeats'] += 1

        elif isinstance(constraint, InvertedRepeatConstraint):
            constraint_counts['Inverted repeats'] += 1

        elif isinstance(constraint, TandemRepeatConstraint):
            constraint_counts['Tandem repeats'] += 1

        else:
            other_count += 1

    lines = []

    if forbidden_motifs:
        motifs = format_forbidden_motifs(
            forbidden_motifs,
            rna=rna,
            max_motifs=max_motifs,
        )
        lines.append(f'    Forbidden motifs: {motifs}')

    if homopolymer_lengths:
        lines.append(f'    Maximum homopolymer length: {min(homopolymer_lengths)}')

    for label in ('Direct repeats', 'Inverted repeats', 'Tandem repeats', 'Hairpins'):
        count = constraint_counts[label]

        if count:
            suffix = 'constraint' if count == 1 else 'constraints'
            lines.append(f'    {label}: {count} {suffix}')

    if other_count:
        suffix = 'constraint' if other_count == 1 else 'constraints'
        lines.append(f'    Other: {other_count} {suffix}')

    return lines or ['    None']


def normalise_motif(seq: str, rna: bool) -> str:
    """
    Normalise a motif by uppercasing and converting to RNA/DNA as specified.
    """
    seq = seq.upper()
    return seq.replace('T', 'U') if rna else seq.replace('U', 'T')


def format_forbidden_motif(motif: 'Motif', rna: bool) -> str:
    """
    Format a forbidden motif for display.
    """
    if isinstance(motif, RestrictionSite):
        sequences = [normalise_motif(seq, rna) for seq in motif.motifs]
        return f'{motif.name} ({", ".join(sequences)})'

    return normalise_motif(motif, rna)


def format_positions(positions) -> str:
    """
    Format positions nicely.
    """
    positions = sorted(positions)
    if not positions:
        return 'None'

    ranges = []
    start = prev = positions[0]

    for pos in positions[1:]:
        if pos == prev + 1:
            prev = pos
        else:
            ranges.append((start, prev))
            start = prev = pos

    ranges.append((start, prev))

    return ', '.join(
        str(start) if start == end else f'{start}-{end}'
        for start, end in ranges
    )
