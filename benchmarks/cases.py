"""
Benchmark case definitions for Codeine.
"""

from copy import deepcopy
from dataclasses import dataclass
from typing import List, Optional

from codeine import CodingSpace, CodonWeights
from codeine.constraints import (
    DirectRepeats,
    ForbiddenMotifs,
    Hairpins,
    InvertedRepeats,
    MaxHomopolymer,
    TandemRepeats,
)

from benchmarks.proteins import DIFFICULT_PROTEINS, LARGE_PROTEINS, NORMAL_PROTEINS


CONSTRAINT_SETS = {
    'none': [],

    # Individual constraints at useful levels of stringency.
    'motifs-small': [ForbiddenMotifs(['GAATTC', 'GGATCC', 'GGTCTC'])],
    'motifs-large': [ForbiddenMotifs(['GAATTC', 'GGATCC', 'GGTCTC', 'GAGACC', 'CTCGAG', 'AAGCTT'])],
    'homopolymer-6': [MaxHomopolymer(6)],
    'homopolymer-4': [MaxHomopolymer(4)],
    'tandem-6x2': [TandemRepeats(6, 2)],
    'tandem-3x3': [TandemRepeats(3, 3)],
    'direct-repeat-18': [DirectRepeats(18, max_distance=120)],
    'direct-repeat-12': [DirectRepeats(12, max_distance=120)],
    'inverted-repeat-15': [InvertedRepeats(15, 3, 90)],
    'inverted-repeat-10': [InvertedRepeats(10, 3, 90)],
    'hairpin-15': [Hairpins(15, 3, 90)],
    'hairpin-10': [Hairpins(10, 3, 90)],

    # Representative combinations used in ordinary biotech workflows.
    'cloning': [ForbiddenMotifs(['GAATTC', 'GGATCC', 'GGTCTC']), MaxHomopolymer(6)],
    'gene-synthesis': [
        ForbiddenMotifs(['GAATTC', 'GGATCC', 'GGTCTC', 'GAGACC', 'CTCGAG', 'AAGCTT']),
        MaxHomopolymer(5),
        TandemRepeats(6, 2),
    ],
    'repeat-sensitive': [
        MaxHomopolymer(4),
        TandemRepeats(3, 3),
        DirectRepeats(18, max_distance=120),
        Hairpins(12, 3, 90),
    ],
}


@dataclass(frozen=True)
class Case:
    """
    A reproducible benchmark workload.

    Parameters
    ----------
    name
        Stable name used in benchmark output and result files.
    group
        Broad workload category used to organise results.
    sequence
        Protein sequence used to construct the coding space.
    constraint_set
        Name of a constraint set in ``CONSTRAINT_SETS``.
    codon_weights
        Optional codon weights passed to ``CodingSpace``.
    quick
        Whether to include the case in the quick benchmark suite.
    mutation_fraction
        Fraction of positions left free in a mutation-space benchmark.
    max_codons
        Optional maximum codon distance for a mutation-space benchmark.
    """

    name: str
    group: str
    sequence: str
    constraint_set: str = 'none'
    codon_weights: Optional[CodonWeights] = None
    quick: bool = False
    mutation_fraction: Optional[float] = None
    max_codons: Optional[int] = None

    def build(self):
        """
        Construct a fresh space for this benchmark case.
        """
        # Constraints are stateful once linked to a graph, so each space gets its own copy.
        constraints = deepcopy(CONSTRAINT_SETS[self.constraint_set])
        space = CodingSpace(self.sequence, constraints=constraints, codon_weights=self.codon_weights, seed=1)

        if self.mutation_fraction is None:
            return space

        coding_sequence = space.sample()
        free_count = max(1, int(len(self.sequence) * self.mutation_fraction))
        return space.mutants(coding_sequence, free_positions=range(1, free_count + 1), max_codons=self.max_codons)


def synthetic_sequence(length: int, unit: str = 'ACDEFGHIKLMNPQRSTVWY') -> str:
    """
    Build a deterministic synthetic protein sequence.
    """
    repeats = (length + len(unit) - 1) // len(unit)
    return (unit * repeats)[:length]


def synthetic_cases() -> List[Case]:
    """
    Return synthetic scaling benchmarks.
    """
    cases = []

    for length in (10, 25, 50, 100, 250, 500, 1000):
        cases.append(Case(f'synthetic/{length:04d}/plain', 'scaling', synthetic_sequence(length), quick=length <= 250))

    for length in (50, 100, 250, 500):
        cases.append(Case(f'synthetic/{length:04d}/gene-synthesis', 'scaling', synthetic_sequence(length), 'gene-synthesis', quick=length <= 100))

    return cases


def real_protein_cases() -> List[Case]:
    """
    Return benchmarks based on representative real proteins.
    """
    cases = []

    for name, sequence in NORMAL_PROTEINS.items():
        quick = name in ('ubiquitin', 'gfp')
        cases.append(Case(f'real/{name}/plain', 'real', sequence, quick=quick))
        cases.append(Case(f'real/{name}/weighted', 'weights', sequence, codon_weights=CodonWeights.ecoli(), quick=quick))

        for constraint_set in ('motifs-small', 'homopolymer-6', 'tandem-6x2', 'direct-repeat-18', 'cloning'):
            cases.append(Case(f'real/{name}/{constraint_set}', 'constraints', sequence, constraint_set, quick=quick))

        for constraint_set in (
            'motifs-large', 'homopolymer-4', 'tandem-3x3', 'direct-repeat-12',
            'inverted-repeat-15', 'inverted-repeat-10', 'hairpin-15', 'hairpin-10',
            'gene-synthesis', 'repeat-sensitive',
        ):
            cases.append(Case(f'real/{name}/{constraint_set}', 'constraints', sequence, constraint_set))

    return cases


def pathological_cases() -> List[Case]:
    """
    Return deliberately repetitive and difficult benchmarks.
    """
    cases = []

    for name, sequence in DIFFICULT_PROTEINS.items():
        cases.append(Case(f'pathological/{name}/plain', 'pathological', sequence))
        cases.append(Case(f'pathological/{name}/tandem-3x3', 'pathological', sequence, 'tandem-3x3'))
        cases.append(Case(f'pathological/{name}/direct-repeat-18', 'pathological', sequence, 'direct-repeat-18'))
        cases.append(Case(f'pathological/{name}/direct-repeat-12', 'pathological', sequence, 'direct-repeat-12'))
        cases.append(Case(f'pathological/{name}/repeat-sensitive', 'pathological', sequence, 'repeat-sensitive'))

    return cases


def large_cases() -> List[Case]:
    """
    Return benchmarks for large protein sequences.
    """
    cases = []

    for name, sequence in LARGE_PROTEINS.items():
        cases.append(Case(f'large/{name}/plain', 'large', sequence))
        cases.append(Case(f'large/{name}/gene-synthesis', 'large', sequence, 'gene-synthesis'))

    return cases


def mutation_cases() -> List[Case]:
    """
    Return mutation-space benchmarks.
    """
    cases = []

    for fraction in (0.05, 0.20, 0.50):
        percent = int(fraction * 100)
        quick = fraction <= 0.20
        sequence = NORMAL_PROTEINS['gfp']
        max_codons = max(1, int(len(sequence) * fraction / 2))
        cases.append(Case(f'mutation/gfp/free-{percent:02d}pct', 'mutation', sequence, quick=quick, mutation_fraction=fraction))
        cases.append(Case(f'mutation/gfp/free-{percent:02d}pct-distance', 'mutation', sequence, mutation_fraction=fraction, max_codons=max_codons))

    return cases


def all_cases() -> List[Case]:
    """
    Return every benchmark case in stable display order.
    """
    return synthetic_cases() + real_protein_cases() + pathological_cases() + large_cases() + mutation_cases()
