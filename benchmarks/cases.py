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

from benchmarks.proteins import ANTIBODIES, DIFFICULT_PROTEINS, LARGE_PROTEINS, NORMAL_PROTEINS


CONSTRAINT_SETS = {
    'none': [],
    'artificial-repeat-4': [DirectRepeats(4)],
    'artificial-repeat-9': [DirectRepeats(9)],
    'artificial-repeat-15': [DirectRepeats(15)],

    # Individual constraints at useful levels of stringency.
    'motifs-small': [ForbiddenMotifs(['GAATTC', 'GGATCC', 'GGTCTC'])],
    'motifs-large': [ForbiddenMotifs(['GAATTC', 'GGATCC', 'GGTCTC', 'CGTCTC', 'GAAGAC', 'AAGCTT'])],
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

    # Representative gene-synthesis constraint combinations.
    'basic': [
        ForbiddenMotifs(['GGTCTC', 'CGTCTC', 'GAAGAC', 'GAATTC', 'GGATCC', 'AAGCTT']),
        MaxHomopolymer(6),
        TandemRepeats(2, 4), TandemRepeats(3, 3), TandemRepeats(4, 3),
        TandemRepeats(5, 3), TandemRepeats(6, 3),
    ],
    'gene-synthesis': [
        ForbiddenMotifs(['GGTCTC', 'CGTCTC', 'GAAGAC', 'GAATTC', 'GGATCC', 'AAGCTT']),
        MaxHomopolymer(6),
        TandemRepeats(2, 4), TandemRepeats(3, 3), TandemRepeats(4, 3),
        TandemRepeats(5, 3), TandemRepeats(6, 3),
        DirectRepeats(18),
        InvertedRepeats(18),
        Hairpins(12, 3, 8),
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
    mutation_fraction
        Fraction of positions left free in a mutation-space benchmark.
    quick
        Whether to include the case in the quick benchmark suite.
    max_codons
        Optional maximum codon distance for a mutation-space benchmark.
    """

    name: str
    group: str
    sequence: str
    constraint_set: str = 'none'
    codon_weights: Optional[CodonWeights] = None
    mutation_fraction: Optional[float] = None
    quick: bool = True
    max_codons: Optional[int] = None

    def build(self, compiler='flat'):
        """
        Construct a fresh space for this benchmark case.
        """
        # Constraints are stateful once linked to a graph, so each space gets its own copy.
        constraints = deepcopy(CONSTRAINT_SETS[self.constraint_set])
        space = CodingSpace(self.sequence, constraints=constraints, codon_weights=self.codon_weights, seed=1, compiler=compiler)

        if self.mutation_fraction is None:
            return space

        # Seeded samples differ between engines. Use the same ordered reference
        # so mutation benchmarks describe identical spaces for both compilers.
        coding_sequence = space[0]
        free_count = max(1, int(len(self.sequence) * self.mutation_fraction))
        return space.mutants(coding_sequence, free_positions=range(1, free_count + 1), max_codons=self.max_codons)


def all_cases() -> List[Case]:
    """
    Return the standard benchmark suite in stable display order.
    """
    ubiquitin = NORMAL_PROTEINS['ubiquitin']
    gfp = NORMAL_PROTEINS['gfp']
    sfgfp = NORMAL_PROTEINS['sfgfp']
    mcherry = NORMAL_PROTEINS['mcherry']
    luciferase = NORMAL_PROTEINS['luciferase']
    caplacizumab = ANTIBODIES['caplacizumab']
    spcas9 = LARGE_PROTEINS['spcas9']

    return [
        # Unconstrained baselines across proteins of different sizes and character.
        Case('baseline/ubiquitin', 'baseline', ubiquitin),
        Case('baseline/gfp', 'baseline', gfp),
        Case('baseline/mcherry', 'baseline', mcherry),
        Case('baseline/luciferase', 'baseline', luciferase),
        Case('baseline/caplacizumab', 'baseline', caplacizumab),
        Case('baseline/spcas9', 'baseline', spcas9),

        # Individual constraints at two useful levels. Ordinary proteins cover the
        # simpler cases; naturally repetitive proteins exercise repeat constraints.
        Case('constraints/gfp/motifs-small', 'constraints', gfp, 'motifs-small'),
        Case('constraints/gfp/motifs-large', 'constraints', gfp, 'motifs-large'),
        Case('constraints/gfp/homopolymer-6', 'constraints', gfp, 'homopolymer-6'),
        Case('constraints/luciferase/homopolymer-4', 'constraints', luciferase, 'homopolymer-4'),
        Case('constraints/collagen/tandem-6x2', 'constraints', DIFFICULT_PROTEINS['collagen'], 'tandem-6x2'),
        Case('constraints/collagen/tandem-3x3', 'constraints', DIFFICULT_PROTEINS['collagen'], 'tandem-3x3'),
        Case('constraints/elastin/direct-repeat-18', 'constraints', DIFFICULT_PROTEINS['elastin'], 'direct-repeat-18'),
        Case('constraints/elastin/direct-repeat-12', 'constraints', DIFFICULT_PROTEINS['elastin'], 'direct-repeat-12', quick=False),
        Case('constraints/luciferase/inverted-repeat-15', 'constraints', luciferase, 'inverted-repeat-15'),
        Case('constraints/luciferase/inverted-repeat-10', 'constraints', luciferase, 'inverted-repeat-10', quick=False),
        Case('constraints/luciferase/hairpin-15', 'constraints', luciferase, 'hairpin-15'),
        Case('constraints/luciferase/hairpin-10', 'constraints', luciferase, 'hairpin-10', quick=False),

        # Full practical constraint stacks on several proteins. sfGFP and SpCas9
        # reproduce the documentation examples; the others broaden the workload.
        Case('full/sfgfp/gene-synthesis', 'full', sfgfp, 'gene-synthesis', CodonWeights.ecoli()),
        Case('full/luciferase/gene-synthesis', 'full', luciferase, 'gene-synthesis'),
        Case('full/caplacizumab/gene-synthesis', 'full', caplacizumab, 'gene-synthesis'),
        Case('full/spcas9/basic', 'full', spcas9, 'basic', quick=False),

        # Artificial repeats from the original factorisation experiments.
        Case('artificial/mikey-repeat-9', 'artificial', 'MIKEYMIKEY', 'artificial-repeat-9'),
        Case('artificial/mikey-spacer-repeat-9', 'artificial', 'MIKEYAAAAAMIKEY', 'artificial-repeat-9'),
        Case('artificial/mikey-spacer-repeat-15', 'artificial', 'MIKEYAAAAAMIKEY', 'artificial-repeat-15'),
        Case('artificial/rnykqt-repeat-4', 'artificial', 'RNYKQT', 'artificial-repeat-4'),
        Case('artificial/iherqw-repeat-4', 'artificial', 'IHERQW', 'artificial-repeat-4'),
        Case('artificial/mikey-sassafras-repeat-15', 'artificial',
             'MIKEYSASSAFRASMIKEYSASSAFRAS', 'artificial-repeat-15'),
        Case('artificial/mikey-sassafras-repeat-15-weighted', 'artificial',
             'MIKEYSASSAFRASMIKEYSASSAFRAS', 'artificial-repeat-15', CodonWeights.ecoli()),

        # A deliberately awkward repeat-heavy sequence retained as a stress case.
        Case('pathological/hard/direct-repeat-18', 'pathological', DIFFICULT_PROTEINS['hard'], 'direct-repeat-18'),

        # Representative local redesign workloads.
        Case('mutation/sfgfp/free-20pct', 'mutation', sfgfp, mutation_fraction=0.20),
        Case('mutation/sfgfp/gene-synthesis-distance', 'mutation', sfgfp, 'gene-synthesis', mutation_fraction=0.20, max_codons=25),
    ]
