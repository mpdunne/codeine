from itertools import product

from codeine.constraints.count import CountConstraint


CODONS = tuple("".join(codon) for codon in product("ACGT", repeat=3))


class GCConstraint(CountConstraint):
    """
    Constrain the total number or proportion of G/C nucleotides.
    """
    _positions_per_choice = 3
    _choice_counts = {
        codon: sum(nt in 'GC' for nt in codon)
        for codon in CODONS
    }


class GC3Constraint(CountConstraint):
    """
    Constrain the number or proportion of codons with G/C as their last nucleotide.
    """
    _positions_per_choice = 1
    _choice_counts = {
        codon: int(codon[2] in "GC")
        for codon in CODONS
    }
