from typing import Dict, List

from Bio.Data import CodonTable as TranslationTable


class CodonTable:
    """
    Basic CodonTable class, using the standard translation table with DNA and uniform probabilites.
    """

    def __init__(self) -> None:
        """
        Constructor for the codon table class.
        """
        dna_to_aa = TranslationTable.unambiguous_dna_by_name["Standard"].forward_table
        aa_to_dna = {}
        for codon, aa in dna_to_aa.items():
            aa_to_dna[aa] = aa_to_dna.get(aa, []) + [codon]

        codon_probabilities = {aa: {codon: 1 / len(codons) for codon in codons} for aa, codons in aa_to_dna.items()}

        self._dna_to_aa = dna_to_aa
        self._aa_to_dna = aa_to_dna
        self._codon_probabilities = codon_probabilities

    @property
    def codons_to_aa(self) -> Dict[str, str]:
        return self._dna_to_aa

    @property
    def aa_to_codons(self) -> Dict[str, List[str]]:
        return self._aa_to_dna

    @property
    def codon_probabilities(self) -> Dict[str, Dict[str, float]]:
        return self._codon_probabilities
