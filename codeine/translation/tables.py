from types import MappingProxyType
from typing import Mapping

from Bio.Data import CodonTable as TranslationTable


class CodonTable:
    """
    Basic codon table using the standard translation table.

    Uses DNA by default. If rna=True, codons are represented with U instead of T.
    Codon probabilities are uniform within each amino acid.
    """

    def __init__(self, rna: bool = False) -> None:
        self.rna = rna

        dna_to_aa = dict(
            TranslationTable.unambiguous_dna_by_name["Standard"].forward_table
        )

        aa_to_dna: dict[str, list[str]] = {}
        for codon, aa in dna_to_aa.items():
            aa_to_dna.setdefault(aa, []).append(codon)

        if rna:
            codons_to_aa = {
                self._normalise_codon(codon, rna=True): aa
                for codon, aa in dna_to_aa.items()
            }
            aa_to_codons = {
                aa: tuple(
                    self._normalise_codon(codon, rna=True)
                    for codon in codons
                )
                for aa, codons in aa_to_dna.items()
            }
        else:
            codons_to_aa = {
                self._normalise_codon(codon, rna=False): aa
                for codon, aa in dna_to_aa.items()
            }
            aa_to_codons = {
                aa: tuple(
                    self._normalise_codon(codon, rna=False)
                    for codon in codons
                )
                for aa, codons in aa_to_dna.items()
            }

        codon_probabilities = {
            aa: {
                codon: 1 / len(codons)
                for codon in codons
            }
            for aa, codons in aa_to_codons.items()
        }

        self.codons_to_aa = self._freeze_dict(codons_to_aa)
        self.aa_to_codons = self._freeze_dict(aa_to_codons)
        self.codon_probabilities = self._freeze_nested_dict(codon_probabilities)

    @staticmethod
    def _normalise_codon(codon: str, rna: bool = False) -> str:
        codon = codon.upper()
        return codon.replace("T", "U") if rna else codon.replace("U", "T")

    @staticmethod
    def _freeze_dict(d: Mapping):
        return MappingProxyType(dict(d))

    @staticmethod
    def _freeze_nested_dict(d: Mapping):
        return MappingProxyType({
            key: MappingProxyType(dict(value))
            for key, value in d.items()
        })
