from typing import Any

from Bio.Data import CodonTable

from codeine.utils.dict import FrozenDict


class TranslationTable:
    """
    Translation table. Organises information from BioPython's CodonTable. Uses DNA by
    default. If rna=True, uses RNA, i.e. codons are represented with U instead of T.
    """

    __slots__ = (
        "_locked",
        "table_id",
        "rna",
        "codons_to_aa",
        "aa_to_codons",
    )

    def __init__(self, table_id: int = 1, rna: bool = False) -> None:
        """
        Constructor for the TranslationTable class.

        Parameters
        ----------
        table_id
            Which translation table to use. Default is 1. See https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
        rna
            Whether to use RNA. Default is no (False), i.e. use DNA.
        """
        object.__setattr__(self, "_locked", False)

        try:
            biopython_table = CodonTable.unambiguous_dna_by_id[table_id]
        except KeyError:
            raise ValueError(f"Unknown NCBI translation table: {table_id}")

        dna_to_aa = biopython_table.forward_table
        dna_to_aa = {**dna_to_aa, **{stop: '*' for stop in biopython_table.stop_codons}}

        aa_to_dna = {}
        for codon, aa in dna_to_aa.items():
            aa_to_dna.setdefault(aa, []).append(codon)

        codons_to_aa = {
            self.normalise_codon(codon, rna=rna): aa
            for codon, aa in dna_to_aa.items()
        }

        aa_to_codons = {
            aa: tuple(
                self.normalise_codon(codon, rna=rna)
                for codon in codons
            )
            for aa, codons in aa_to_dna.items()
        }

        object.__setattr__(self, "table_id", table_id)
        object.__setattr__(self, "rna", rna)
        object.__setattr__(self, "codons_to_aa", FrozenDict(codons_to_aa))
        object.__setattr__(self, "aa_to_codons", FrozenDict(aa_to_codons))

        object.__setattr__(self, "_locked", True)

    def __setattr__(self, name: str, value: Any) -> None:
        if getattr(self, "_locked", False):
            raise AttributeError(f"{type(self).__name__} is immutable")
        object.__setattr__(self, name, value)

    def __repr__(self) -> str:
        return f"{type(self).__name__}(table_id={self.table_id}, rna={self.rna})"

    @staticmethod
    def normalise_codon(codon: str, rna: bool = False) -> str:
        codon = codon.upper()
        return codon.replace("T", "U") if rna else codon.replace("U", "T")
