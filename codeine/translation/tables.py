from typing import Any

from Bio.Data import CodonTable

from codeine.utils.dict import FrozenDict


class TranslationTable:
    """
    Translation table. Organises information from BioPython's CodonTable. Here we use the
    NCBI translation table IDs (see https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi)
    Uses DNA by default. If rna=True, uses RNA, i.e. codons are represented with U instead of T.
    """

    def __init__(self, table_id: int = 1, rna: bool = False) -> None:
        """
        Constructor for the TranslationTable class. This class is immutable after construction!

        Parameters
        ----------
        table_id
            Which translation table to use. Default is 1.
        rna
            Whether to use RNA. Default is no (False), i.e. use DNA.
        """
        self._locked = False

        self.table_id = table_id
        self.rna = rna

        try:
            biopython_table = CodonTable.unambiguous_dna_by_id[table_id]
        except KeyError:
            raise ValueError(f'Unknown NCBI translation table ID: {table_id}.')

        dna_to_aa = biopython_table.forward_table
        dna_to_aa = {**dna_to_aa, **{stop: '*' for stop in biopython_table.stop_codons}}

        aa_to_dna = {}
        for codon, aa in dna_to_aa.items():
            aa_to_dna.setdefault(aa, []).append(codon)

        codons_to_aa = {
            self.normalise_codon(codon): aa
            for codon, aa in dna_to_aa.items()
        }

        aa_to_codons = {
            aa: tuple(
                self.normalise_codon(codon)
                for codon in codons
            )
            for aa, codons in aa_to_dna.items()
        }

        self.codons_to_aa = FrozenDict(codons_to_aa)
        self.aa_to_codons = FrozenDict(aa_to_codons)

        self._locked = True

    def __setattr__(self, name: str, value: Any) -> None:
        if getattr(self, "_locked", False) and name != "_locked":
            raise AttributeError(f"{type(self).__name__} is immutable")
        object.__setattr__(self, name, value)

    def __repr__(self) -> str:
        return f"{type(self).__name__}(table_id={self.table_id}, rna={self.rna})"

    def normalise_codon(self, codon: str) -> str:
        """
        Format a codon in the format specified by this codon table, i.e. convert to
        RNA/DNA and cast to upper case.

        Parameters
        ----------
        codon
            The inputted codon.

        Returns
        -------
        The normalised codon.
        """
        codon = codon.upper()
        return codon.replace("T", "U") if self.rna else codon.replace("U", "T")
