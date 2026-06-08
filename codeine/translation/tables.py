from collections.abc import Mapping
from typing import Any, Iterator

from Bio.Data import CodonTable as TranslationTable


class FrozenDict(Mapping):
    """
    Immutable mapping with a normal dict-like repr.
    """

    def __init__(self, data: Mapping) -> None:
        self._data = dict(data)

    def __getitem__(self, key: Any) -> Any:
        return self._data[key]

    def __iter__(self) -> Iterator:
        return iter(self._data)

    def __len__(self) -> int:
        return len(self._data)

    def __repr__(self) -> str:
        return repr(self._data)


class CodonTable:
    """
    Basic codon table using the standard translation table.

    Uses DNA by default. If rna=True, codons are represented with U instead of T.
    Codon probabilities are uniform within each amino acid.
    """

    __slots__ = (
        "_locked",
        "rna",
        "codons_to_aa",
        "aa_to_codons",
        "codon_probabilities",
    )

    def __init__(self, rna: bool = False) -> None:
        """
        Constructor for the CodonTable class.

        Parameters
        ----------
        rna
            Whether to use RNA. Default is no (False), i.e. use DNA.
        """
        object.__setattr__(self, "_locked", False)

        dna_to_aa = TranslationTable.unambiguous_dna_by_name["Standard"].forward_table

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

        codon_probabilities = {
            aa: {codon: 1 / len(codons) for codon in codons}
            for aa, codons in aa_to_codons.items()
        }

        object.__setattr__(self, "rna", rna)
        object.__setattr__(self, "codons_to_aa", FrozenDict(codons_to_aa))
        object.__setattr__(self, "aa_to_codons", FrozenDict(aa_to_codons))
        object.__setattr__(self, "codon_probabilities",
                           FrozenDict({
                               aa: FrozenDict(probs)
                               for aa, probs in codon_probabilities.items()
                           }),
        )

        object.__setattr__(self, "_locked", True)

    def __setattr__(self, name: str, value: Any) -> None:
        if getattr(self, "_locked", False):
            raise AttributeError(f"{type(self).__name__} is immutable")
        object.__setattr__(self, name, value)

    @staticmethod
    def normalise_codon(codon: str, rna: bool = False) -> str:
        codon = codon.upper()
        return codon.replace("T", "U") if rna else codon.replace("U", "T")
