from typing import Any, Dict, Optional, Union

from codeine.translation.tables import TranslationTable
from codeine.utils.dict import FrozenDict


WeightDict = Dict[str, Dict[str, Union[float, int]]]


class CodonWeights:
    """
    A class to store codon weights, for example codon usage information for a
    specific organism.

    Input weights are grouped by amino acid:

        {
            'A': {'GCT': 1.0, 'GCC': 1.0, ...},
            'R': {'CGT': 1.0, 'CGC': 1.0, ...},
            ...
        }

    Stored weights are flat:

        {
            'GCT': 0.25,
            'GCC': 0.25,
            ...
        }

    Weights are normalised per amino acid.
    """

    def __init__(
        self,
        weights: WeightDict,
        table: Optional[TranslationTable] = None,
    ) -> None:
        """
            Constructor for the CodonWeights class.

            Parameters
            ----------
            weights
                Codon weights grouped by amino acid, for codons in the TranslationTable.

                Example:

                    {
                        'A': {'GCT': 1.0, 'GCC': 1.0, 'GCA': 1.0, 'GCG': 1.0},
                        'R': {'CGT': 1.0, 'CGC': 1.0, 'CGA': 1.0,
                              'CGG': 1.0, 'AGA': 1.0, 'AGG': 1.0},
                        ...
                    }

            table
                Translation table defining the codon-AA mapping.
            """
        self._locked = False

        if table is None:
            table = TranslationTable()
            
        flat_weights: Dict[str, float] = {}

        expected_aas = set(table.aa_to_codons)
        actual_aas = {aa.upper() for aa in weights}

        missing_aas = expected_aas - actual_aas
        if missing_aas:
            raise ValueError(f'Missing weights for amino acids: {sorted(missing_aas)}')

        extra_aas = actual_aas - expected_aas
        if extra_aas:
            raise ValueError(f'Unexpected amino acids in weights: {sorted(extra_aas)}')

        for aa, codons in table.aa_to_codons.items():
            codon_weights = weights[aa]

            expected_codons = set(codons)
            actual_codons = set(codon_weights)

            missing_codons = expected_codons - actual_codons
            if missing_codons:
                raise ValueError(f'Missing codon weights for amino acid {aa}: {sorted(missing_codons)}')

            extra_codons = actual_codons - expected_codons
            if extra_codons:
                raise ValueError(f'Unexpected codons for amino acid {aa}: {sorted(extra_codons)}')

            total = sum(codon_weights.values())
            if total <= 0:
                raise ValueError(f'Weights for amino acid {aa} must sum to > 0')

            for codon, weight in codon_weights.items():
                if weight < 0:
                    raise ValueError(f'Weight for codon {codon} cannot be negative')

                flat_weights[codon] = float(weight) / total

        self.table = table
        self.weights = FrozenDict(flat_weights)

        self._locked = True

    def __setattr__(self, name: str, value: Any) -> None:
        if getattr(self, '_locked', False):
            raise AttributeError(f'{type(self).__name__} is immutable')
        super().__setattr__(name, value)

    def __repr__(self) -> str:
        return f'{type(self).__name__}(table={self.table!r})'

    def __getitem__(self, codon: str) -> float:
        return self.weights[codon]

    def by_aa(self, aa: str) -> Dict[str, float]:
        """
        Return the codon weights corresponding to a particular AA.

        Parameters
        ----------
        aa
            The amino acid of interest.

        Returns
        -------
        A set of codon weights keyed by codon.
        """
        codons = self.table.aa_to_codons[aa.upper()]
        weights = {codon: self.weights[codon] for codon in codons}
        return weights

    @classmethod
    def uniform(cls, table: Optional[TranslationTable] = None) -> 'CodonWeights':
        """
        Construct a CodonWeights object with uniform codon weights for a given translation table.

        Parameters
        ----------
        table
            The reference table. If blank, use the standard genetic code.

        Returns
        -------
        A uniform CodonWeights object.
        """
        table = table or TranslationTable()

        uniform_weights: WeightDict = {
            aa: {codon: 1.0 for codon in codons}
            for aa, codons in table.aa_to_codons.items()
        }

        return cls(uniform_weights, table=table)
