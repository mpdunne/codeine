from typing import Any, Dict, Optional, Tuple, Union

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
        rna: Optional[bool] = None,
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

            rna
                Whether to use rna (default is False, i.e. use DNA).
            """
        self._locked = False

        if table is None:
            table = TranslationTable(rna=False if rna is None else rna)

        elif rna is not None and table.rna != rna:
            raise ValueError('table and rna specify inconsistent molecule types.')

        expected_aa = set(table.aa_to_codons)
        observed_aa = {aa.upper() for aa in weights}

        missing_aa = expected_aa - observed_aa
        if missing_aa:
            raise ValueError(f'Missing weights for amino acid(s): {sorted(missing_aa)}')

        extra_aa = observed_aa - expected_aa
        if extra_aa:
            raise ValueError(f'Unknown amino acid(s): {sorted(extra_aa)}')

        weights_flat: Dict[str, float] = {}
        aa_to_codons: Dict[str, Tuple[str, ...]] = {}

        for aa, codon_weights in weights.items():
            aa = aa.upper()

            normalised_codon_weights = {
                table.normalise_sequence(codon): weight
                for codon, weight in codon_weights.items()
            }

            expected_codons = set(table.aa_to_codons[aa])
            observed_codons = set(normalised_codon_weights.keys())

            missing_codons = expected_codons - observed_codons
            if missing_codons:
                raise ValueError(f'Missing weights for codon(s) for amino acid {aa}: {sorted(missing_codons)}')

            extra_codons = observed_codons - expected_codons
            if extra_codons:
                raise ValueError(f'Unknown codon(s) for amino acid {aa}: {sorted(extra_codons)}')

            total = sum(normalised_codon_weights.values())
            if total <= 0:
                raise ValueError(f'Weights for amino acid {aa} must sum to > 0')

            codons = []
            for codon, weight in normalised_codon_weights.items():
                if weight < 0:
                    raise ValueError(f'Weight for codon {codon} cannot be negative')

                codons.append(codon)
                weights_flat[codon] = float(weight) / total

            aa_to_codons[aa] = tuple(codons)

        self.rna = table.rna
        self.aa_to_codons = FrozenDict(aa_to_codons)
        self.weights = FrozenDict(weights_flat)

        self._locked = True

    def __setattr__(self, name: str, value: Any) -> None:
        if getattr(self, '_locked', False):
            raise AttributeError(f'{type(self).__name__} is immutable')
        super().__setattr__(name, value)

    def __repr__(self) -> str:
        molecule = 'RNA' if self.rna else 'DNA'

        lines = [
            f'CodonWeights',
            f'Molecule type: {molecule}',
            '',
            'Weights:',
        ]

        for aa in sorted(self.aa_to_codons):
            weights = ' '.join(
                f'{codon}={self.weights[codon]:.3f}'
                for codon in self.aa_to_codons[aa]
            )
            lines.append(f'    {aa}: {weights}')

        return '\n'.join(lines)

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
        codons = self.aa_to_codons[aa.upper()]
        weights = {codon: self.weights[codon] for codon in codons}
        return weights

    @classmethod
    def uniform(cls, table: Optional[TranslationTable] = None, rna: Optional[bool] = None) -> 'CodonWeights':
        """
        Construct a CodonWeights object with uniform codon weights for a given translation table.

        Parameters
        ----------
        table
            The reference table. If blank, use the standard genetic code.

        rna
            Whether to use RNA.

        Returns
        -------
        A uniform CodonWeights object.
        """
        if table is None:
            table = TranslationTable(rna=False if rna is None else rna)

        elif rna is not None and table.rna != rna:
            raise ValueError('table and rna specify inconsistent molecule types.')

        uniform_weights: WeightDict = {
            aa: {codon: 1.0 for codon in codons}
            for aa, codons in table.aa_to_codons.items()
        }

        return cls(uniform_weights, table=table)

    @classmethod
    def ecoli(cls, rna: Optional[bool] = None) -> 'CodonWeights':
        """
        Construct a CodonWeights object with codon probabilities corresponding to E. coli.

        Weights are obtained from GenScript https://www.genscript.com/tools/codon-frequency-table

        Parameters
        ----------
        rna
            Whether to use RNA.

        Returns
        -------
        A CodonWeights object corresponding to E. Coli
        """
        from codeine.translation.data import ECOLI_WEIGHTS
        return cls(ECOLI_WEIGHTS, rna=rna)

    @classmethod
    def yeast(cls, rna: Optional[bool] = None) -> 'CodonWeights':
        """
        Construct a CodonWeights object with codon probabilities corresponding to 'yeast'.

        Weights are obtained from GenScript https://www.genscript.com/tools/codon-frequency-table

        Parameters
        ----------
        rna
            Whether to use RNA.

        Returns
        -------
        A CodonWeights object corresponding to S. cerevisiea
        """
        from codeine.translation.data import YEAST_WEIGHTS
        return cls(YEAST_WEIGHTS, rna=rna)

    @classmethod
    def human(cls, rna: Optional[bool] = None) -> 'CodonWeights':
        """
        Construct a CodonWeights object with codon probabilities corresponding to Human.

        Weights are obtained from GenScript https://www.genscript.com/tools/codon-frequency-table

        Parameters
        ----------
        rna
            Whether to use RNA.

        Returns
        -------
        A CodonWeights object corresponding to Human.
        """
        from codeine.translation.data import HUMAN_WEIGHTS
        return cls(HUMAN_WEIGHTS, rna=rna)

    @classmethod
    def mouse(cls, rna: Optional[bool] = None) -> 'CodonWeights':
        """
        Construct a CodonWeights object with codon probabilities corresponding to Mouse.

        Weights are obtained from GenScript https://www.genscript.com/tools/codon-frequency-table

        Parameters
        ----------
        rna
            Whether to use RNA.

        Returns
        -------
        A CodonWeights object corresponding to Mouse.
        """
        from codeine.translation.data import MOUSE_WEIGHTS
        return cls(MOUSE_WEIGHTS, rna=rna)

    @classmethod
    def arabidopsis(cls, rna: Optional[bool] = None) -> 'CodonWeights':
        """
        Construct a CodonWeights object with codon probabilities corresponding to A. thaliana.

        Weights are obtained from GenScript https://www.genscript.com/tools/codon-frequency-table

        Parameters
        ----------
        rna
            Whether to use RNA.

        Returns
        -------
        A CodonWeights object corresponding to Arabidopsis thaliana.
        """
        from codeine.translation.data import ARABIDOPSIS_WEIGHTS
        return cls(ARABIDOPSIS_WEIGHTS, rna=rna)

    @classmethod
    def drosophila(cls, rna: Optional[bool] = None) -> 'CodonWeights':
        """
        Construct a CodonWeights object with codon probabilities corresponding to D. melanogaster.

        Weights are obtained from GenScript https://www.genscript.com/tools/codon-frequency-table

        Parameters
        ----------
        rna
            Whether to use RNA.

        Returns
        -------
        A CodonWeights object corresponding to Drosophila melanogaster.
        """
        from codeine.translation.data import DROSOPHILA_WEIGHTS
        return cls(DROSOPHILA_WEIGHTS, rna=rna)
