from typing import Any, Dict, Optional, Tuple, Union

from codeine.translation.tables import TranslationTable
from codeine.utils.dict import FrozenDict


WeightDict = Dict[str, Dict[str, Union[float, int]]]


class CodonWeights:
    """
    A class to store codon weights, for example codon usage information for a
    specific organism.

    Input weights are grouped by amino acid:

    .. code-block:: python

        {
            'A': {'GCT': 1.0, 'GCC': 1.0, ...},
            'R': {'CGT': 1.0, 'CGC': 1.0, ...},
            ...
        }

    Stored weights are flat:

    .. code-block:: python

        {
            'GCT': 0.25,
            'GCC': 0.25,
            ...
        }

    Weights are normalised per amino acid. Weights can be inputted as either
    RNA or DNA, they will be normalised later.
    """

    def __init__(
            self,
            weights: WeightDict,
    ) -> None:
        """
            Parameters
            ----------
            weights
                Codon weights grouped by amino acid.

                Example:

                .. code-block:: python

                    {
                        'A': {'GCT': 1.0, 'GCC': 1.0, 'GCA': 1.0, 'GCG': 1.0},
                        'R': {'CGT': 1.0, 'CGC': 1.0, 'CGA': 1.0,
                              'CGG': 1.0, 'AGA': 1.0, 'AGG': 1.0},
                        ...
                    }

            """
        self._locked = False

        weights_flat: Dict[str, float] = {}
        aa_to_codons: Dict[str, Tuple[str, ...]] = {}
        observed_codons = set()

        for aa, codon_weights in weights.items():
            aa = aa.upper()

            if aa in aa_to_codons:
                raise ValueError(f'Duplicate weights for amino acid {aa}')

            normalised_codon_weights: Dict[str, Union[float, int]] = {}

            for codon, weight in codon_weights.items():
                codon = codon.upper()

                if codon in normalised_codon_weights:
                    raise ValueError(f'Duplicate weight for codon {codon} for amino acid {aa}')

                if codon in observed_codons:
                    raise ValueError(f'Duplicate weight for codon {codon}')

                if weight < 0:
                    raise ValueError(f'Weight for codon {codon} cannot be negative')

                normalised_codon_weights[codon] = weight
                observed_codons.add(codon)

            total = sum(normalised_codon_weights.values())
            if total <= 0:
                raise ValueError(f'Weights for amino acid {aa} must sum to > 0')

            codons = tuple(normalised_codon_weights.keys())

            for codon, weight in normalised_codon_weights.items():
                weights_flat[codon] = float(weight) / total

            aa_to_codons[aa] = codons

        self.aa_to_codons = FrozenDict(aa_to_codons)
        self.weights = FrozenDict(weights_flat)

        self._locked = True

    def __setattr__(self, name: str, value: Any) -> None:
        if getattr(self, '_locked', False):
            raise AttributeError(f'{type(self).__name__} is immutable')
        super().__setattr__(name, value)

    def __repr__(self) -> str:
        lines = [
            'CodonWeights',
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

    def for_table(self, table: TranslationTable) -> 'CodonWeights':
        """
        Resolve this CodonWeights class against a translation table, i.e. check that the
        codons match and normalise to RNA/DNA, and return the resolved version.

        Parameters
        ----------
        table
            The translation table to resolve against.

        Returns
        -------
        A resolved ``CodonWeights`` object.
        """
        expected_aas = set(table.aa_to_codons.keys())
        observed_aas = set(self.aa_to_codons.keys())

        missing_aas = expected_aas - observed_aas
        if missing_aas:
            raise ValueError(
                f'Missing weights for amino acid(s): {sorted(missing_aas)}'
            )

        extra_aas = observed_aas - expected_aas
        if extra_aas:
            raise ValueError(
                f'Unknown amino acid(s): {sorted(extra_aas)}'
            )

        resolved_weights: WeightDict = {}
        resolved_codons = set()
        changed = False

        for aa, codons in self.aa_to_codons.items():
            codon_weights: Dict[str, float] = {}

            for codon in codons:
                resolved_codon = table.normalise_sequence(codon)

                if resolved_codon != codon:
                    changed = True

                if resolved_codon in resolved_codons:
                    raise ValueError(
                        f'Duplicate weight for codon {resolved_codon}'
                    )

                resolved_codons.add(resolved_codon)
                codon_weights[resolved_codon] = self.weights[codon]

            expected_codons = set(table.aa_to_codons[aa])
            observed_codons = set(codon_weights.keys())

            missing_codons = expected_codons - observed_codons
            if missing_codons:
                raise ValueError(
                    f'Missing weights for codon(s) for amino acid {aa}: '
                    f'{sorted(missing_codons)}'
                )

            extra_codons = observed_codons - expected_codons
            if extra_codons:
                raise ValueError(
                    f'Unknown codon(s) for amino acid {aa}: '
                    f'{sorted(extra_codons)}'
                )

            resolved_weights[aa] = codon_weights

        if not changed:
            return self

        return type(self)(resolved_weights)

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
    def uniform(cls, table: Optional[TranslationTable] = None) -> 'CodonWeights':
        """
        Construct a ``CodonWeights`` object with uniform codon weights for a given translation table.

        Parameters
        ----------
        table
            The reference table. If blank, use the standard genetic code.

        Returns
        -------
        A uniform CodonWeights object.
        """
        if table is None:
            table = TranslationTable()

        uniform_weights: WeightDict = {
            aa: {codon: 1.0 for codon in codons}
            for aa, codons in table.aa_to_codons.items()
        }

        return cls(uniform_weights)

    @classmethod
    def ecoli(cls) -> 'CodonWeights':
        """
        Construct a ``CodonWeights`` object with codon probabilities corresponding to E. coli.

        Weights are based on codon usage counts from GenScript:

        https://www.genscript.com/tools/codon-frequency-table

        Returns
        -------
        A ``CodonWeights`` object corresponding to E. Coli
        """
        from codeine.translation.data.weights import ECOLI_WEIGHTS
        return cls(ECOLI_WEIGHTS)

    @classmethod
    def yeast(cls) -> 'CodonWeights':
        """
        Construct a ``CodonWeights`` object with codon probabilities corresponding to 'yeast'.

        Weights are based on codon usage counts from GenScript:

        https://www.genscript.com/tools/codon-frequency-table

        Returns
        -------
        A ``CodonWeights`` object corresponding to S. cerevisiea
        """
        from codeine.translation.data.weights import YEAST_WEIGHTS
        return cls(YEAST_WEIGHTS)

    @classmethod
    def human(cls) -> 'CodonWeights':
        """
        Construct a ``CodonWeights`` object with codon probabilities corresponding to Human.

        Weights are based on codon usage counts from GenScript:

        https://www.genscript.com/tools/codon-frequency-table

        Returns
        -------
        A ``CodonWeights`` object corresponding to Human.
        """
        from codeine.translation.data.weights import HUMAN_WEIGHTS
        return cls(HUMAN_WEIGHTS)

    @classmethod
    def mouse(cls) -> 'CodonWeights':
        """
        Construct a ``CodonWeights`` object with codon probabilities corresponding to Mouse.

        Weights are based on codon usage counts from GenScript:

        https://www.genscript.com/tools/codon-frequency-table

        Returns
        -------
        A ``CodonWeights`` object corresponding to Mouse.
        """
        from codeine.translation.data.weights import MOUSE_WEIGHTS
        return cls(MOUSE_WEIGHTS)

    @classmethod
    def arabidopsis(cls) -> 'CodonWeights':
        """
        Construct a ``CodonWeights`` object with codon probabilities corresponding to A. thaliana.

        Weights are based on codon usage counts from GenScript:

        https://www.genscript.com/tools/codon-frequency-table

        Returns
        -------
        A ``CodonWeights`` object corresponding to Arabidopsis thaliana.
        """
        from codeine.translation.data.weights import ARABIDOPSIS_WEIGHTS
        return cls(ARABIDOPSIS_WEIGHTS)

    @classmethod
    def drosophila(cls) -> 'CodonWeights':
        """
        Construct a ``CodonWeights`` object with codon probabilities corresponding to D. melanogaster.

        Weights are based on codon usage counts from GenScript:

        https://www.genscript.com/tools/codon-frequency-table

        Returns
        -------
        A ``CodonWeights`` object corresponding to Drosophila melanogaster.
        """
        from codeine.translation.data.weights import DROSOPHILA_WEIGHTS
        return cls(DROSOPHILA_WEIGHTS)
