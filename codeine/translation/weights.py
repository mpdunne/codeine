from typing import Any, Dict, Mapping, Optional

from codeine.translation.tables import TranslationTable
from codeine.utils.dict import FrozenDict


class CodonWeights:
    """
    Codon usage weights for a TranslationTable.

    Accepts weights as:

        {
            'A': {'GCT': 0.25, 'GCC': 0.25, ...},
            'R': {'CGT': 0.1, 'CGC': 0.2, ...},
            ...
        }

    Stores weights as:

        {
            'GCT': 0.25,
            'GCC': 0.25,
            ...
        }

    Weights are normalised per amino acid.
    """

    __slots__ = (
        '_locked',
        'table',
        'weights',
    )

    def __init__(
        self,
        weights: Mapping[str, Mapping[str, float]],
        table: Optional[TranslationTable] = None,
    ) -> None:
        object.__setattr__(self, '_locked', False)

        table = table or TranslationTable()

        normalised = self._normalise_weights(weights, table)
        self._validate_weights(normalised, table)

        flat = {
            codon: weight
            for codon_weights in normalised.values()
            for codon, weight in codon_weights.items()
        }

        object.__setattr__(self, 'table', table)
        object.__setattr__(self, 'weights', FrozenDict(flat))

        object.__setattr__(self, '_locked', True)

    def __setattr__(self, name: str, value: Any) -> None:
        if getattr(self, '_locked', False):
            raise AttributeError(f'{type(self).__name__} is immutable')
        object.__setattr__(self, name, value)

    def __repr__(self) -> str:
        return f'{type(self).__name__}(table={self.table!r})'

    def __getitem__(self, codon: str) -> float:
        codon = self.table.normalise_codon(codon, rna=self.table.rna)
        return self.weights[codon]

    def for_amino_acid(self, aa: str) -> Dict[str, float]:
        aa = aa.upper()
        return {
            codon: self.weights[codon]
            for codon in self.table.aa_to_codons[aa]
        }

    @classmethod
    def uniform(cls, table: Optional[TranslationTable] = None) -> 'CodonWeights':
        table = table or TranslationTable()

        weights = {}
        for aa, codons in table.aa_to_codons.items():
            weight = 1.0 / len(codons)
            weights[aa] = {codon: weight for codon in codons}

        return cls(weights, table=table)

    @classmethod
    def from_codon_weights(
        cls,
        codon_weights: Mapping[str, float],
        table: Optional[TranslationTable] = None,
    ) -> 'CodonWeights':
        table = table or TranslationTable()

        grouped: Dict[str, Dict[str, float]] = {}

        for codon, weight in codon_weights.items():
            codon = table.normalise_codon(codon, rna=table.rna)
            aa = table.codons_to_aa[codon]
            grouped.setdefault(aa, {})[codon] = weight

        return cls(grouped, table=table)

    @staticmethod
    def _normalise_weights(
        weights: Mapping[str, Mapping[str, float]],
        table: TranslationTable,
    ) -> Dict[str, Dict[str, float]]:
        normalised = {}

        for aa, codon_weights in weights.items():
            aa = aa.upper()

            total = sum(codon_weights.values())
            if total <= 0:
                raise ValueError(f'Weights for amino acid {aa} must sum to > 0')

            normalised[aa] = {
                table.normalise_codon(codon, rna=table.rna): float(weight) / total
                for codon, weight in codon_weights.items()
            }

        return normalised

    @staticmethod
    def _validate_weights(
        weights: Mapping[str, Mapping[str, float]],
        table: TranslationTable,
    ) -> None:
        expected_aas = set(table.aa_to_codons)
        actual_aas = set(weights)

        missing_aas = expected_aas - actual_aas
        extra_aas = actual_aas - expected_aas

        if missing_aas:
            raise ValueError(f'Missing weights for amino acids: {sorted(missing_aas)}')

        if extra_aas:
            raise ValueError(f'Unexpected amino acids in weights: {sorted(extra_aas)}')

        for aa, expected_codons in table.aa_to_codons.items():
            expected = set(expected_codons)
            actual = set(weights[aa])

            missing_codons = expected - actual
            extra_codons = actual - expected

            if missing_codons:
                raise ValueError(f'Missing codon weights for amino acid {aa!r}: {sorted(missing_codons)}')

            if extra_codons:
                raise ValueError(f'Unexpected codons for amino acid {aa!r}: {sorted(extra_codons)}')

            for codon, weight in weights[aa].items():
                if weight < 0:
                    raise ValueError(f'Weight for codon {codon!r} cannot be negative')