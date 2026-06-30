import json
import re

from pathlib import Path
from typing import Any, Dict

from codeine.utils.dict import FrozenDict


class TranslationTable:
    """
    Translation table. Here we use the NCBI translation table IDs
    (see https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi)
    Uses DNA by default, but can switch to RNA by setting ``rna=True``.
    """

    def __init__(self, table_id: int = 1, rna: bool = False) -> None:
        """
        Parameters
        ----------
        table_id
            Which translation table to use. Default is 1, which is the standard genetic code.
        rna
            Whether to use RNA. Default is no (False), i.e. use DNA.
        """
        self._locked = False

        self.table_id = table_id
        self.rna = rna

        tables = self._load_tables()

        try:
            table = tables[str(table_id)]
        except KeyError:
            raise ValueError(f'Unknown NCBI translation table ID: {table_id}.')

        self.name = table['name']
        self.names = tuple(table.get('names', (self.name,)))

        self.stop_codons = tuple(
            self.normalise_sequence(codon)
            for codon in table['stop_codons']
        )

        dna_to_aa = table['codon_to_aa']

        aa_to_dna = {}
        for codon, aa in dna_to_aa.items():
            aa_to_dna.setdefault(aa, []).append(codon)

        codons_to_aa = {
            self.normalise_sequence(codon): aa
            for codon, aa in dna_to_aa.items()
        }

        aa_to_codons = {
            aa: tuple(
                self.normalise_sequence(codon)
                for codon in codons
            )
            for aa, codons in aa_to_dna.items()
        }

        self.codons_to_aa = FrozenDict(codons_to_aa)
        self.aa_to_codons = FrozenDict(aa_to_codons)

        self._locked = True

    @classmethod
    def custom(
            cls,
            codons_to_aa: Dict[str, str],
            name: str = 'Custom',
            rna: bool = False,
    ) -> 'TranslationTable':
        """
        Create a custom translation table.

        This is useful for synthetic, partial, experimental, or externally
        defined translation tables, or for tables that have been discovered
        but not yet added to Codeine. Codons are normalised to DNA or RNA
        according to the ``rna`` argument.

        Parameters
        ----------
        codons_to_aa
            Mapping from codons to amino-acid symbols.
        name
            Name of the custom translation table. Default is ``'Custom'``.
        rna
            Whether the table should use RNA codons.

        Returns
        -------
        TranslationTable
            A custom translation table.
        """
        obj = cls.__new__(cls)
        obj._locked = False

        obj.table_id = None
        obj.rna = rna
        obj.name = name
        obj.names = (name,)

        normalised = {
            obj.normalise_sequence(codon): aa
            for codon, aa in codons_to_aa.items()
        }

        if any(len(codon) != 3 for codon in normalised):
            raise ValueError('All codons must have length 3.')

        if any(len(aa) != 1 for aa in normalised.values()):
            raise ValueError('Amino-acid symbols must be single characters.')

        aa_to_codons = {}
        for codon, aa in normalised.items():
            aa_to_codons.setdefault(aa, []).append(codon)

        obj.stop_codons = tuple(codon for codon, aa in normalised.items() if aa == '*')

        obj.codons_to_aa = FrozenDict(normalised)
        obj.aa_to_codons = FrozenDict({aa: tuple(codons) for aa, codons in aa_to_codons.items()})

        obj._locked = True
        return obj

    def __setattr__(self, name: str, value: Any) -> None:
        if getattr(self, '_locked', False) and name != '_locked':
            raise AttributeError(f'{type(self).__name__} is immutable')
        object.__setattr__(self, name, value)

    def __repr__(self) -> str:
        molecule = 'RNA' if self.rna else 'DNA'

        lines = [
            'TranslationTable',
            f'Table ID: {self.table_id} ({self.name})',
            f'Molecule type: {molecule}',
            '',
            'Table:',
        ]

        for aa in sorted(self.aa_to_codons):
            codons = ' '.join(self.aa_to_codons[aa])
            lines.append(f'    {aa}: {codons}')

        return '\n'.join(lines)

    def __getitem__(self, codon: str) -> str:
        try:
            codon = self.normalise_sequence(codon)
            return self.codons_to_aa[codon]

        except (AttributeError, TypeError, ValueError, KeyError) as e:
            raise KeyError(f'Invalid codon: {codon}.') from e

    def normalise_sequence(self, seq: str) -> str:
        """
        Format a sequence in the format specified by this codon table, i.e. convert to
        RNA/DNA and cast to upper case.

        Parameters
        ----------
        codon
            The inputted codon.

        Returns
        -------
        The normalised codon.
        """
        seq = seq.upper().replace(' ', '')
        seq = seq.replace('T', 'U') if self.rna else seq.replace('U', 'T')
        regex = r'^[ACGU]*$' if self.rna else r'^[ACGT]*$'
        if not re.match(regex, seq):
            raise ValueError('Sequence to normalise must be a nucleotide sequence.')
        else:
            return seq

    def translate(self, seq: str) -> str:
        """
        Translate a DNA/RNA coding sequence into its amino-acid sequence.
        """
        if len(seq) % 3 != 0:
            raise ValueError('Sequence length must be a multiple of 3')

        seq = self.normalise_sequence(seq)

        return ''.join(self.codons_to_aa[seq[i:i + 3]] for i in range(0, len(seq), 3))

    @staticmethod
    def _load_tables() -> dict:
        """
        Load tables from stored JSON entries.
        """
        path = Path(__file__).parent / 'data' / 'tables.json'
        with path.open() as f:
            return json.load(f)
