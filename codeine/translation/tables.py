import re

from Bio.Data import CodonTable
from typing import Any

from codeine.utils.dict import FrozenDict


class TranslationTable:
    """
    Translation table. Here we use the NCBI translation table IDs
    (see https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi)
    Uses DNA by default, but can switch to RNA by setting rna=True.
    """

    def __init__(self, table_id: int = 1, rna: bool = False) -> None:
        """
        Constructor for the TranslationTable class.

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

        self.name = biopython_table.names[0]

        dna_to_aa = biopython_table.forward_table
        dna_to_aa = {**dna_to_aa, **{stop: '*' for stop in biopython_table.stop_codons}}

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

    def __setattr__(self, name: str, value: Any) -> None:
        if getattr(self, '_locked', False) and name != '_locked':
            raise AttributeError(f'{type(self).__name__} is immutable')
        object.__setattr__(self, name, value)

    def __repr__(self) -> str:
        molecule = 'RNA' if self.rna else 'DNA'

        lines = [
            f'TranslationTable',
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
        return self.codons_to_aa[codon]

    def normalise_sequence(self, seq: str) -> str:
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
