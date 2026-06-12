from typing import Dict, Generator, List, Optional, Sequence, Union

from codeine.motifs.restriction import RestrictionSite
from codeine.sequence.display import format_banned_sequences, format_count, format_restrictions
from codeine.sequence.graph import CodonGraph
from codeine.translation.tables import TranslationTable
from codeine.translation.weights import CodonWeights


ForbiddenMotif = Union[str, RestrictionSite]


class CodingSpace:
    def __init__(
        self,
        aa_seq: str,
        codon_restrictions: Optional[Dict[int, str]] = None,
        forbidden_motifs: Optional[Sequence[ForbiddenMotif]] = None,
        translation_table: TranslationTable = None,
        codon_weights: CodonWeights = None,
        context_l: str = '',
        context_r: str = '',
    ) -> None:
        if translation_table is None:
            rna = codon_weights.rna if codon_weights is not None else False
        else:
            rna = translation_table.rna

        self.forbidden_sequences = self._expand_and_validate_forbidden_motifs(forbidden_motifs, rna=rna)

        self.view = CodonGraph(
            aa_seq,
            codon_restrictions=codon_restrictions,
            banned_sequences=self.forbidden_sequences,
            translation_table=translation_table,
            weights=codon_weights,
            context_l=context_l,
            context_r=context_r,
        ).view()

    @staticmethod
    def _expand_and_validate_forbidden_motifs(forbidden_motifs: Sequence[Union[str, ForbiddenMotif]], rna: bool) -> List[str]:
        all_sequences = []

        if forbidden_motifs is None:
            return []

        if isinstance(forbidden_motifs, (str, RestrictionSite)):
            forbidden_motifs = [forbidden_motifs]

        for motif in forbidden_motifs:
            if isinstance(motif, RestrictionSite):
                sequences = [*motif.motifs]

            elif isinstance(motif, str):
                if len(motif) == 0:
                    raise ValueError('Forbidden motifs cannot be empty.')

                sequences = [motif]

            else:
                raise TypeError('Forbidden motifs must be strings or codeine.RestrictionSite.')

            sequences = [seq.upper() for seq in sequences]
            sequences = [seq.replace('T', 'U') if rna else seq.replace('U', 'T') for seq in sequences]

            all_sequences += sequences

        return sorted(set(all_sequences))

    @classmethod
    def from_graph(cls, graph: CodonGraph) -> 'CodingSpace':
        return cls.from_view(graph.view())

    @classmethod
    def from_view(cls, view) -> 'CodingSpace':
        obj = cls.__new__(cls)
        obj.view = view
        obj.forbidden_motifs = []
        obj.forbidden_sequences = list(view.graph.banned_sequences)
        return obj

    def __getitem__(self, index: int) -> str:
        return self.view[index]

    def __iter__(self) -> Generator[str, None, None]:
        yield from self.view

    def __contains__(self, seq: str) -> bool:
        return seq in self.view

    def __repr__(self) -> str:
        molecule = 'RNA' if self.translation_table.rna else 'DNA'

        lines = [
            f'{type(self).__name__}',
            '',
            f'Translation table: {self.translation_table.table_id} ({self.translation_table.name})',
            f'Molecule type: {molecule}',
            '',
            f'Amino acid sequence ({len(self.view.aa_seq)} aa)',
            f'{self.view.aa_seq}',
            '',
        ]

        if self.view.graph.codon_restrictions:
            lines += [
                'Codon restrictions:',
                *format_restrictions(
                    self.view.graph.codon_restrictions,
                    label='restricted positions',
                    max_lines=4,
                ),
                '',
            ]

        if self.forbidden_motifs:
            lines += [
                'Forbidden motifs:',
                *format_banned_sequences(
                    self.forbidden_sequences,
                    max_lines=4,
                ),
                '',
            ]

        elif self.forbidden_sequences:
            lines += [
                'Forbidden motifs:',
                *format_banned_sequences(
                    self.forbidden_sequences,
                    max_lines=4,
                ),
                '',
            ]

        if self.view.pinned_codons:
            lines += [
                'Temporary pins:',
                *format_restrictions(
                    self.view.pinned_codons,
                    label='pinned positions',
                    max_lines=4,
                ),
                '',
            ]

        lines.append(f'Num. valid coding sequences: {format_count(self.n_valid_sequences)}')

        return '\n'.join(lines)

    def sample(self) -> str:
        return self.view.sample()

    def pin_codons(self, pinned_codons):
        self.view.pin_codons(pinned_codons)

    def unpin_codons(self, positions):
        self.view.unpin_codons(positions)

    def clear_pins(self):
        self.view.clear_pins()

    def contains(self, seq: str) -> bool:
        return self.view.contains(seq)

    @property
    def n_valid_sequences(self) -> int:
        return self.view.n_valid_sequences

    @property
    def translation_table(self) -> TranslationTable:
        return self.view.graph.tt

    @property
    def codon_weights(self) -> CodonWeights:
        return self.view.graph.cw

    def enumerate(self) -> Generator[str, None, None]:
        yield from self.view.enumerate()

    def mutants(
        self,
        seq: str,
        positions: Sequence[int],
    ) -> 'CodingSpace':
        seq = seq.upper()

        if not self.contains(seq):
            raise ValueError('Parent sequence is not contained in this coding space.')

        positions = set(positions)

        if any(pos < 1 or pos > len(self.view.aa_seq) for pos in positions):
            raise ValueError('Mutation positions out of range.')

        mutation_pins = {}

        for pos in range(1, len(self.view.aa_seq) + 1):
            if pos not in positions:
                start = (pos - 1) * 3
                mutation_pins[pos] = seq[start:start + 3]

        view = self.view.copy()
        view.pin_codons(mutation_pins)

        return CodingSpace.from_view(view)