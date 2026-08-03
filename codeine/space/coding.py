from typing import Dict, Optional, Sequence, Tuple, Union, TYPE_CHECKING

if TYPE_CHECKING:
    from codeine.space.mutation import MutationSpace

from codeine.constraints.motifs import ForbiddenMotifConstraint, Motifs
from codeine.constraints.homopolymers import HomopolymerConstraint
from codeine.constraints.base import Constraint
from codeine.graph.base import CodonGraph, CodonRestriction
from codeine.space.base import Space
from codeine.translation.tables import TranslationTable
from codeine.translation.weights import CodonWeights
from codeine.utils.display import format_constraints, format_count, format_restrictions,\
    format_sequence
from codeine.utils.sampling import Seedable


class CodingSpace(Space):
    """
    Represents a space of valid coding sequences for a protein under constraints.
    """
    def __init__(
        self,
        aa_seq: str,
        *,
        translation_table: Optional[TranslationTable] = None,
        rna: Optional[bool] = None,
        fixed_codons: Optional[Dict[int, CodonRestriction]] = None,
        constraints: Optional[Sequence[Constraint]] = None,
        forbidden_motifs: Optional[Motifs] = None,
        max_homopolymer: Optional[int] = None,
        context_l: str = '',
        context_r: str = '',
        codon_weights: Optional[CodonWeights] = None,
        seed: Optional[Seedable] = None,
    ) -> None:
        """
        Parameters
        ----------
        aa_seq
            The amino acid sequence.
        translation_table
            The translation table to use. Leave blank to use standard table.
        rna
            Whether to use RNA. If false or blank, use DNA.
        fixed_codons
            Any codon restrictions in the format e.g. ``{4: 'TCC'}`` or ``{5: ['AGT', 'AGC']}``. Positions are 1-based.
        constraints
            Constraints to apply to this coding space.
        forbidden_motifs
            Forbidden motifs, either as strings or as ``codeine.RestrictionSite``.
        max_homopolymer
            The maximum allowed length of nucleotide homopolymer
        context_l
            The context sequence to the left of the coding sequence.
        context_r
            The context sequence to the right of the coding sequence.
        codon_weights
            The codon weights to use. Leave blank to sample uniformly.
        seed
            Seed used to initialise the random number generator for sampling.
        """

        translation_table, codon_weights = self._resolve_tables(translation_table, codon_weights, rna)

        graph = CodonGraph(
            aa_seq,
            fixed_codons=fixed_codons,
            translation_table=translation_table,
            context_l=context_l,
            context_r=context_r,
        )

        view = graph.view(seed=seed, weights=codon_weights)
        self.view = view

        self.constraints = tuple(constraints or ())
        self.forbidden_motifs = forbidden_motifs
        self.max_homopolymer = max_homopolymer

        self._update_constraints()

    def __repr__(self) -> str:
        molecule = 'RNA' if self.translation_table.rna else 'DNA'

        lines = [
            f'{type(self).__name__}',
            '',
            f'Translation table: {self.translation_table.table_id} ({self.translation_table.name})',
            f'Molecule type: {molecule}',
            '',
            f'Amino acid sequence ({len(self.aa_seq)} aa):',
            *format_sequence(self.aa_seq, max_lines=6),
            '',
        ]

        if self.fixed_codons:
            lines += [
                'Codon restrictions:',
                *format_restrictions(
                    self.fixed_codons,
                    label='restricted positions',
                    max_lines=4,
                ),
                '',
            ]

        if self.context_l:
            lines += [
                f'Left context ({len(self.context_l)} nt):',
                *format_sequence(self.context_l, max_lines=2),
                '',
            ]

        if self.context_r:
            lines += [
                f'Right context ({len(self.context_r)} nt):',
                *format_sequence(self.context_r, max_lines=2),
                '',
            ]

        if self.view.constraints:
            lines += [
                'Constraints:',
                *format_constraints(
                    self.view.constraints,
                    rna=self.translation_table.rna,
                    max_motifs=4,
                ),
                '',
            ]

        if self.pinned_codons:
            lines += [
                'Temporary pins:',
                *format_restrictions(
                    self.pinned_codons,
                    label='pinned positions',
                    max_lines=4,
                ),
                '',
            ]

        lines.append(f'Num. valid coding sequences: {format_count(self.n_valid_sequences)}')

        return '\n'.join(lines)

    def mutants(
            self,
            cds: str,
            free_positions: Optional[Sequence[int]] = None,
            min_nts: Optional[int] = None,
            max_nts: Optional[int] = None,
            min_codons: Optional[int] = None,
            max_codons: Optional[int] = None,
    ) -> 'MutationSpace':
        """
        Return a space of mutants relative to a given coding sequence, i.e. a space derived
        from this one but which fixes the sequence on all but the specified positions.

        Parameters
        ----------
        cds
            The sequence to mutate.
        free_positions
            The positions that are allowed to vary.
        min_nts
            The min nucleotide (Hamming) distance relative to the reference sequence.
        max_nts
            The max nucleotide (Hamming) distance relative to the reference sequence.
        min_codons
            The min number of changed codons relative to the reference sequence.
        max_codons
            The max number of changed codons relative to the reference sequence.
        """
        cds = self.translation_table.normalise_sequence(cds)

        if not self.contains(cds):
            raise ValueError('CDS is not contained in this coding space.')

        from codeine.space.mutation import MutationSpace

        return MutationSpace(
            space=self,
            cds=cds,
            free_positions=free_positions,
            min_nts=min_nts,
            max_nts=max_nts,
            min_codons=min_codons,
            max_codons=max_codons,
        )

    def pin_codons(self, pinned_codons: Dict[int, str]) -> None:
        """
        Pin temporary codons in this coding space.

        Parameters
        ----------
        pinned_codons
            A dict specifying which codons to pin, by position.
        """
        self.view.pin_codons(pinned_codons)

    def unpin_codons(self, positions: Sequence[int]) -> None:
        """
        Remove temporary codon pins by position.

        Parameters
        ----------
        positions
            Positions to unpin.
        """
        self.view.unpin_codons(positions)

    def set_pinned_codons(self, pinned_codons: Dict[int, str]) -> None:
        """
        Replace all temporary codon pins on this coding space.

        Parameters
        ----------
        pinned_codons
            A dict specifying which codons to pin, by position.
        """
        self.view.set_pinned_codons(pinned_codons)

    def clear_pins(self) -> None:
        """
        Remove all temporary codon pins from this coding space.
        """
        self.view.clear_pins()

    def set_forbidden_motifs(self, forbidden_motifs: Optional[Motifs]) -> None:
        """
        Set the forbidden motifs for this coding space.

        Parameters
        ----------
        forbidden_motifs
            Motifs that should be forbidden in generated sequences.
        """
        self.forbidden_motifs = forbidden_motifs
        self._update_constraints()

    def clear_forbidden_motifs(self) -> None:
        """
        Remove all forbidden motifs from this coding space.
        """
        self.set_forbidden_motifs(None)

    def set_max_homopolymer(self, max_homopolymer: Optional[int]) -> None:
        """
        Set the maximum allowed homopolymer length.

        Parameters
        ----------
        max_homopolymer
            The longest allowed repeated run of one nucleotide, or None for no limit.
        """
        self.max_homopolymer = max_homopolymer
        self._update_constraints()

    def clear_max_homopolymer(self) -> None:
        """
        Remove the maximum homopolymer constraint from this coding space.
        """
        self.set_max_homopolymer(None)

    def add_constraints(self, constraints: Union[Constraint, Sequence[Constraint]]) -> None:
        """
        Add one or more constraints to this coding space.

        Parameters
        ----------
        constraints
            Constraints to add.
        """
        if isinstance(constraints, Constraint):
            constraints = [constraints]

        constraints = tuple(constraints)

        if not constraints:
            return

        self.constraints += constraints
        self.view.add_constraints(constraints)

    def set_constraints(self, constraints: Sequence[Constraint]) -> None:
        """
        Replace the additional constraints for this coding space.

        Parameters
        ----------
        constraints
            The constraints to set, as ``Constraint`` objects.
        """
        self.constraints = tuple(constraints)
        self._update_constraints()

    def clear_constraints(self) -> None:
        """
        Remove all additional constraints. Constraints configured through
        ``forbidden_motifs`` and ``max_homopolymer`` are unaffected.
        """
        self.set_constraints(())

    def _update_constraints(self) -> None:
        """
        Rebuild and apply all constraints for this coding space.
        """
        constraints = self.constraints

        if self.forbidden_motifs is not None:
            constraints = (
                ForbiddenMotifConstraint(self.forbidden_motifs),
                *constraints,
            )

        if self.max_homopolymer is not None:
            constraints = (
                HomopolymerConstraint(self.max_homopolymer),
                *constraints,
            )

        self.view.set_constraints(constraints)

    @staticmethod
    def _resolve_tables(
            translation_table: Optional[TranslationTable],
            codon_weights: Optional[CodonWeights],
            rna: Optional[bool],
    ) -> Tuple[TranslationTable, CodonWeights]:
        """
        Resolve user-submitted (or not) translation table, codon weights, and RNA flag.
        """

        if translation_table is None:
            translation_table = TranslationTable(table_id=1, rna=False if rna is None else rna)

        elif rna is not None and translation_table.rna != rna:
            raise ValueError('Value for rna is inconsistent with the provided translation table.')

        if codon_weights is None:
            codon_weights = CodonWeights.uniform(table=translation_table)
        else:
            codon_weights = codon_weights.for_table(translation_table)

        return translation_table, codon_weights
