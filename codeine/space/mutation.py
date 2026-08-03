from typing import Collection, FrozenSet, Optional, Set

from codeine.constraints.mutations import MutationDistanceConstraint
from codeine.space.base import Space
from codeine.space.coding import CodingSpace
from codeine.utils.display import format_constraints, format_count, format_restrictions,\
    format_positions, format_sequence


class MutationSpace(Space):
    """
    Represents the subset of valid coding sequences reachable from a reference CDS
    by mutation under constraints.

    A ``MutationSpace`` is defined by:

    - A ``CodingSpace`` containing the global sequence constraints.
    - A reference CDS.
    - A set of codon positions that are free to mutate.

    Positions that are not free are temporarily considered frozen and
    will remain identical to the reference CDS.
    """

    def __init__(self,
                 space: CodingSpace,
                 cds: str,
                 *,
                 free_positions: Optional[Collection[int]] = None,
                 min_nts: Optional[int] = None,
                 max_nts: Optional[int] = None,
                 min_codons: Optional[int] = None,
                 max_codons: Optional[int] = None,
                 ):
        """
        Parameters
        ----------
        space
            The ``CodingSpace`` object to which this given sequence should belong.
        cds
            The parent/reference CDS.
        free_positions
            Which positions are allowed to change?
        min_nts
            Minimum number of nucleotide differences from the reference CDS.
        max_nts
            Maximum number of nucleotide differences from the reference CDS.
        min_codons
            Minimum number of codon differences from the reference CDS.
        max_codons
            Maximum number of codon differences from the reference CDS.
        """
        self.forbidden_motifs = space.forbidden_motifs
        self.max_homopolymer = space.max_homopolymer

        self._base_view = space.view.copy()
        self._base_pins = dict(self._base_view.pinned_codons)

        self.cds = self._validate_cds(cds)

        self.view = self._base_view.copy()

        if free_positions is None:
            free_positions = range(1, len(self.view.aa_seq) + 1)

        self._free_positions: Set[int] = set()
        self.set_free_positions(free_positions)

        self.min_nts: Optional[int] = None
        self.max_nts: Optional[int] = None
        self.min_codons: Optional[int] = None
        self.max_codons: Optional[int] = None

        self.set_distance_constraints(
            min_nts=min_nts,
            max_nts=max_nts,
            min_codons=min_codons,
            max_codons=max_codons,
        )

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
            f'Reference CDS ({len(self.cds)} nt):',
            *format_sequence(self.cds, max_lines=8),
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

        if self._base_view.constraints:
            lines += [
                'Constraints:',
                *format_constraints(
                    self._base_view.constraints,
                    rna=self.translation_table.rna,
                    max_motifs=4,
                ),
                '',
            ]

        if self._base_pins:
            lines += [
                'Inherited pins:',
                *format_restrictions(
                    self._base_pins,
                    label='pinned positions',
                    max_lines=4,
                ),
                '',
            ]

        lines += [
            'Free positions:',
            f'    {format_positions(self.free_positions)}',
            '',
        ]

        if self.has_distance_constraints:
            lines += [
                'Mutation distance:',
                f'    nts: {self._format_distance(self.min_nts, self.max_nts)}',
                f'    codons: {self._format_distance(self.min_codons, self.max_codons)}',
                '',
            ]

        lines.append(f'Num. valid variants: {format_count(self.n_valid_sequences)}')

        return '\n'.join(lines)

    def set_free_positions(self, positions: Collection[int]) -> None:
        """
        Replace the current set of free positions.
        """
        self._free_positions = self._validate_positions(positions)
        self._update_pins()

    def freeze_positions(self, positions: Collection[int]) -> None:
        """
        Freeze the given codon positions.
        """
        positions = self._validate_positions(positions)

        self._free_positions -= positions
        self._update_pins()

    def unfreeze_positions(self, positions: Collection[int]) -> None:
        """
        Unfreeze the given codon positions.
        """
        positions = self._validate_positions(positions)

        self._free_positions |= positions
        self._update_pins()

    def freeze_all(self) -> None:
        """
        Freeze all codon positions.
        """
        self._free_positions.clear()
        self._update_pins()

    def unfreeze_all(self) -> None:
        """
        Unfreeze all codon positions.
        """
        self._free_positions = set(range(1, len(self.view.aa_seq) + 1))
        self._update_pins()

    def set_distance_constraints(self,
                                 min_nts: Optional[int] = None,
                                 max_nts: Optional[int] = None,
                                 min_codons: Optional[int] = None,
                                 max_codons: Optional[int] = None,
                                 ) -> None:
        """
        Set mutation distance constraints.

        Distances are measured from the reference CDS and can be either nucleotide (Hamming)
        distances, i.e. the number of nucleotides that are different from the reference CDS,
        or codon distances, i.e. the number of codons that are different.
        """

        self._validate_distance('min_nts', min_nts)
        self._validate_distance('max_nts', max_nts)
        self._validate_distance('min_codons', min_codons)
        self._validate_distance('max_codons', max_codons)

        if min_nts is not None and max_nts is not None and min_nts > max_nts:
            raise ValueError('min_nts cannot be greater than max_nts.')

        if min_codons is not None and max_codons is not None and min_codons > max_codons:
            raise ValueError('min_codons cannot be greater than max_codons.')

        self.min_nts = min_nts
        self.max_nts = max_nts
        self.min_codons = min_codons
        self.max_codons = max_codons

        self._update_distance_constraint()

    def clear_distance_constraints(self) -> None:
        """
        Remove all mutation distance constraints.
        """
        self.set_distance_constraints()

    @property
    def free_positions(self) -> FrozenSet[int]:
        """
        Codon positions that are currently free to mutate.
        """
        return frozenset(self._free_positions)

    @property
    def frozen_positions(self) -> FrozenSet[int]:
        """
        Codon positions that are currently fixed to the reference CDS.
        """
        all_positions = set(range(1, len(self.view.aa_seq) + 1))
        return frozenset(all_positions - self._free_positions)

    @property
    def has_distance_constraints(self) -> bool:
        """
        Whether this mutation space has mutation distance constraints.
        """
        distance_constraints = [self.min_nts, self.max_nts, self.min_codons, self.max_codons]
        return any(value is not None for value in distance_constraints)

    def _update_distance_constraint(self) -> None:
        """
        Rebuild the mutation layer over the compiled base view.
        """
        pins = self.view.pinned_codons

        self.view = self._base_view.copy()
        self.view.set_pinned_codons(pins)

        if self.has_distance_constraints:
            self.view.add_constraints(
                MutationDistanceConstraint(
                    reference_cds=self.cds,
                    min_nts=self.min_nts,
                    max_nts=self.max_nts,
                    min_codons=self.min_codons,
                    max_codons=self.max_codons,
                )
            )

    def _validate_positions(self, positions: Collection[int]) -> Set[int]:
        """
        Validate a collection of codon positions.
        """
        positions = set(positions)
        invalid = [pos for pos in positions if pos < 1 or pos > len(self.view.aa_seq)]
        if invalid:
            raise ValueError(f'Invalid codon positions: {sorted(invalid)}')

        return positions

    @staticmethod
    def _validate_distance(name: str, value: Optional[int]) -> None:
        """
        Validate a mutation distance value.
        """
        if value is None:
            return

        if not isinstance(value, int):
            raise TypeError(f'{name} must be an integer.')

        if value < 0:
            raise ValueError(f'{name} must be non-negative.')

    @staticmethod
    def _format_distance(minimum: Optional[int], maximum: Optional[int]) -> str:
        """
        Format a distance range for display.
        """
        if minimum is None and maximum is None:
            return 'any'

        if minimum == maximum:
            return str(minimum)

        if minimum is None:
            return f'up to {maximum}'

        if maximum is None:
            return f'at least {minimum}'

        return f'{minimum}..{maximum}'

    def _validate_cds(self, cds: str) -> str:
        """
        Check that the CDS belongs to the underlying space.

        Parameters
        ----------
        cds
            The inputted CDS.

        Returns
        -------
        A normalised and validated version of the inputted CDS.
        """
        cds = self._base_view.translation_table.normalise_sequence(cds)

        if not self._base_view.contains(cds):
            raise ValueError('CDS is not contained in this coding space.')

        return cds

    def _codon_at_position(self, pos: int) -> str:
        """
        Get the codon of the reference CDS at the specified position.

        Parameters
        ----------
        pos
            The position in the AA sequence.

        Returns
        -------
        A codon.
        """
        return self.cds[3 * (pos - 1): 3 * pos]

    def _update_pins(self) -> None:
        """
        Update the pins on the underlying view from inherited and frozen pins.
        """
        frozen_pins = {pos: self._codon_at_position(pos) for pos in self.frozen_positions}
        pins = {**self._base_pins, **frozen_pins}

        if pins == self.view.pinned_codons:
            return

        self.view.set_pinned_codons(pins)
