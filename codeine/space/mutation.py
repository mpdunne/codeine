from typing import Collection, Dict, FrozenSet, List, Generator, Optional, Set, Union

from codeine.constraints.mutations import MutationDistanceConstraint
from codeine.graph.base import CodonRestriction
from codeine.motifs.restriction import RestrictionSite
from codeine.space.coding import CodingSpace
from codeine.translation.tables import TranslationTable
from codeine.translation.weights import CodonWeights
from codeine.utils.display import format_forbidden_motifs, format_forbidden_motif,\
    format_count, format_restrictions, format_positions


class MutationSpace:
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

        self.view = space.view.copy()
        self._base_pins = dict(space.view.pinned_codons)

        self.cds = self._validate_cds(cds)

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

    def __getitem__(self, index: Union[int, slice]) -> Union[str, List[str]]:
        """
        Return one or more valid sequences.

        Parameters
        ----------
        index
            Zero-based sequence index or slice.

        Returns
        -------
        str or List[str]
            The indexed sequence, or a list of sequences for a slice.
        """
        return self.view[index]

    def __iter__(self) -> Generator[str, None, None]:
        """
        Iterate over all valid sequences in this mutation space.
        Be aware that "all valid sequences" can be astronomically many!

        Yields
        ----------
        All valid sequences in the coding space, in order.
        """
        yield from self.view

    def __contains__(self, seq: str) -> bool:
        """
        Does the given seq exist in this space?

        Returns
        ----------
        True if and only if this is a valid sequence in this space.
        """
        return seq in self.view

    def __repr__(self) -> str:
        molecule = 'RNA' if self.translation_table.rna else 'DNA'

        lines = [
            f'{type(self).__name__}',
            '',
            f'Translation table: {self.translation_table.table_id} ({self.translation_table.name})',
            f'Molecule type: {molecule}',
            '',
            f'Amino acid sequence ({len(self.aa_seq)} aa):',
            f'{self.aa_seq}',
            '',
            'Reference CDS:',
            self.cds,
            '',
        ]

        if self.codon_restrictions:
            lines += [
                'Codon restrictions:',
                *format_restrictions(
                    self.codon_restrictions,
                    label='restricted positions',
                    max_lines=4,
                ),
                '',
            ]

        if self.forbidden_motifs:
            motifs = self.forbidden_motifs

            if isinstance(motifs, (str, RestrictionSite)):
                motifs = [motifs]

            lines += [
                'Forbidden motifs:',
                *format_forbidden_motifs(
                    [
                        format_forbidden_motif(
                            motif,
                            rna=self.translation_table.rna,
                        )
                        for motif in motifs
                    ],
                    max_lines=4,
                ),
                '',
            ]

        if self.max_homopolymer is not None:
            lines += [
                'Maximum homopolymer length:',
                f'    {self.max_homopolymer}',
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
        else:
            lines.append(
                f'Num. valid variants: {format_count(self.n_valid_variants)}'
            )

        return '\n'.join(lines)

    def sample(self, n: Optional[int] = None) -> str:
        """
        Sample one or more variants from this mutation space.

        Parameters
        ----------
        n
            Number of sequences to sample. If omitted, return a single sequence.

        Returns
        -------
        A sampled string sequence from this mutation space.
        """
        return self.view.sample(n=n)

    def enumerate(self) -> Generator[str, None, None]:
        """
        Generate all sequences in this mutation space.

        Yields
        ------
        str
            A valid coding sequence.
        """
        yield from self.view.enumerate()

    def contains(self, seq: str) -> bool:
        """
        Check whether a coding sequence is contained in this mutation space.

        Parameters
        ----------
        seq
            The sequence to check

        Returns
        -------
        True if and only if the sequence is contained in this mutation space.
        """
        return self.view.contains(seq)

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

        Distances are measured from the reference CDS and can be either nucleotide (Hammming)
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
    def aa_seq(self) -> str:
        """
        The amino acid sequence for this mutation space.
        """
        return self.view.aa_seq

    @property
    def translation_table(self) -> TranslationTable:
        """
        The translation table from the underlying graph.
        """
        return self.view.translation_table

    @property
    def codon_weights(self) -> CodonWeights:
        """
        The codon weights from the underlying graph.
        """
        return self.view.codon_weights

    @property
    def codon_restrictions(self) -> Dict[int, CodonRestriction]:
        """
        The fixed codon restrictions from the underlying graph.
        """
        return self.view.codon_restrictions

    @property
    def context_l(self) -> str:
        """
        The left context sequence from the underlying graph.
        """
        return self.view.context_l

    @property
    def context_r(self) -> str:
        """
        The right context sequence from the underlying graph.
        """
        return self.view.context_r

    @property
    def pinned_codons(self) -> Dict[int, List[str]]:
        """
        Pins currently applied to the mutation view.

        This includes inherited pins and pins used internally to freeze positions.
        """
        return self.view.pinned_codons

    @property
    def n_valid_variants(self) -> int:
        """
        Number of valid variants under the current mutation constraints.
        """
        return self.view.n_valid_sequences

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
        Replace the mutation-distance constraint on the underlying view.
        """
        constraints = tuple(
            constraint
            for constraint in self.view.constraints
            if not isinstance(constraint, MutationDistanceConstraint)
        )

        if self.has_distance_constraints:
            constraints += (
                MutationDistanceConstraint(
                    reference_cds=self.cds,
                    min_nts=self.min_nts,
                    max_nts=self.max_nts,
                    min_codons=self.min_codons,
                    max_codons=self.max_codons,
                ),
            )

        self.view.set_constraints(constraints)

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
        cds = cds.upper()

        if not self.view.contains(cds):
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
