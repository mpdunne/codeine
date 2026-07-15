import pickle

from pathlib import Path
from typing import Dict, Iterator, List, Optional, Sequence, Tuple, Union, TYPE_CHECKING

if TYPE_CHECKING:
    from codeine.space.mutation import MutationSpace

from codeine.constraints.banned import BannedSequenceConstraint
from codeine.constraints.base import Constraint
from codeine.graph.base import CodonGraph, CodonRestriction
from codeine.motifs.validate import expand_and_validate_sequence_constraints, ForbiddenMotifs
from codeine.motifs.restriction import RestrictionSite
from codeine.translation.tables import TranslationTable
from codeine.translation.weights import CodonWeights
from codeine.utils.display import format_forbidden_motifs, format_forbidden_motif,\
    format_count, format_restrictions
from codeine.utils.sampling import Seedable


class CodingSpace:
    """
    Represents a space of valid coding sequences for a protein under constraints.
    """
    def __init__(
        self,
        aa_seq: str,
        *,
        translation_table: Optional[TranslationTable] = None,
        rna: Optional[bool] = None,
        codon_restrictions: Optional[Dict[int, CodonRestriction]] = None,
        forbidden_motifs: Optional[ForbiddenMotifs] = None,
        max_homopolymer: Optional[int] = None,
        constraints: Optional[Sequence[Constraint]] = None,
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
        codon_restrictions
            Any codon restrictions in the format e.g. ``{4: 'TCC'}`` or ``{5: ['AGT', 'AGC']}``. Positions are 1-based.
        forbidden_motifs
            Forbidden motifs, either as strings or as ``codeine.RestrictionSite``.
        max_homopolymer
            The maximum allowed length of nucleotide homopolymer
        constraints
            Any additional constraints.
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
            codon_restrictions=codon_restrictions,
            translation_table=translation_table,
            context_l=context_l,
            context_r=context_r,
        )

        view = graph.view(seed=seed, weights=codon_weights)
        self.view = view

        self.constraints = tuple(constraints or ())
        self.forbidden_motifs = self._normalise_forbidden_motifs(forbidden_motifs)
        self.max_homopolymer = max_homopolymer

        self._update_constraints()

    @classmethod
    def load(cls, path) -> 'CodingSpace':
        """
        Load a coding space from disc.
        """
        with Path(path).open('rb') as f:
            return pickle.load(f)

    def save(self, path) -> None:
        """
        Save this coding space to disc.
        """
        with Path(path).open('wb') as f:
            pickle.dump(self, f)

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

    def __iter__(self) -> Iterator[str]:
        """
        Iterate over all valid sequences in this coding space.
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
                        format_forbidden_motif(motif, rna=self.translation_table.rna)
                        for motif in motifs
                    ],
                    max_lines=4,
                ),
                '',
            ]

        if self.max_homopolymer is not None:
            lines += [
                'Maximum homopolymer length:',
                f'  {self.max_homopolymer}',
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

    def sample(self, n: Optional[int] = None) -> Union[str, List[str]]:
        """
        Sample one or more variants from this coding space.

        Parameters
        ----------
        n
            Number of sequences to sample. If omitted, return a single sequence.

        Returns
        -------
        A sampled sequence, or a list of sampled sequences.
        """
        return self.view.sample(n=n)

    def enumerate(self) -> Iterator[str]:
        """
        Generate all sequences in this space. If there are many (and often there are
        astronomically many), one would not expect to reach the 'end'. However for smaller
        sequence spaces, such as mutation spaces, it's quite possible to get there.

        Yields
        ------
        str
            A valid coding sequence.
        """
        yield from self.view.enumerate()

    def contains(self, seq: str) -> bool:
        """
        Check whether a coding sequence is contained in this coding space.

        Parameters
        ----------
        seq
            The sequence to check.

        Returns
        -------
        True if and only if the sequence is contained in this coding space.
        """
        return self.view.contains(seq)

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

    def set_forbidden_motifs(self, forbidden_motifs: Optional[ForbiddenMotifs]) -> None:
        """
        Set the forbidden motifs for this coding space.

        Parameters
        ----------
        forbidden_motifs
            Motifs that should be forbidden in generated sequences.
        """
        self.forbidden_motifs = self._normalise_forbidden_motifs(forbidden_motifs)
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

    def set_constraints(self, constraints: Sequence[Constraint]) -> None:
        """
        Set any additional constraints.

        Parameters
        ----------
        constraints
            The constraints to set, as ``Constraint`` objects.
        """
        self.constraints = tuple(constraints)
        self._update_constraints()

    def clear_constraints(self) -> None:
        """
        Remove custom constraints.
        """
        self.set_constraints(())

    def set_codon_weights(self, codon_weights: CodonWeights) -> None:
        """
        Set the codon weights used when sampling from this coding space.
        """
        self.view.set_weights(codon_weights)

    @property
    def n_valid_sequences(self) -> int:
        """
        The number of valid sequences in this space.
        """
        return self.view.n_valid_sequences

    @property
    def aa_seq(self) -> str:
        """
        The amino acid sequence for this coding space.
        """
        return self.view.aa_seq

    @property
    def translation_table(self) -> TranslationTable:
        """
        The translation table being used in this space.
        """
        return self.view.translation_table

    @property
    def codon_weights(self) -> CodonWeights:
        """
        The codon weights being used in this space.
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
        Temporary codon pins currently applied to this coding space.
        """
        return self.view.pinned_codons

    def _update_constraints(self) -> None:
        """
        Rebuild and apply all constraints for this coding space.
        """
        forbidden_sequences = expand_and_validate_sequence_constraints(
            forbidden_motifs=self.forbidden_motifs,
            max_homopolymer=self.max_homopolymer,
            rna=self.translation_table.rna,
        )

        constraints = self.constraints

        if forbidden_sequences:
            constraints = (
                BannedSequenceConstraint(forbidden_sequences),
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

    def _normalise_forbidden_motifs(
            self,
            forbidden_motifs: Optional[ForbiddenMotifs],
    ) -> Optional[ForbiddenMotifs]:
        """
        Normalise string forbidden motifs to the molecule type used by this coding
        space. RestrictionSite objects are left unchanged.
        """
        if forbidden_motifs is None:
            return None

        if isinstance(forbidden_motifs, str):
            return self.translation_table.normalise_sequence(forbidden_motifs)

        if isinstance(forbidden_motifs, RestrictionSite):
            return forbidden_motifs

        return [
            self.translation_table.normalise_sequence(motif)
            if isinstance(motif, str)
            else motif
            for motif in forbidden_motifs
        ]
