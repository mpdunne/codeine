from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from codeine.space.mutation import MutationSpace

import pickle
import random

from pathlib import Path
from typing import Dict, Generator, List, Optional, Sequence, Union

from codeine.motifs.restriction import RestrictionSite
from codeine.motifs.constraints import expand_and_validate_sequence_constraints, ForbiddenMotifs
from codeine.utils.display import format_banned_sequences, format_forbidden_motif,\
    format_count, format_restrictions
from codeine.graph.base import CodonGraph
from codeine.translation.tables import TranslationTable
from codeine.translation.weights import CodonWeights
from codeine.utils.sampling import Seedable


class CodingSpace:
    """
    Class representing coding sequence space, for sampling and mutating CDS coding sequences.
    """
    def __init__(
        self,
        aa_seq: str,
        codon_restrictions: Optional[Dict[int, str]] = None,
        forbidden_motifs: ForbiddenMotifs = None,
        max_homopolymer: Optional[int] = None,
        translation_table: TranslationTable = None,
        codon_weights: CodonWeights = None,
        context_l: str = '',
        context_r: str = '',
        seed: Optional[Seedable] = None,
        rng: Optional[random.Random] = None
    ) -> None:
        """
        Constructor for the CodingSpace class.

        Parameters
        ----------
        aa_seq
            The amino acid sequence.
        codon_restrictions
            Any codon restrictions in the format e.g. {4: 'TCC'} or {5: ['AGT', 'AGC']}.
        translation_table
            The translation table to use. Leave blank to use standard table.
        codon_weights
            The codon weights to use. Leave blank to sample uniformly.
        context_l
            The context sequence to the left of the coding sequence.
        context_r
            The context sequence to the right of the coding sequence.
        seed
            Seed used to initialise the view's random number generator.
        rng
            Random number generator used by the view for sampling.
        """
        self.forbidden_motifs = forbidden_motifs
        self.max_homopolymer = max_homopolymer

        graph = CodonGraph(
            aa_seq,
            codon_restrictions=codon_restrictions,
            translation_table=translation_table,
            weights=codon_weights,
            context_l=context_l,
            context_r=context_r,
        )

        view = graph.view(
            seed=seed,
            rng=rng,
        )

        self.view = view
        self._update_forbidden_sequences()

    @classmethod
    def from_view(cls, view) -> 'CodingSpace':
        """
        Build a CodingSpace wrapper around an existing graph view.
        """
        obj = cls.__new__(cls)
        obj.view = view
        obj.forbidden_motifs = []
        obj.max_homopolymer = None
        return obj

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

    def __iter__(self) -> Generator[str, None, None]:
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
                *format_banned_sequences(
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

    def sample(self) -> str:
        """
        Sample a DNA sequence from this coding space.

        Returns
        -------
        A sampled sequence. By default, only the coding sequence is returned.
        """
        return self.view.sample()

    def pin_codons(self, pinned_codons: Dict[int, str]):
        """
        Pin temporary codons in this coding space.

        Parameters
        ----------
        pinned_codons
            A dict specifying which codons to pin, by position.
        """
        self.view.pin_codons(pinned_codons)

    def unpin_codons(self, positions: Sequence[int]):
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

    def clear_pins(self):
        """
        Remove all temporary codon pins from this coding space.
        """
        self.view.clear_pins()

    def set_forbidden_motifs(self, forbidden_motifs: ForbiddenMotifs) -> None:
        """
        Set the forbidden motifs for this coding space.

        Parameters
        ----------
        forbidden_motifs
            Motifs that should be forbidden in generated sequences.
        """
        self.forbidden_motifs = forbidden_motifs
        self._update_forbidden_sequences()

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
        self._update_forbidden_sequences()

    def clear_max_homopolymer(self) -> None:
        """
        Remove the maximum homopolymer constraint from this coding space.
        """
        self.set_max_homopolymer(None)

    def _update_forbidden_sequences(self) -> None:
        """
        Rebuild concrete forbidden sequences and apply them to the view.
        """
        forbidden_sequences = expand_and_validate_sequence_constraints(
            forbidden_motifs=self.forbidden_motifs,
            max_homopolymer=self.max_homopolymer,
            rna=self.translation_table.rna,
        )
        self.view.set_banned_sequences(forbidden_sequences)

    def contains(self, seq: str) -> bool:
        """
        Check whether a DNA sequence is contained in this coding space.

        Parameters
        ----------
        seq
            The sequence to check.

        Returns
        -------
        True if and only if the sequence is contained in this coding space.
        """
        return self.view.contains(seq)

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
    def codon_restrictions(self):
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
    def pinned_codons(self):
        """
        Temporary codon pins currently applied to this coding space.
        """
        return self.view.pinned_codons

    @property
    def forbidden_sequences(self):
        """
        Concrete forbidden nucleotide sequences currently applied to this coding space.
        """
        return tuple(self.view.banned_sequences)

    def enumerate(self) -> Generator[str, None, None]:
        """
        Generate all sequences in this space. If there are many (and often there are
        astronomically many), one would not expect to reach the 'end'. However for smaller
        sequence spaces, such as mutation spaces, it's quite possible to get there.

        Yields
        ------
        str
            A valid DNA sequence.
        """
        yield from self.view.enumerate()

    def mutants(
        self,
        seq: str,
        free_positions: Sequence[int] = None,
    ) -> 'MutationSpace':
        """
        Return a space of mutants relative to a given coding sequence, i.e. a space derived
        from this one but which fixes the sequence on all but the specified positions.

        Parameters
        ----------
        seq
            The sequence to mutate.
        free_positions
            The positions that are allowed to vary.
        """
        cds = seq.upper()

        if not self.contains(cds):
            raise ValueError('CDS is not contained in this coding space.')

        from codeine.space.mutation import MutationSpace

        return MutationSpace(
            space=self,
            cds=cds,
            free_positions=free_positions,
        )
