from typing import Dict, Generator, Optional, Sequence

from codeine.sequence.graph import CodonGraph
from codeine.translation.tables import TranslationTable
from codeine.translation.weights import CodonWeights


class CodingSpace:
    """
    Class representing coding sequence space, for sampling and mutating CDS coding sequences.
    """

    def __init__(
            self,
            aa_seq: str,
            codon_restrictions: Optional[Dict[int, str]] = None,
            translation_table: TranslationTable = None,
            codon_weights: CodonWeights = None,
            context_l: str = '',
            context_r: str = '',
    ) -> None:
        """
        Constructor for the CodingSpace class.

        Parameters
        ----------
        aa_seq
            The amino acid sequence.
        codon_restrictions
            Any codon restrictions in the format e.g. {4: 'TCC'} or {5: ['AGT', 'AGC']}
        translation_table
            The translation table to use. Leave blank to use standard table.
        codon_weights
            The codon weights to use. Leave blank to sample uniformly.
        context_l
            The context sequence to the left of the coding sequence
        context_r
            The context sequence to the right of the coding sequence
        """
        self.view = CodonGraph(
            aa_seq,
            codon_restrictions=codon_restrictions,
            translation_table=translation_table,
            weights=codon_weights,
            context_l=context_l,
            context_r=context_r,
        ).view()

    @classmethod
    def from_graph(cls, graph: CodonGraph) -> 'CodingSpace':
        return cls.from_view(graph.view())

    @classmethod
    def from_view(cls, view) -> 'CodingSpace':
        obj = cls.__new__(cls)
        obj.view = view
        return obj

    def __getitem__(self, index: int) -> str:
        """
        Return the valid sequence at a given index.

        Parameters
        ----------
        index
            Zero-based sequence index.

        Returns
        -------
        str
            The indexed valid DNA sequence.
        """
        return self.view[index]

    def __len__(self):
        """
        The number of valid sequences in this graph.

        Returns
        -------
        The number of valid sequences in this graph.
        """
        return self.n_valid_sequences

    def __iter__(self) -> Generator[str, None, None]:
        """
        Iterate over all valid sequences in this coding space.
        Be aware that "all valid sequences" can be astronomically many!

        Yields
        ----------
        All valid sequences in the coding space, in order.
        """
        yield from self.view

    def sample(self) -> str:
        """
        Sample a DNA sequence from this coding space.

        Returns
        -------
        A sampled sequence. By default, only the coding sequence is returned.
        """
        return self.view.sample()

    def pin_codons(self, pinned_codons):
        """
        Pin (temporarily fix) a codon in the codon graph.

        Parameters
        ----------
        pinned_codons:
            A dict specifying which codons to pin, by pos: codon
        """
        self.view.pin_codons(pinned_codons)

    def unpin_codons(self, positions):
        """
        Unpin codon nodes by pos.

        Parameters
        ----------
        positions:
            A list of positions
        """
        self.view.unpin_codons(positions)

    def clear_pins(self):
        """
        Remove all codon pins from the generator.
        """
        self.view.clear_pins()

    def contains(self, seq: str) -> bool:
        """
        Check whether a DNA sequence is contained in this coding space.

        Parameters
        ----------
        seq
            The sequence to check

        Returns
        -------
        True if and only if the sequence is contained in this coding space.
        """
        return self.view.contains(seq)

    @property
    def n_valid_sequences(self) -> int:
        """
        The number of valid sequences in this space.

        Returns
        -------
        The number of valid sequences in this space.
        """
        return self.view.n_valid_sequences

    @property
    def translation_table(self) -> TranslationTable:
        """
        The translation table being used in this space.

        Returns
        -------
        The TranslationTable being used.
        """
        return self.view.graph.tt

    @property
    def codon_weights(self) -> CodonWeights:
        """
        The codon weights being used in this space.

        Returns
        -------
        The CodonWeights being used.
        """
        return self.view.graph.cw

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

    def mutants(self,
                seq: str,
                positions: Sequence[int],
                ) -> 'CodingSpace':
        """
        Return a space of mutants relative to a given coding sequence, i.e. a space derived
        from this one but which fixes the sequence on all but the specified positions.

        Parameters
        ----------
        seq
            The sequence to mutate.
        positions
            The positions that are allowed to vary.
        """
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
