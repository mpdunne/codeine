from typing import Collection, Generator, Optional, Set

from codeine.sequence.space import CodingSpace


class MutationSpace:
    """
    Represents the set of valid coding sequences that can be reached by
    mutating a reference CDS in a given CodingSpace, subject to constraints.

    A MutationSpace is defined by:

    - A CodingSpace containing the global sequence constraints.
    - A reference CDS.
    - A set of codon positions that are free to mutate.

    Positions that are not free are temporarily considered frozen and
    will remain identical to the reference CDS.
    """

    def __init__(self,
                 space: CodingSpace,
                 cds: str,
                 free_positions: Optional[Collection[int]] = None,
                 ):
        """
        Constructor for the MutationSpace class.

        Parameters
        ----------
        space
            The CodingSpace object to which this given sequence should belong.
        cds
            The parent/reference CDS.
        free_positions
            Which positions are allowed to change?
        """
        self.space = space
        self.cds = cds

        if free_positions is None:
            free_positions = range(1, len(space.view.aa_seq) + 1)

        self._free_positions: Set[int] = set()
        self.set_free_positions(free_positions)

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
        return self.space[index]

    def __iter__(self) -> Generator[str, None, None]:
        """
        Iterate over all valid sequences in this mutation space.
        Be aware that "all valid sequences" can be astronomically many!

        Yields
        ----------
        All valid sequences in the coding space, in order.
        """
        yield from self.space

    def __contains__(self, seq: str) -> bool:
        """
        Does the given seq exist in this space?

        Returns
        ----------
        True if and only if this is a valid sequence in this space.
        """
        return seq in self.space

    @property
    def free_positions(self):
        """
        Codon positions that are currently free to mutate.
        """
        return frozenset(self._free_positions)

    @property
    def frozen_positions(self):
        """
        Codon positions that are currently fixed to the reference CDS.
        """
        all_positions = set(range(1, len(self.space.view.aa_seq) + 1))
        return frozenset(all_positions - self._free_positions)

    def _validate_positions(self, positions: Collection[int]) -> Set[int]:
        """
        Validate a collection of codon positions.
        """
        positions = set(positions)
        invalid = [pos for pos in positions if pos < 1 or pos > len(self.space.view.aa_seq)]
        if invalid:
            raise ValueError(f'Invalid codon positions: {sorted(invalid)}')

        return positions

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
        Update the pins on the underlying space.
        """
        pins = {pos: self._codon_at_position(pos) for pos in self.frozen_positions}
        self.space.set_pinned_codons(pins)

    def contains(self, seq: str) -> bool:
        """
        Check whether a DNA sequence is contained in this mutation space.

        Parameters
        ----------
        seq
            The sequence to check

        Returns
        -------
        True if and only if the sequence is contained in this mutation space.
        """
        return self.space.contains(seq)

    @property
    def n_valid_variants(self) -> int:
        """
        Number of valid variants under the current mutation constraints.
        """
        return self.space.n_valid_sequences

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
        self._free_positions = set(range(1, len(self.space.view.aa_seq) + 1))
        self._update_pins()

    def sample(self) -> str:
        """
        Sample a variant from this mutation space.

        Returns
        -------
        A sampled string sequence from this mutation space.
        """
        return self.space.sample()

    def enumerate(self) -> Generator[str, None, None]:
        """
        Generate all sequences in this mutation space.

        Yields
        ------
        str
            A valid DNA sequence.
        """
        yield from self.space.enumerate()
