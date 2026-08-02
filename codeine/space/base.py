import pickle

from pathlib import Path
from typing import Dict, Iterator, List, Optional, TYPE_CHECKING, Type, TypeVar, Union

if TYPE_CHECKING:
    from codeine.graph.view import CodonGraphView

from codeine.graph.base import CodonRestriction
from codeine.translation.tables import TranslationTable
from codeine.translation.weights import CodonWeights


T = TypeVar('T', bound='Space')


class Space:
    """
    Base class for spaces of valid coding sequences.
    """

    view: 'CodonGraphView'

    @classmethod
    def load(cls: Type[T], path) -> T:
        """
        Load a space from disc.
        """
        with Path(path).open('rb') as f:
            return pickle.load(f)

    def save(self, path) -> None:
        """
        Save this space to disc.
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
        Iterate over all valid sequences in this space.
        Be aware that "all valid sequences" can be astronomically many!

        Yields
        ----------
        All valid sequences in the space, in order.
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

    def __len__(self) -> int:
        """
        This method exists only to provide a helpful error message. Being able to call
        len(space) is a totally reasonable thing to expect, but python's ``len`` hits a limit
        for very large spaces.

        Use ``count()`` or ``n_valid_sequences`` instead.
        """
        raise TypeError(
            f'len() is not supported for {type(self).__name__}; '
            f'use .count() or .n_valid_sequences instead.'
        )

    def compile(self) -> None:
        """
        Compile this space.

        If the space is still uncompiled at the point of counting, sampling, or enumeration, ``compile()``
        will be called automatically. Explicit compilation is useful for benchmarking or for preparing
        a space before repeated sampling.
        """
        self.view.compile()

    def sample(self, n: Optional[int] = None) -> Union[str, List[str]]:
        """
        Sample one or more sequences from this space.

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
        Check whether a coding sequence is contained in this space.

        Parameters
        ----------
        seq
            The sequence to check.

        Returns
        -------
        True if and only if the sequence is contained in this space.
        """
        return self.view.contains(seq)

    def count(self) -> int:
        """
        Return the number of valid sequences in this space.

        Returns
        -------
        int
            The number of valid sequences.
        """
        return self.n_valid_sequences

    def set_codon_weights(self, codon_weights: CodonWeights) -> None:
        """
        Set the codon weights used when sampling from this space.
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
        The amino acid sequence for this space.
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
    def fixed_codons(self) -> Dict[int, CodonRestriction]:
        """
        The fixed codon restrictions in this space.
        """
        return self.view.fixed_codons

    @property
    def context_l(self) -> str:
        """
        The left context sequence.
        """
        return self.view.context_l

    @property
    def context_r(self) -> str:
        """
        The right context sequence.
        """
        return self.view.context_r

    @property
    def pinned_codons(self) -> Dict[int, Union[str, List[str]]]:
        """
        Temporary codon pins currently applied to this space.
        """
        return self.view.pinned_codons
