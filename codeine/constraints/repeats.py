from abc import ABC

from codeine.constraints.base import Constraint


class RepeatConstraint(Constraint, ABC):
    """
    Base class for direct and inverted repeat constraints.
    """

    def __init__(
        self,
        repeat_length: int,
        min_distance: int = 0,
        max_distance: int = None,
        inverted: bool = False,
    ):
        """
        Parameters
        ----------
        repeat_length
            Length of each repeated sequence in nucleotides.
        min_distance
            Minimum number of nucleotides between the repeated sequences.
        max_distance
            Maximum number of nucleotides between the repeated sequences, or
            ``None`` for no maximum.
        inverted
            Whether the second sequence is reverse-complementary to the first.
        """
        super().__init__()

        if not isinstance(repeat_length, int) or isinstance(repeat_length, bool):
            raise TypeError('repeat_length must be an integer')

        if repeat_length < 1:
            raise ValueError('repeat_length must be at least 1')

        if not isinstance(min_distance, int) or isinstance(min_distance, bool):
            raise TypeError('min_distance must be an integer')

        if min_distance < 0:
            raise ValueError('min_distance must be at least 0')

        if max_distance is not None:
            if not isinstance(max_distance, int) or isinstance(max_distance, bool):
                raise TypeError('max_distance must be an integer or None')

            if max_distance < min_distance:
                raise ValueError('max_distance must be at least min_distance')

        if not isinstance(inverted, bool):
            raise TypeError('inverted must be a boolean')

        self.repeat_length = repeat_length
        self.min_distance = min_distance
        self.max_distance = max_distance
        self.inverted = inverted

    def initial_state(self):
        raise NotImplementedError

    def link(self, graph):
        raise NotImplementedError

    def advance(self, state, pos, choice):
        raise NotImplementedError


class DirectRepeatConstraint(RepeatConstraint):
    """
    Forbid identical nucleotide sequences separated by a specified distance.
    """

    def __init__(
        self,
        repeat_length: int,
        min_distance: int = 0,
        max_distance: int = None,
    ):
        """
        Parameters
        ----------
        repeat_length
            Length of each repeated sequence in nucleotides.
        min_distance
            Minimum number of nucleotides between the repeated sequences.
        max_distance
            Maximum number of nucleotides between the repeated sequences, or
            ``None`` for no maximum.
        """
        super().__init__(
            repeat_length=repeat_length,
            min_distance=min_distance,
            max_distance=max_distance,
            inverted=False,
        )


class InvertedRepeatConstraint(RepeatConstraint):
    """
    Forbid reverse-complementary nucleotide sequences separated by a specified
    distance.
    """

    def __init__(
        self,
        repeat_length: int,
        min_distance: int = 0,
        max_distance: int = None,
    ):
        """
        Parameters
        ----------
        repeat_length
            Length of each reverse-complementary sequence in nucleotides.
        min_distance
            Minimum number of nucleotides between the repeated sequences.
        max_distance
            Maximum number of nucleotides between the repeated sequences, or
            ``None`` for no maximum.
        """
        super().__init__(
            repeat_length=repeat_length,
            min_distance=min_distance,
            max_distance=max_distance,
            inverted=True,
        )
