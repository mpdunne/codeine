from codeine.constraints.subpaths import SubPathConstraint


class TandemRepeatConstraint(SubPathConstraint):
    """
    Forbid exact tandem repeats of a specified repeat-unit length.
    """

    def __init__(
        self,
        repeat_length: int,
        min_copies: int = 2,
    ):
        """
        Parameters
        ----------
        repeat_length
            Length of the repeated unit in nucleotides.
        min_copies
            Minimum number of consecutive copies to forbid.
        """
        super().__init__()

        if not isinstance(repeat_length, int) or isinstance(repeat_length, bool):
            raise TypeError('repeat_length must be an integer')

        if repeat_length < 1:
            raise ValueError('repeat_length must be at least 2')

        if not isinstance(min_copies, int) or isinstance(min_copies, bool):
            raise TypeError('min_copies must be an integer')

        if min_copies < 2:
            raise ValueError('min_copies must be at least 2')

        self.repeat_length = repeat_length
        self.min_copies = min_copies

    def _find_paths(self):
        """
        Find tandem repeats sequences that can possibly appear in this graph.

        Returns
        -------
        Any paths that can appear but which we wish to ban.
        """
        return ()
