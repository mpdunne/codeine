from codeine.constraints.motifs import ForbiddenMotifConstraint


class HomopolymerConstraint(ForbiddenMotifConstraint):
    """
    Exclude homopolymers longer than a specified maximum length.
    """

    def __init__(self, max_length: int) -> None:
        """
        Parameters
        ----------
        max_length
            Maximum permitted number of consecutive identical nucleotides.
        """
        if not isinstance(max_length, int):
            raise TypeError('max_length must be an integer.')

        if max_length < 1:
            raise ValueError('max_length must be at least 1.')

        self.max_length = max_length

        banned_homopolymers = [nt * (max_length + 1) for nt in 'ACGT']
        super().__init__(banned_homopolymers)
