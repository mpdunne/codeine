from codeine.constraints.motifs import ForbiddenMotifs


class MaxHomopolymer(ForbiddenMotifs):
    """
    Exclude homopolymers longer than a specified maximum length.
    """

    def __init__(
        self,
        max_length: int,
        allow_single_interruption: bool = True,
    ) -> None:
        """
        Parameters
        ----------
        max_length
            Maximum permitted number of identical nucleotides.
        allow_single_interruption
            If True, a homopolymer may contain one interrupting nucleotide
            different from the repeated nucleotide. The interrupting
            nucleotide does not count toward the homopolymer length.
        """
        if not isinstance(max_length, int):
            raise TypeError('max_length must be an integer.')

        if max_length < 1:
            raise ValueError('max_length must be at least 1.')

        if not isinstance(allow_single_interruption, bool):
            raise TypeError('allow_single_interruption must be a boolean.')

        self.max_length = max_length
        self.allow_single_interruption = allow_single_interruption

        banned_homopolymers = [nt * (max_length + 1) for nt in 'ACGT']

        if allow_single_interruption:
            for nt in 'ACGT':
                for interruption in 'ACGT':
                    if interruption == nt:
                        continue

                    for split in range(1, max_length + 1):
                        banned_homopolymers.append(
                            nt * split
                            + interruption
                            + nt * (max_length + 1 - split)
                        )

        super().__init__(banned_homopolymers)
