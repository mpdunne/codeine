from codeine.constraints.repeats import InvertedRepeats


class Hairpins(InvertedRepeats):
    """
    Forbid nucleotide sequences capable of forming specified hairpins.
    """

    def __init__(
        self,
        stem_length: int,
        min_loop_length: int = 3,
        max_loop_length: int = None,
    ):
        """
        Parameters
        ----------
        stem_length
            Length of each complementary stem in nucleotides.
        min_loop_length
            Minimum number of nucleotides in the hairpin loop.
        max_loop_length
            Maximum number of nucleotides in the hairpin loop, or ``None`` for
            no maximum.
        """
        super().__init__(
            repeat_length=stem_length,
            min_distance=min_loop_length,
            max_distance=max_loop_length,
        )

        self.stem_length = stem_length
        self.min_loop_length = min_loop_length
        self.max_loop_length = max_loop_length
