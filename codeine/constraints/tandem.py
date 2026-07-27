from codeine.constraints.subpaths import SubPath, SubPathConstraint


class TandemRepeatConstraint(SubPathConstraint):
    """
    Forbid exact tandem repeats of a specified repeat-unit length.
    """

    def __init__(self, repeat_length: int, min_copies: int = 2):
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
            raise ValueError('repeat_length must be at least 1')

        if not isinstance(min_copies, int) or isinstance(min_copies, bool):
            raise TypeError('min_copies must be an integer')

        if min_copies < 2:
            raise ValueError('min_copies must be at least 2')

        self.repeat_length = repeat_length
        self.min_copies = min_copies

    def _find_paths(self):
        """
        Find tandem repeat paths that can possibly appear in this graph.

        Returns
        -------
        Any paths that can appear but which we wish to ban.
        """
        repeat_span = self.repeat_length * self.min_copies
        full_sequence_length = len(self.context_l) + 3 * len(self.aa_seq) + len(self.context_r)

        paths = []

        for start in range(full_sequence_length - repeat_span + 1):
            paths.extend(self._find_paths_at(start))

        return tuple(set(paths))

    def _nucleotide_to_pos(self, nt_ix):
        """
        Get the graph index and offset corresponding to a given nucleotide ix.

        Parameters
        ----------
        nt_ix
            The nucleotide ix relative, starting at the beginning of the left context sequence

        Returns
        -------
        Graph position, offset.
        """
        if nt_ix < len(self.context_l):
            return 0, nt_ix

        coding_end = len(self.context_l) + 3 * len(self.aa_seq)

        if nt_ix < coding_end:
            nt_ix -= len(self.context_l)
            return 1 + nt_ix // 3, nt_ix % 3

        return len(self.aa_seq) + 1, nt_ix - coding_end

    def _find_paths_at(self, start):
        """
        Find tandem repeat paths beginning at one nucleotide position.
        """
        repeat_span = self.repeat_length * self.min_copies

        start_pos, start_offset = self._nucleotide_to_pos(start)
        end_pos, _ = self._nucleotide_to_pos(start + repeat_span - 1)

        paths = [('', ())]

        for pos in range(start_pos, end_pos + 1):
            next_paths = []

            for sequence, steps in paths:
                for choice in self.graph.nodes[pos].transitions:
                    next_sequence = sequence + choice

                    # Extract as much of the candidate repeat as has been assembled so far.
                    # The upper part of this list slice may overhang the sequence length and
                    # that's fine, ideal, even.
                    partial_repeat = next_sequence[start_offset:start_offset + repeat_span]

                    # Comparing the sequence with itself shifted by one repeat unit checks every
                    # subsequent copy, including partial copies.
                    #
                    # For example for a repeat length of 2:
                    #
                    #    at AT AT AT       <-- partial_repeat[self.repeat_length:]
                    #       AT AT AT at    <-- partial_repeat[:-self.repeat_length]
                    #
                    if partial_repeat[self.repeat_length:] != partial_repeat[:-self.repeat_length]:
                        continue

                    next_paths.append((next_sequence, steps + ((pos, choice),)))

            paths = next_paths

            if not paths:
                return []

        return [
            SubPath(
                sequence=sequence[start_offset:start_offset + repeat_span],
                steps=steps,
                offset=start_offset,
            )
            for sequence, steps in paths
        ]