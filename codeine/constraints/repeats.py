from abc import ABC
from typing import Optional

from codeine.constraints.base import Constraint
from codeine.graph.base import CodonGraph
from codeine.utils.bitmasks import choices_to_nt_bitmasks, pack_nt_bitmasks


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

        self.graph: Optional[CodonGraph] = None
        self.aa_seq = None
        self.context_l = None
        self.context_r = None
        self.translation_table = None

        self.nt_trans = None
        self.full_sequence_length = None

        self.reference_choices = None
        self.compare_choices = None

        self.nt_bitmasks_reference = None
        self.nt_bitmasks_compare = None

        self.nt_bitmasks_reference_packed = None
        self.nt_bitmasks_compare_packed = None
        self.nt_position_mask = None

        self.nt_positions_reference = None
        self.nt_positions_compare = None

        self.repeats = None

    def initial_state(self):
        raise NotImplementedError

    def advance(self, state, pos, choice):
        raise NotImplementedError

    def link(self, graph):
        self.graph = graph
        self.aa_seq = graph.aa_seq
        self.context_l = graph.context_l
        self.context_r = graph.context_r
        self.translation_table = graph.tt

        self.nt_trans = str.maketrans('ACGU', 'UGCA') if graph.tt.rna else str.maketrans('ACGT', 'TGCA')

        self.reference_choices = [
            tuple(node.transitions)
            for node in self.graph.nodes
            if node is not self.graph.end_node
        ]

        self.compare_choices = (
            [
                tuple(self._reverse_complement(choice) for choice in choice_set)
                for choice_set in self.reference_choices[::-1]
            ]
            if self.inverted
            else self.reference_choices
        )

        self.nt_bitmasks_reference = choices_to_nt_bitmasks(self.reference_choices)
        self.nt_bitmasks_compare = choices_to_nt_bitmasks(self.compare_choices)

        self.nt_positions_reference = self._get_nt_positions(self.reference_choices)
        self.nt_positions_compare = self._get_nt_positions(self.compare_choices)

        self.full_sequence_length = len(self.nt_bitmasks_reference)

        # Pack each 4-bit nucleotide mask into one nibble of a single integer.
        self.nt_bitmasks_reference_packed = pack_nt_bitmasks(self.nt_bitmasks_reference)
        self.nt_bitmasks_compare_packed = pack_nt_bitmasks(self.nt_bitmasks_compare)

        # 0001 0001 0001 ...: one marker bit per nucleotide position.
        self.nt_position_mask = ((1 << (4 * self.full_sequence_length)) - 1) // 15

        self.repeats = []

        min_shift = -(self.full_sequence_length - self.repeat_length) if self.inverted \
            else self.repeat_length + self.min_distance
        max_shift = self.full_sequence_length - self.repeat_length

        # Cap the max shift in the direct case. For inverted repeats, distance
        # depends on where in the alignment the repeat occurs and is checked later.
        if not self.inverted and self.max_distance is not None:
            max_shift = min(max_shift, self.repeat_length + self.max_distance)

        for shift in range(min_shift, max_shift + 1):
            self._find_repeats_at_shift(shift)

    @property
    def is_trivial(self) -> bool:
        return not self.repeats

    def _reverse_complement(self, sequence: str):
        """
        Return the reverse complement of a nucleotide sequence.
        """
        return sequence.translate(self.nt_trans)[::-1]

    @staticmethod
    def _get_nt_positions(choices):
        """
        Get the graph position and choice offset at each nucleotide position.
        """
        return [(pos, choice_offset) for pos, choice_set in enumerate(choices)
                for choice_offset in range(len(choice_set[0]))]

    @staticmethod
    def _get_run_starts(matches, length):
        """
        Find starts of runs containing at least ``length`` matching positions.

        The input has one possible marker bit every four bits, corresponding to
        nucleotide positions. The returned integer contains a marker bit at each
        position where a sufficiently long run begins.
        """
        blocks = {1: matches}
        size = 1

        # Build runs of lengths 1, 2, 4, 8, ...
        while size * 2 <= length:
            blocks[size * 2] = blocks[size] & (blocks[size] >> (4 * size))
            size *= 2

        starts = blocks[size]
        covered = size
        remaining = length - size

        # Combine powers of two to make the requested run length.
        while remaining:
            size = 1 << (remaining.bit_length() - 1)
            starts &= blocks[size] >> (4 * covered)
            covered += size
            remaining -= size

        return starts

    def _find_repeats_at_shift(self, shift):
        """
        Find candidate repeats at one alignment of the reference and compare sequences.
        """

        # Compare all nucleotide positions at this shift in a few bigint operations.
        compare = (
            self.nt_bitmasks_compare_packed >> (4 * shift)
            if shift >= 0
            else self.nt_bitmasks_compare_packed << (4 * abs(shift))
        )

        nt_intersections = self.nt_bitmasks_reference_packed & compare

        # Collapse each non-zero 4-bit nucleotide intersection to its marker bit.
        matches = (
            nt_intersections
            | (nt_intersections >> 1)
            | (nt_intersections >> 2)
            | (nt_intersections >> 3)
        ) & self.nt_position_mask

        starts = self._get_run_starts(matches, self.repeat_length)

        # For inverted repeats, later start positions would map back onto or before
        # the reference repeat, so they can never form a valid pair.
        if self.inverted:
            max_start_l = (self.full_sequence_length - shift - 2 * self.repeat_length - self.min_distance) // 2

            if max_start_l < 0:
                return

            starts &= (1 << (4 * (max_start_l + 1))) - 1

        # Only visit actual candidate repeat starts.
        while starts:
            bit = starts & -starts
            start_l = (bit.bit_length() - 1) // 4
            starts ^= bit

            start_compare = shift + start_l

            start_r = self.full_sequence_length - start_compare - self.repeat_length \
                if self.inverted else start_compare

            if start_r <= start_l:
                continue

            distance = start_r - start_l - self.repeat_length

            if distance < self.min_distance:
                continue

            if self.max_distance is not None and distance > self.max_distance:
                continue

            # Check which whole choices can actually participate in the repeat.
            # The nucleotide-level filter above treats positions independently;
            # this step restores dependencies within codons/context choices.
            filtered_choices = self._filter_choices(start_l, start_compare)

            if filtered_choices is None:
                continue

            reference_choices_filtered, compare_choices_filtered = filtered_choices

            requirements = self._get_compatible_choices(
                start_l,
                start_compare,
                reference_choices_filtered,
                compare_choices_filtered,
            )

            self.repeats.append((start_l, start_r, requirements))

    def _get_offsets_by_position_pairs(self, reference_start, compare_start):
        """
        Map the nucleotide comparison onto pairs of graph positions, and determine what pairs of
        offsets within them will be required to have at least one match.

        A  pair of graph positions can have multiple comparison points, because contexts can be
        longer than a single codon.

        Returns
        -------
        aligned_positions
            Dictionary mapping each ``(reference_pos, compare_pos)`` pair to the
            nucleotide offsets within those choices that must be capable of matching
        """
        comparisons = {}

        for offset in range(self.repeat_length):
            reference_pos, reference_offset = self.nt_positions_reference[reference_start + offset]
            compare_pos, compare_offset = self.nt_positions_compare[compare_start + offset]

            comparisons.setdefault((reference_pos, compare_pos), []).append(
                (reference_offset, compare_offset)
            )

        return comparisons

    def _filter_choices(self, reference_start, compare_start):
        """
        Find the codon choices that can realise a candidate repeat.

        Align the ``repeat_length`` nucleotides starting at ``reference_start`` and
        ``compare_start``, and remove codon choices that prohibit a match.

        Filtering continues until no further choices can be removed. By construction, any
        non-empty result guarantees that at least one match exists.

        Returns
        -------
        reference_choices, compare_choices
            Dictionaries mapping each involved position to the set of choices that
            remain possible, or ``None`` if no repeat is possible here.
        """
        comparisons = self._get_offsets_by_position_pairs(reference_start, compare_start)

        reference_choices = {pos: set(self.reference_choices[pos]) for pos, _ in comparisons}
        compare_choices = {pos: set(self.compare_choices[pos]) for _, pos in comparisons}

        # TODO: We can probably use bitmasks for this. We love bitmasks.
        # Remove choices that have no compatible choice at a compared position. Due to offsets,
        # a single position will usually need to be compared against multiple neighbouring positions,
        # Filtering in one of them can affect the choices in another: proceed iteratively until stable.
        while True:
            changed = False

            for (reference_pos, compare_pos), offsets in comparisons.items():
                next_reference_choices = set()
                next_compare_choices = set()

                # Find choices on either side that have at least one compatible partner.
                for reference_choice in reference_choices[reference_pos]:
                    for compare_choice in compare_choices[compare_pos]:
                        if all(
                            reference_choice[reference_offset] == compare_choice[compare_offset]
                            for reference_offset, compare_offset in offsets
                        ):
                            next_reference_choices.add(reference_choice)
                            next_compare_choices.add(compare_choice)

                if not next_reference_choices or not next_compare_choices:
                    return None

                if next_reference_choices != reference_choices[reference_pos]:
                    reference_choices[reference_pos] = next_reference_choices
                    changed = True

                if next_compare_choices != compare_choices[compare_pos]:
                    compare_choices[compare_pos] = next_compare_choices
                    changed = True

            if not changed:
                return reference_choices, compare_choices

    def _get_compatible_choices(
        self,
        reference_start,
        compare_start,
        reference_choices,
        compare_choices,
    ):
        """
        Get the compare choices permitted by each reference choice.

        Returns
        -------
        requirements
            Dictionary mapping each reference position to its compared positions and,
            for each reference choice, the compare choices compatible with it.
        """
        comparisons = self._get_offsets_by_position_pairs(reference_start, compare_start)
        requirements = {}

        for (reference_pos, compare_pos), offsets in comparisons.items():
            allowed_choices = {}

            for reference_choice in reference_choices[reference_pos]:
                allowed_choices[reference_choice] = {
                    compare_choice
                    for compare_choice in compare_choices[compare_pos]
                    if all(
                        reference_choice[reference_offset] == compare_choice[compare_offset]
                        for reference_offset, compare_offset in offsets
                    )
                }

            requirements.setdefault(reference_pos, []).append(
                (compare_pos, allowed_choices)
            )

        return requirements


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
