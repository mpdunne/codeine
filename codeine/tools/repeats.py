from bisect import bisect_right


_NT_TRANS = str.maketrans('ACGT', 'TGCA')


def _reverse_complement(sequence):
    return sequence.translate(_NT_TRANS)[::-1]


def _contains_repeat(sequence, repeat_length, min_distance, max_distance, inverted):
    sequence = sequence.upper().replace('U', 'T')
    positions = {}

    for start in range(len(sequence) - repeat_length + 1):
        repeat = sequence[start:start + repeat_length]
        previous = positions.get(_reverse_complement(repeat) if inverted else repeat)

        if previous:
            latest_start = start - repeat_length - min_distance
            previous_ix = bisect_right(previous, latest_start) - 1

            if previous_ix >= 0:
                distance = start - previous[previous_ix] - repeat_length

                if max_distance is None or distance <= max_distance:
                    return True

        positions.setdefault(repeat, []).append(start)

    return False


def contains_direct_repeat(sequence, repeat_length, min_distance=0, max_distance=None):
    """
    Determine whether a sequence contains an exact direct repeat.

    Distance is the number of nucleotides between the two repeated sequences.
    """
    if len(sequence) < 2 * repeat_length + min_distance:
        return False

    return _contains_repeat(sequence, repeat_length, min_distance, max_distance, inverted=False)


def contains_inverted_repeat(sequence, repeat_length, min_distance=0, max_distance=None):
    """
    Determine whether a sequence contains an exact inverted repeat.

    Distance is the number of nucleotides between the two reverse-complementary
    sequences.
    """
    if len(sequence) < 2 * repeat_length + min_distance:
        return False

    return _contains_repeat(sequence, repeat_length, min_distance, max_distance, inverted=True)


def contains_hairpin(sequence, stem_length, min_loop_length=0, max_loop_length=None):
    """
    Return whether a sequence contains an exact hairpin.

    ``stem_length`` is the length of each reverse-complementary stem and loop
    length is the number of nucleotides between them.
    """
    return contains_inverted_repeat(
        sequence,
        repeat_length=stem_length,
        min_distance=min_loop_length,
        max_distance=max_loop_length,
    )
