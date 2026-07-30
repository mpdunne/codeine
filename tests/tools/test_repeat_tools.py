import random
import pytest

from itertools import product

from codeine.tools.repeats import contains_direct_repeat, contains_hairpin, contains_inverted_repeat


################
# Helpers!
################


def reverse_complement(sequence):
    return sequence.translate(str.maketrans('ACGT', 'TGCA'))[::-1]


def brute_contains_direct_repeat(sequence, repeat_length, min_distance=0, max_distance=None):
    for start_l in range(len(sequence) - repeat_length + 1):
        repeat = sequence[start_l:start_l + repeat_length]

        min_start_r = start_l + repeat_length + min_distance
        max_start_r = len(sequence) - repeat_length

        if max_distance is not None:
            max_start_r = min(max_start_r, start_l + repeat_length + max_distance)

        for start_r in range(min_start_r, max_start_r + 1):
            if sequence[start_r:start_r + repeat_length] == repeat:
                return True

    return False


def brute_contains_inverted_repeat(sequence, repeat_length, min_distance=0, max_distance=None):
    for start_l in range(len(sequence) - repeat_length + 1):
        repeat = reverse_complement(sequence[start_l:start_l + repeat_length])

        min_start_r = start_l + repeat_length + min_distance
        max_start_r = len(sequence) - repeat_length

        if max_distance is not None:
            max_start_r = min(max_start_r, start_l + repeat_length + max_distance)

        for start_r in range(min_start_r, max_start_r + 1):
            if sequence[start_r:start_r + repeat_length] == repeat:
                return True

    return False


################
# Direct repeats
################


@pytest.mark.parametrize(
    'sequence,repeat_length',
    [
        ('ATGATG', 3),
        ('ATCGACCAATTGATCGACCAATTG', 12),
        ('GTCAGGATCCGATGCAAT' * 2, 18),
        ('A' * 2000 + 'CG' * 500 + 'A' * 2000, 1000),
    ],
)
def test_contains_direct_repeat(sequence, repeat_length):
    assert contains_direct_repeat(sequence, repeat_length)


@pytest.mark.parametrize(
    'sequence,repeat_length',
    [
        ('ATGATA', 3),
        ('ATCGACCAATTGATCGACCAATTA', 12),
        ('ACGTGCACTGATCAGTACGATCGT', 12),
    ],
)
def test_does_not_contain_direct_repeat(sequence, repeat_length):
    assert not contains_direct_repeat(sequence, repeat_length)


def test_direct_repeat_distance():
    sequence = 'ATG' + 'ACGTAC' + 'ATG'

    assert not contains_direct_repeat(sequence, 3, min_distance=7)
    assert contains_direct_repeat(sequence, 3, min_distance=6)
    assert contains_direct_repeat(sequence, 3, max_distance=6)
    assert not contains_direct_repeat(sequence, 3, max_distance=5)


##################
# Inverted repeats
##################


@pytest.mark.parametrize(
    'repeat',
    [
        'ATG',
        'ATCGACCAATTG',
        'ATGATCAAAGAATAC',
        'GTCAGGATCCGATGCAAT' * 50,
    ],
)
def test_contains_inverted_repeat(repeat):
    sequence = repeat + 'AGAGAG' + reverse_complement(repeat)

    assert contains_inverted_repeat(sequence, len(repeat))


@pytest.mark.parametrize(
    'repeat',
    [
        'ATG',
        'ATCGACCAATTG',
        'ATGATCAAAGAATAC',
    ],
)
def test_does_not_contain_modified_inverted_repeat(repeat):
    inverted = reverse_complement(repeat)
    replacement = 'A' if inverted[0] != 'A' else 'C'
    sequence = repeat + 'AGAGAG' + replacement + inverted[1:]

    assert not contains_inverted_repeat(sequence, len(repeat))


def test_inverted_repeat_distance():
    repeat = 'ATCGACCAATTG'
    sequence = repeat + 'ACGTAC' + reverse_complement(repeat)

    assert not contains_inverted_repeat(sequence, len(repeat), min_distance=7)
    assert contains_inverted_repeat(sequence, len(repeat), min_distance=6)
    assert contains_inverted_repeat(sequence, len(repeat), max_distance=6)
    assert not contains_inverted_repeat(sequence, len(repeat), max_distance=5)


##########
# Hairpins
##########


def test_contains_hairpin():
    stem = 'ATCGACCAATTG'
    sequence = stem + 'AAAA' + reverse_complement(stem)

    assert contains_hairpin(sequence, stem_length=len(stem), min_loop_length=4)


def test_hairpin_loop_length():
    stem = 'ATCGACCAATTG'
    sequence = stem + 'AAAAAA' + reverse_complement(stem)

    assert not contains_hairpin(sequence, len(stem), min_loop_length=7)
    assert contains_hairpin(sequence, len(stem), min_loop_length=6)
    assert contains_hairpin(sequence, len(stem), max_loop_length=6)
    assert not contains_hairpin(sequence, len(stem), max_loop_length=5)


def test_inverted_repeat_in_right_half():
    # An early implementation missed this case, so test for it.
    sequence = 'AAAAATGAACGCAT'

    assert contains_inverted_repeat(sequence, 3)


##########################
# Bitmask implementation
##########################


@pytest.mark.parametrize(
    'repeat_length,min_distance,max_distance',
    [
        (3, 0, None),
        (6, 0, None),
        (6, 3, 12),
        (12, 0, 30),
    ],
)
def test_direct_repeat_matches_brute_force(
    repeat_length,
    min_distance,
    max_distance,
):
    for seed in range(20):
        rng = random.Random(seed)
        sequence = ''.join(rng.choice('ACGT') for _ in range(120))

        assert contains_direct_repeat(sequence, repeat_length, min_distance, max_distance) == \
               brute_contains_direct_repeat(sequence, repeat_length, min_distance, max_distance)


@pytest.mark.parametrize(
    'repeat_length,min_distance,max_distance',
    [
        (3, 0, None),
        (6, 0, None),
        (6, 3, 12),
        (12, 0, 30),
    ],
)
def test_inverted_repeat_matches_brute_force(
    repeat_length,
    min_distance,
    max_distance,
):
    for seed in range(20):
        rng = random.Random(seed)
        sequence = ''.join(rng.choice('ACGT') for _ in range(120))

        assert contains_inverted_repeat(sequence, repeat_length, min_distance, max_distance) == \
               brute_contains_inverted_repeat(sequence, repeat_length, min_distance, max_distance)


@pytest.mark.parametrize(
    'repeat_length,min_distance,max_distance',
    [
        (2, 0, None),
        (2, 1, 2),
        (3, 0, None),
    ],
)
def test_direct_repeat_matches_brute_force_exhaustive(repeat_length, min_distance, max_distance):
    for nts in product('ACGT', repeat=7):
        sequence = ''.join(nts)

        assert contains_direct_repeat(sequence, repeat_length, min_distance, max_distance) == \
               brute_contains_direct_repeat(sequence, repeat_length, min_distance, max_distance)


@pytest.mark.parametrize(
    'repeat_length,min_distance,max_distance',
    [
        (2, 0, None),
        (2, 1, 2),
        (3, 0, None),
    ],
)
def test_inverted_repeat_matches_brute_force_exhaustive(repeat_length, min_distance, max_distance):
    for nts in product('ACGT', repeat=7):
        sequence = ''.join(nts)

        assert contains_inverted_repeat(sequence, repeat_length, min_distance, max_distance) == \
               brute_contains_inverted_repeat(sequence, repeat_length, min_distance, max_distance)


#######
# RNA
#######


def test_direct_repeat_accepts_rna():
    assert contains_direct_repeat('AUGCCCAUG', 3)


def test_inverted_repeat_accepts_rna():
    assert contains_inverted_repeat('AUGCCCAU', 3)
