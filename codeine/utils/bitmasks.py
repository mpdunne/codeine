NT_TO_MASK = {
    'A': 1,
    'C': 2,
    'G': 4,
    'T': 8,
    'U': 8,
}


def sequence_to_nt_bitmasks(sequence):
    """
    Encode a nucleotide sequence as one bitmask per position.
    """
    return [NT_TO_MASK[nt] for nt in sequence]


def choices_to_nt_bitmasks(choices):
    """
    Get the possible nucleotides at each sequence position, encoded as a bitmask.

    Parameters
    ----------
    choices
        Sequence of choice sets, e.g. context sequence or codons at each position.

    Returns
    -------
    Bitmask of possible nucleotides at each nucleotide position.
    """
    nt_masks = []

    for step_choices in choices:
        for choice_offset in range(len(step_choices[0])):

            if len(set([len(choice) for choice in step_choices])) != 1:
                raise ValueError('Choices within each choice set must have equal length')

            mask = 0

            for choice in step_choices:
                mask |= NT_TO_MASK[choice[choice_offset]]

            nt_masks.append(mask)

    return nt_masks


def pack_nt_bitmasks(nt_bitmasks):
    """
    Pack 4-bit nucleotide masks into consecutive nibbles of one integer.
    """
    packed = bytearray((len(nt_bitmasks) + 1) // 2)

    for ix, mask in enumerate(nt_bitmasks):
        packed[ix // 2] |= mask << (4 * (ix % 2))

    return int.from_bytes(packed, 'little')
