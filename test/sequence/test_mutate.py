import pytest
from unittest.mock import MagicMock

from codeine.sequence.space import CodingSpace
from codeine.sequence.mutate import MutationSpace


def test_mutation_space_defaults_to_all_positions_free():
    space = CodingSpace('MIKEY')
    cds = space[0]
    muts = MutationSpace(space, cds)

    assert muts.free_positions == frozenset(range(1, len(space.view.aa_seq) + 1))
    assert muts.frozen_positions == frozenset()


def test_mutation_space_can_set_free_positions():
    space = CodingSpace('MIKEY')
    cds = space[0]
    muts = MutationSpace(space, cds)

    muts.set_free_positions([1, 3])

    assert muts.free_positions == frozenset({1, 3})
    assert muts.frozen_positions == frozenset({2, 4, 5})


def test_mutation_space_can_freeze_positions():
    space = CodingSpace('MIKEY')
    cds = space[0]
    muts = MutationSpace(space, cds)

    muts.freeze_positions([1])

    assert 1 not in muts.free_positions
    assert 1 in muts.frozen_positions


def test_mutation_space_can_unfreeze_positions():
    space = CodingSpace('MIKEY')
    cds = space[0]
    muts = MutationSpace(space, cds, free_positions=[3])

    muts.unfreeze_positions([1])

    assert muts.free_positions == frozenset({1, 3})


def test_mutation_space_can_freeze_all():
    space = CodingSpace('MIKEY')
    cds = space[0]
    muts = MutationSpace(space, cds)

    muts.freeze_all()

    assert muts.free_positions == frozenset()
    assert muts.frozen_positions == frozenset(range(1, len(space.view.aa_seq) + 1))


def test_mutation_space_can_unfreeze_all():
    space = CodingSpace('MIKEY')
    cds = space[0]
    muts = MutationSpace(space, cds, free_positions=[])

    muts.unfreeze_all()

    assert muts.free_positions == frozenset(range(1, len(space.view.aa_seq) + 1))
    assert muts.frozen_positions == frozenset()


def test_mutation_space_rejects_invalid_positions():
    space = CodingSpace('MIKEY')
    cds = space[0]
    muts = MutationSpace(space, cds)

    with pytest.raises(ValueError):
        muts.set_free_positions([-1])

    with pytest.raises(ValueError):
        muts.set_free_positions([len(space.view.aa_seq) + 1])


def test_mutation_space_updates_pins_when_free_positions_change():
    space = CodingSpace('MIKEY')
    cds = space[0]

    muts = MutationSpace(space, cds)
    muts.set_free_positions([2])

    expected_pins = {1: [cds[0:3]], 3: [cds[6:9]], 4: [cds[9:12]], 5: [cds[12:15]]}
    assert space.view.pinned_codons == expected_pins


def test_mutation_space_freeze_all_pins_every_codon():
    space = CodingSpace('MIKEY')
    cds = space.sample()
    muts = MutationSpace(space, cds)

    muts.freeze_all()

    expected_pins = {
        pos: [cds[3 * (pos - 1): 3 * pos]]
        for pos in range(1, len(space.view.aa_seq) + 1)
    }

    assert space.view.pinned_codons == expected_pins


def test_mutation_space_unfreeze_all_clears_pins():
    space = CodingSpace('MIKEY')
    cds = space.sample()
    muts = MutationSpace(space, cds, free_positions=[])

    muts.unfreeze_all()
    assert space.view.pinned_codons == {}


def test_mutation_space_sample_uses_current_pins():
    space = CodingSpace('MIKEY')
    cds = space.sample()

    muts = MutationSpace(space, cds, free_positions=[2])
    variants = [muts.sample() for _ in range(1000)]
    assert len(set(variants)) == muts.n_valid_variants == 3
    for variant in variants:
        assert variant[0:3] == cds[0:3]
        assert variant[6:9] == cds[6:9]
        assert variant[9:12] == cds[9:12]
        assert variant[12:15] == cds[12:15]

    muts = MutationSpace(space, cds, free_positions=[2, 3])
    variants = [muts.sample() for _ in range(1000)]
    assert len(set(variants)) == muts.n_valid_variants == 6
    for variant in variants:
        assert variant[0:3] == cds[0:3]
        assert variant[9:12] == cds[9:12]
        assert variant[12:15] == cds[12:15]

    muts = MutationSpace(space, cds, free_positions=[5])
    variants = [muts.sample() for _ in range(1000)]
    assert len(set(variants)) == muts.n_valid_variants == 2
    for variant in variants:
        assert variant[0:3] == cds[0:3]
        assert variant[3:6] == cds[3:6]
        assert variant[6:9] == cds[6:9]
        assert variant[9:12] == cds[9:12]


def test_mutation_space_does_not_modify_original_space():
    space = CodingSpace('MIKEY')
    cds = space[0]

    muts = space.mutants(cds, free_positions=[2])
    muts.freeze_all()

    assert space.view.pinned_codons == {}
    assert muts.space.view.pinned_codons != {}


def test_mutation_space_rejects_cds_not_in_space():
    space = CodingSpace('MIKEY')

    with pytest.raises(ValueError):
        space.mutants('ATG' * 5)