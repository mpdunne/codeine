import pytest

from codeine.space.coding import CodingSpace
from codeine.space.mutation import MutationSpace


def test_mutation_space_requires_cds_in_space():
    space = CodingSpace('MIKEY')

    with pytest.raises(ValueError):
        MutationSpace(space=space, cds='AAAAAAAAAAAA')

    with pytest.raises(ValueError):
        MutationSpace(space=space, cds='')

    with pytest.raises(ValueError):
        MutationSpace(space=space, cds='MIKEY')

    for seq in space:
        _ = MutationSpace(space=space, cds=seq)
        _ = MutationSpace(space=space, cds=seq.lower())
        _ = MutationSpace(space=space, cds=seq.upper())


def test_mutation_space_inherits_existing_pins():
    space = CodingSpace('MIKEY')
    cds = space[0]
    space.pin_codons({3: cds[6:9]})
    mutation_space = space.mutants(cds)
    assert mutation_space.view.pinned_codons[3] == [cds[6:9]]
    assert mutation_space.contains(cds)

    space = CodingSpace('MIKEY')
    cds = space[0]
    space.pin_codons({3: cds[6:9]})
    mutation_space = MutationSpace(space, cds)
    assert mutation_space.view.pinned_codons[3] == [cds[6:9]]
    assert mutation_space.contains(cds)


def test_mutation_space_remains_unchanged_if_space_changes():
    space = CodingSpace('MIKEY')
    cds = space[0]
    space.pin_codons({3: cds[6:9]})
    mutation_space = space.mutants(cds)
    space.pin_codons({4: cds[9:12]})
    assert 4 not in mutation_space.view.pinned_codons

    space = CodingSpace('MIKEY')
    cds = space[0]
    space.pin_codons({3: cds[6:9]})
    mutation_space = MutationSpace(space, cds)
    space.pin_codons({4: cds[9:12]})
    assert 4 not in mutation_space.view.pinned_codons


def test_mutation_space_combines_original_pins_with_frozen_positions():
    space = CodingSpace('MIKEY')
    cds = space[0]
    space.pin_codons({3: cds[6:9]})
    mutation_space = space.mutants(cds, free_positions=[2, 3])
    expected_pins = {1: [cds[0:3]], 3: [cds[6:9]], 4: [cds[9:12]], 5: [cds[12:15]]}
    assert mutation_space.view.pinned_codons == expected_pins

    space = CodingSpace('MIKEY')
    cds = space[0]
    space.pin_codons({3: cds[6:9]})
    mutation_space = MutationSpace(space, cds, free_positions=[2, 3])
    expected_pins = {1: [cds[0:3]], 3: [cds[6:9]], 4: [cds[9:12]], 5: [cds[12:15]]}
    assert mutation_space.view.pinned_codons == expected_pins


def test_mutation_space_original_space_untouched_by_freezing():
    space = CodingSpace('MIKEY')
    cds = space[0]
    space.pin_codons({3: cds[6:9]})
    original_pins = dict(space.view.pinned_codons)
    mutation_space = space.mutants(cds, free_positions=[1, 3, 4])
    assert space.view.pinned_codons == original_pins
    assert mutation_space.view.pinned_codons != original_pins

    space = CodingSpace('MIKEY')
    cds = space[0]
    space.pin_codons({3: cds[6:9]})
    original_pins = dict(space.view.pinned_codons)
    mutation_space = MutationSpace(space, cds, free_positions=[1, 3, 4])
    assert space.view.pinned_codons == original_pins
    assert mutation_space.view.pinned_codons != original_pins


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
    assert muts.view.pinned_codons == expected_pins


def test_mutation_space_freeze_all_pins_every_codon():
    space = CodingSpace('MIKEY')
    cds = space.sample()
    muts = MutationSpace(space, cds)

    muts.freeze_all()

    expected_pins = {
        pos: [cds[3 * (pos - 1): 3 * pos]]
        for pos in range(1, len(space.view.aa_seq) + 1)
    }

    assert muts.view.pinned_codons == expected_pins


def test_mutation_space_unfreeze_all_clears_pins():
    space = CodingSpace('MIKEY')
    cds = space.sample()
    muts = MutationSpace(space, cds, free_positions=[])

    muts.unfreeze_all()
    assert space.view.pinned_codons == {}
    assert muts.view.pinned_codons == {}


def test_mutation_space_samples_correctly():
    space = CodingSpace('MIKEY')
    cds = space.sample()

    muts = MutationSpace(space, cds, free_positions=[2])
    assert space.view.pinned_codons == {}
    assert set(muts.view.pinned_codons.keys()) == {1, 3, 4, 5}
    variants = [muts.sample() for _ in range(1000)]
    assert len(set(variants)) == muts.n_valid_variants == 3
    for variant in variants:
        assert variant[0:3] == cds[0:3]
        assert variant[6:9] == cds[6:9]
        assert variant[9:12] == cds[9:12]
        assert variant[12:15] == cds[12:15]

    muts = MutationSpace(space, cds, free_positions=[2, 3])
    assert space.view.pinned_codons == {}
    assert set(muts.view.pinned_codons.keys()) == {1, 4, 5}
    variants = [muts.sample() for _ in range(1000)]
    assert len(set(variants)) == muts.n_valid_variants == 6
    for variant in variants:
        assert variant[0:3] == cds[0:3]
        assert variant[9:12] == cds[9:12]
        assert variant[12:15] == cds[12:15]

    muts = MutationSpace(space, cds, free_positions=[5])
    assert space.view.pinned_codons == {}
    assert set(muts.view.pinned_codons.keys()) == {1, 2, 3, 4}
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
    assert muts.view.pinned_codons != {}


def test_mutation_space_has_no_distance_constraints_by_default():
    space = CodingSpace('MIKEY')
    cds = space.sample()
    muts = MutationSpace(space, cds)

    assert muts.min_nts is None
    assert muts.max_nts is None
    assert muts.min_codons is None
    assert muts.max_codons is None
    assert not muts.has_distance_constraints


def test_mutation_space_accepts_distance_constraints_in_constructor():
    space = CodingSpace('MIKEY')
    cds = space.sample()
    muts = MutationSpace(space, cds, min_nts=3, max_nts=4, min_codons=1, max_codons=2)

    assert muts.min_nts == 3
    assert muts.max_nts == 4
    assert muts.min_codons == 1
    assert muts.max_codons == 2
    assert muts.has_distance_constraints


def test_mutation_space_accepts_exact_nt_distance_in_constructor():
    space = CodingSpace('MIKEY')
    cds = space.sample()
    muts = MutationSpace(space, cds, min_nts=3, max_nts=5)

    assert muts.min_nts == 3
    assert muts.max_nts == 5


def test_mutation_space_accepts_exact_codon_distance_in_constructor():
    space = CodingSpace('MIKEY')
    cds = space.sample()
    muts = MutationSpace(space, cds, min_codons=3, max_codons=5)

    assert muts.min_codons == 3
    assert muts.max_codons == 5


def test_set_distance_sets_nt_range():
    space = CodingSpace('MIKEY')
    cds = space.sample()
    muts = MutationSpace(space, cds)

    muts.set_distance_constraints(min_nts=3, max_nts=4)

    assert muts.min_nts == 3
    assert muts.max_nts == 4
    assert muts.min_codons is None
    assert muts.max_codons is None


def test_set_distance_sets_codon_range():
    space = CodingSpace('MIKEY')
    cds = space.sample()
    muts = MutationSpace(space, cds)

    muts.set_distance_constraints(min_codons=3, max_codons=10)

    assert muts.min_codons == 3
    assert muts.max_codons == 10
    assert muts.min_nts is None
    assert muts.max_nts is None


def test_set_distance_can_set_nt_and_codon_constraints_together():
    space = CodingSpace('MIKEY')
    cds = space.sample()
    muts = MutationSpace(space, cds)

    muts.set_distance_constraints(min_nts=3, max_nts=6, min_codons=2, max_codons=3)

    assert muts.min_nts == 3
    assert muts.max_nts == 6
    assert muts.min_codons == 2
    assert muts.max_codons == 3


@pytest.mark.parametrize(
    'kwargs',
    [
        {'min_nts': -1},
        {'max_nts': -1},
        {'min_codons': -1},
        {'max_codons': -1},
    ],
)
def test_set_distance_rejects_negative_distances(kwargs):
    space = CodingSpace('MIKEY')
    cds = space.sample()
    muts = MutationSpace(space, cds)
    with pytest.raises(ValueError):
        muts.set_distance_constraints(**kwargs)


def test_set_distance_rejects_min_nts_greater_than_max_nts():
    space = CodingSpace('MIKEY')
    cds = space.sample()
    muts = MutationSpace(space, cds)
    with pytest.raises(ValueError):
        muts.set_distance_constraints(min_nts=5, max_nts=3)


def test_set_distance_rejects_min_codons_greater_than_max_codons():
    space = CodingSpace('MIKEY')
    cds = space.sample()
    muts = MutationSpace(space, cds)
    with pytest.raises(ValueError):
        muts.set_distance_constraints(min_codons=5, max_codons=3)


def test_set_distance_with_no_args_clears_distance_constraints():
    space = CodingSpace('MIKEY')
    cds = space.sample()
    muts = MutationSpace(space, cds, min_nts=1, max_nts=10, min_codons=1, max_codons=5)

    muts.set_distance_constraints()
    assert muts.min_nts is None
    assert muts.max_nts is None
    assert muts.min_codons is None
    assert muts.max_codons is None
    assert not muts.has_distance_constraints


def test_clear_distance_constraints():
    space = CodingSpace('MIKEY')
    cds = space.sample()
    muts = MutationSpace(space, cds, min_nts=1, max_nts=10, min_codons=1, max_codons=5)

    muts.clear_distance_constraints()
    assert muts.min_nts is None
    assert muts.max_nts is None
    assert muts.min_codons is None
    assert muts.max_codons is None
    assert not muts.has_distance_constraints


def test_distance_constraints_raise_on_methods():
    space = CodingSpace('MIKEY')
    cds = space.sample()
    muts = MutationSpace(space, cds, max_nts=5)

    with pytest.raises(NotImplementedError):
        _ = muts.n_valid_variants

    with pytest.raises(NotImplementedError):
        _ = muts.contains(cds)

    with pytest.raises(NotImplementedError):
        _ = muts.sample()

    with pytest.raises(NotImplementedError):
        _ = [*muts.enumerate()]


def test_distance_constraints_raise_for_contains():
    space = CodingSpace('MIKEY')
    cds = space.sample()
    muts = MutationSpace(space, cds, max_nts=5)
    with pytest.raises(NotImplementedError):
        _ = cds in muts


def test_distance_constraints_raise_for_iter():
    space = CodingSpace('MIKEY')
    cds = space.sample()
    muts = MutationSpace(space, cds, max_nts=5)
    with pytest.raises(NotImplementedError):
        _ = [*muts]


def test_distance_constraints_raise_for_getitem():
    space = CodingSpace('MIKEY')
    cds = space.sample()
    muts = MutationSpace(space, cds, max_nts=5)
    with pytest.raises(NotImplementedError):
        _ = muts[0]
