Mutation spaces
===============

A ``MutationSpace`` represents the subset of a coding space centred around a
single reference coding sequence. It can be constrained to mutate only certain
positions or to sample sequences with a given distance from the reference.

Creating a mutation space
-------------------------

To sample mutants, start with a ``CodingSpace`` and a reference coding sequence.

.. code-block:: python

    from codeine import CodingSpace

    space = CodingSpace('CYIQNCPLG')

    muts = space.mutants(
        'TGTTACATACAAAATTGTCCTCTAGGC',
        max_nts=2,
    )

    print(muts.n_valid_sequences)

Calling ``space.mutants(...)`` returns a ``MutationSpace`` object.
Alternatively the ``MutationSpace`` can be constructed directly:


.. code-block:: python

    from codeine import CodingSpace, MutationSpace

    space = CodingSpace('CYIQNCPLG')

    muts = MutationSpace(
        space=space,
        cds='TGTTACATACAAAATTGTCCTCTAGGC',
        max_nts=2,
    )

    print(muts.n_valid_sequences)

Here, ``max_nts=2`` allows coding sequences with at most two nucleotide
differences from the reference sequence.

Distance constraints
--------------------

Mutation spaces can be constrained by nucleotide distance, codon distance, or
both.

.. code-block:: python

    muts = space.mutants(
        reference_cds,
        min_nts=1,
        max_nts=3,
        min_codons=1,
        max_codons=2,
    )

``min_nts`` and ``max_nts`` count individual nucleotide changes (Hamming distance).

``min_codons`` and ``max_codons`` count codons that differ, regardless of how
many nucleotides differ within each codon. For example, ``ATGAAGTTT`` has codon distance ``1`` from ``ATGAAATTT``.

For longer sequences, the number of variants available can grow significantly as the distance constraints are relaxed. For example:

.. code-block:: python

    from codeine import CodingSpace

    space = CodingSpace('MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQANLAGKPE', seed=42)

    muts = space.mutants(cds=space.sample())

    muts.set_distance_constraints(max_nts=0)
    print(muts.n_valid_sequences) # 1

    muts.set_distance_constraints(max_nts=1)
    print(muts.n_valid_sequences) # 86

    muts.set_distance_constraints(max_nts=2)
    print(muts.n_valid_sequences) # 3599

    muts.set_distance_constraints(max_nts=3)
    print(muts.n_valid_sequences) # 97690

Constraining free positions
---------------------------

Use ``free_positions`` to choose which amino acid positions are allowed to vary.

.. code-block:: python

    muts = space.mutants(
        reference_cds,
        free_positions=[4, 5, 6, 7, 8],
        max_codons=2,
    )

    cds = muts.sample()

You can reset free position constraints using ``muts.unfreeze_all()``.

Contiguous regions can be supplied using ``range``.

.. code-block:: python

    muts = space.mutants(
        reference_cds,
        free_positions=range(20, 41),
        max_nts=5,
    )

Don't forget that positions are 1-based!

Enumerating & sampling variants
-------------------------------

Mutation spaces can be counted, enumerated, and sampled just like ordinary coding spaces.

.. code-block:: python

    from codeine import CodingSpace

    space = CodingSpace('CYIQNCPLG', seed=42)
    muts = space.mutants(
        cds=space.sample(),
        min_nts=2,
        max_nts=2,
    )

    print(muts.n_valid_sequences)

    print(space.sample())

    for cds in muts.enumerate():
        print(cds)

Inherited constraints
---------------------

Mutation spaces inherit constraints from the parent ``CodingSpace``.

.. code-block:: python

    from codeine import CodingSpace, RestrictionSite

    space = CodingSpace(
            aa_seq,
            forbidden_motifs=[
                RestrictionSite.EcoRI,
                RestrictionSite.NotI,
            ],
            max_homopolymer=5,
    )

    muts = space.mutants(
        reference_cds,
        free_positions=range(10, 31),
        max_codons=3,
    )

    cds = muts.sample()

The sampled sequences satisfy the mutation-space constraints and the
constraints from the original coding space.
