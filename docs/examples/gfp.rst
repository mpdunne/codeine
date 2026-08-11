Designing with GFP
==================

In this example, we combine several features of Codeine in a single coding sequence
design workflow.

Here we will design coding sequences for expression of superfolder GFP (sfGFP),
a widely used fluorescent reporter, in *E. coli*.

Initialising the coding space
-----------------------------

Let's grab the amino acid sequence for sfGFP, from
`FPbase <https://www.fpbase.org/protein/superfolder-gfp/>`_.

.. code-block:: python

    sfgfp = (
        'MSKGEELFTGVVPILVELDGDVNGHKFSVRGEGEGDATNGKLTLKFICTTGKLPVPWPTL'
        'VTTLTYGVQCFSRYPDHMKRHDFFKSAMPEGYVQERTISFKDDGTYKTRAEVKFEGDTLV'
        'NRIELKGIDFKEDGNILGHKLEYNFNSHNVYITADKQKNGIKANFKIRHNVEDGSVQLAD'
        'HYQQNTPIGDGPVLLPDNHYLSTQSVLSKDPNEKRDHMVLLEFVTAAGITHGMDELYK'
    )

We'll initialise the coding space using codon weights for *E. coli*. These
bias sampling towards more commonly used synonymous codons, without restricting
which sequences belong to the coding space.

.. code-block:: python

    from codeine import CodingSpace, CodonWeights

    weights = CodonWeights.ecoli()

    space = CodingSpace(
        sfgfp,
        codon_weights=weights,
        seed=8675309,
    )

    unconstrained_n = space.count()

Adding constraints
------------------

We'll now restrict the coding space using a set of practical sequence design
constraints.

.. code-block:: python

    from codeine import RestrictionSite
    from codeine.constraints import (
        DirectRepeats,
        ForbiddenMotifs,
        Hairpins,
        InvertedRepeats,
        MaxHomopolymer,
        TandemRepeats,
    )

    constraints = [
        # Forbidden motifs
        ForbiddenMotifs([
            RestrictionSite.BsaI,     # Golden Gate
            RestrictionSite.BsmBI,    # Golden Gate
            RestrictionSite.BbsI,     # Golden Gate
            RestrictionSite.EcoRI,    # Conventional cloning
            RestrictionSite.BamHI,    # Conventional cloning
            RestrictionSite.HindIII,  # Conventional cloning
        ]),

        # Nucleotide homopolymers
        MaxHomopolymer(6),

        # Tandem repeats
        TandemRepeats(repeat_length=2, copies=4),
        TandemRepeats(repeat_length=3, copies=3),
        TandemRepeats(repeat_length=4, copies=3),
        TandemRepeats(repeat_length=5, copies=3),
        TandemRepeats(repeat_length=6, copies=3),

        # Longer direct and inverted repeats
        DirectRepeats(repeat_length=18),
        InvertedRepeats(repeat_length=18),

        # Hairpins
        Hairpins(stem_length=12, min_loop_length=3, max_loop_length=8),
    ]

    space.add_constraints(constraints)
    constrained_n = space.count()

    remaining = 100 * constrained_n / unconstrained_n

    print(f'{remaining:.2f}% of sequences remaining')

Output:

.. code-block:: text

    6.91% of sequences remaining

These constraints remove more than 93% of the synonymous sfGFP coding space!

Sampling sequences
------------------

We can sample coding sequences directly from the constrained space.

.. code-block:: python

    for sequence in space.sample(n=5):
        print(sequence)

Every sampled sequence encodes sfGFP and satisfies the specified constraints.
The *E. coli* codon weights bias sampling towards more commonly used codons,
without excluding less common synonymous codons from the space.

Local redesign
--------------

Once we have selected a valid sequence, we may want to explore nearby alternatives,
either limiting the number of changes or restricting them to particular regions
of the sequence.

Let's start by selecting a reference coding sequence at random:

.. code-block:: python

    reference = space.sample()

We can construct a mutation space containing every valid synonymous sequence
within two nucleotide changes of the reference.

.. code-block:: python

    mutants = space.mutants(reference, max_nts=2)

    print(mutants.count())

Output:

.. code-block:: text

    104358

Even within just two nucleotide changes of a single valid sfGFP coding
sequence, there are 104,358 synonymous variants that still satisfy all of the
original constraints.

Increasing the allowed mutation distance quickly expands the local sequence
space.

.. code-block:: python

    for max_nts in range(1, 6):
        mutants = space.mutants(reference, max_nts=max_nts)
        print(f'{max_nts} nt: {mutants.count():,} variants')

Output:

.. code-block:: text

    1 nt: 458 variants
    2 nt: 104,358 variants
    3 nt: 15,773,006 variants
    4 nt: 1,779,012,929 variants
    5 nt: 159,714,128,610 variants

We can also restrict redesign to positions that are experimentally accessible or
permitted to change. For example, here we allow changes only at positions 60-65.

.. code-block:: python

    mutants = space.mutants(
        reference,
        free_positions=range(60, 65),
        max_nts=5,
    )

    print(mutants.count())

Output:

.. code-block:: text

    1494

With only 1,494 valid variants remaining, this subspace can be enumerated in full rather than sampled:

.. code-block:: python

    for sequence in mutants:
        print(sequence)

The resulting sequences can then be evaluated exhaustively using whatever
objective or predictive model is appropriate for the application.