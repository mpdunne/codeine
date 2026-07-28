Adding constraints
==================

Sequence constraints are central to **Codeine**'s design philosophy, allowing you to easily narrow the space of valid coding sequences to match your requirements.

Currently **Codeine** handles:

* `Forbidden motifs`_
* `Max homopolymer length`_
* `Tandem repeats`_

As well as arbitrary `combinations of the above`_. Let's go through these!

.. _Forbidden motifs:

Forbidden motifs
----------------

In biotechnology we often wish to exclude sequence motifs that interfere with cloning, synthesis, or downstream applications.

**Codeine** can exclude such sequences via ``forbidden_motifs``:

.. code-block:: python

    from codeine import CodingSpace, RestrictionSite

    aa_seq = 'MKTLEFQNGSCPRYKKL'

    unconstrained = CodingSpace(aa_seq)

    constrained = CodingSpace(
       aa_seq,
       forbidden_motifs=[
           'GAATTC',
           'CTCGAG',
       ],
       seed=42,
    )

    print(f'Without constraints: {unconstrained.n_valid_sequences:,} sequences')
    print(f'With constraints: {constrained.n_valid_sequences:,} sequences')
    print(constrained.sample())


For convenience, ``codeine.RestrictionSite`` provides a collection of commonly used restriction
enzyme recognition sequences. For example:

.. code-block:: python

   from codeine import RestrictionSite

   print(RestrictionSite.EcoRI)
   print(RestrictionSite.BsaI)


This outputs a summary of each of these restriction sites:

.. code-block:: text

    EcoRI (GAATTC)
    BsaI (GGTCTC / GAGACC)

The built-in motifs can then be used alongside custom motifs:

.. code-block:: python

   from codeine import CodingSpace, RestrictionSite

   space = CodingSpace(
       aa_seq,
       forbidden_motifs=[
           RestrictionSite.EcoRI,
           RestrictionSite.BsaI,
           'GATTACA',
       ],
   )

Each named restriction site is automatically interpreted to mean both the recognition sequence and its
reverse complement.

The built-in restriction sites are:

* **BioBricks**: ``EcoRI``, ``XbaI``, ``SpeI``, ``PstI``;
* **Common cloning**: ``BamHI``, ``HindIII``, ``XhoI``, ``SalI``, ``KpnI``, ``SacI``, ``NcoI``, ``NdeI``, ``NotI``, ``MluI``, ``AgeI``, ``AvrII``, ``BglII``;
* **Golden Gate**: ``BsaI``, ``BsmBI``, ``BbsI``, ``SapI``.

The recognition sequences for these were taken from the
`New England Biolabs alphabetized list of recognition sequences
<https://www.neb.com/en/tools-and-resources/selection-charts/alphabetized-list-of-recognition-specificities>`_.

.. _Max homopolymer length:


Max homopolymer
---------------

Long homopolymer runs can be difficult to work with because they increase the risk of polymerase slippage and sequencing errors. In **Codeine**, runs of identical nucleotides can be avoided using the ``max_homopolymer`` constraint.

For example, setting ``max_homopolymer=5`` excludes any coding sequence
containing six or more consecutive identical nucleotides.


.. code-block:: python

   from codeine import CodingSpace

   space = CodingSpace(
       aa_seq,
       max_homopolymer=5,
   )

   print(space.n_valid_sequences)


.. _Tandem repeats:

Tandem repeats
--------------

Tandem repeats can increase the risk of polymerase slippage during DNA
replication. In **Codeine**, tandem repeats can be avoided using the
``TandemRepeatConstraint``.

For example, the following constraint forbids tandem repeats with a
repeat unit of 4 nucleotides occurring three or more times consecutively.

.. code-block:: python

   from codeine import CodingSpace
   from codeine.constraints import TandemRepeatConstraint

   space = CodingSpace(
       aa_seq,
       constraints=[
           TandemRepeatConstraint(repeat_length=4, min_copies=3),
       ],
   )

   print(space.n_valid_sequences)

Tandem repeats are, by design, exact. Different repeat lengths are
therefore handled by separate constraints, which can be combined as
desired:


.. code-block:: python

   space = CodingSpace(
       aa_seq,
       constraints=[
           TandemRepeatConstraint(2, 4),  # e.g. ATATATAT
           TandemRepeatConstraint(3, 3),  # e.g. GTCGTCGTC
           TandemRepeatConstraint(4, 3),  # e.g. AGCTAGCTAGCT
       ],
   )

.. _combinations of the above:

Combining constraints
---------------------

Constraints can be combined freely!

.. code-block:: python

   from codeine import CodingSpace, RestrictionSite

   aa_seq = 'MKTLEFQNGSCPRYKKL'

   space = CodingSpace(
       aa_seq,
       codon_restrictions={
           2: ['AAA', 'AAG'],
           5: 'GAA',
       },
       forbidden_motifs=[
           RestrictionSite.EcoRI,
           RestrictionSite.XhoI,
           'TAGATA',
       ],
       max_homopolymer=5,
       seed=42,
   )

   print(space.n_valid_sequences)
   print(space.sample())

Empty spaces
------------

It is possible for constraints to combine in such a way that there exist no valid coding sequences satisfying them. In this case, the coding space is considered to be empty, and cannot be sampled from.

For example:

.. code-block:: python

    from codeine import CodingSpace, RestrictionSite

    space = CodingSpace(
       'MKTLEFQNGSCPRYKKL',
       forbidden_motifs=[
           'ATGAAA',
           'ATGAAG',
       ],
    )

    print(f'Num. valid sequences: {space.n_valid_sequences} ')