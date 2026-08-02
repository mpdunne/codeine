Adding constraints
==================

Sequence constraints are central to **Codeine**'s design philosophy, allowing you to easily narrow the space of valid coding sequences to match your requirements.

Common constraints are available directly as convenience arguments to
``CodingSpace``. More specialised constraints are specified via the
``constraints`` argument, and both approaches can be used together.

Currently **Codeine** handles:

* `Forbidden motifs`_
* `Max homopolymer length`_
* `Tandem repeats`_
* `Direct repeats`_
* `Inverted repeats`_
* `Exact hairpins`_

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

Forbidden motifs can either be passed directly via
``forbidden_motifs``, or as a
``ForbiddenMotifConstraint`` via the ``constraints`` argument.

For convenience, ``codeine.RestrictionSite`` provides a collection of commonly used restriction
enzyme recognition sequences. The built-in motifs can be used alongside custom motifs:

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

Homopolymer limits can either be specified via
``max_homopolymer``, or as a
``HomopolymerConstraint`` via the ``constraints`` argument.

.. _Tandem repeats:
.. _Direct repeats:
.. _Inverted repeats:
.. _Exact hairpins:

Repeat constraints
------------------

Repeat constraints are specified through the ``constraints`` argument to
``CodingSpace``:

.. code-block:: python

   from codeine import CodingSpace
   from codeine.constraints import (
       DirectRepeatConstraint,
       HairpinConstraint,
       InvertedRepeatConstraint,
       TandemRepeatConstraint,
   )

   space = CodingSpace(
       aa_seq,
       constraints=[
           TandemRepeatConstraint(4, 3),
           DirectRepeatConstraint(
               repeat_length=15,
               min_distance=30,
               max_distance=300,
           ),
           InvertedRepeatConstraint(
               repeat_length=15,
               min_distance=30,
               max_distance=300,
           ),
           HairpinConstraint(
               stem_length=12,
               min_loop=3,
               max_loop=20,
           ),
       ],
   )

   print(space.n_valid_sequences)

See :doc:`repeat_constraints` for details and additional examples.

Empty spaces
------------

It is possible for constraints to combine in such a way that there exist no valid coding sequences satisfying them. In this case, the coding space is considered to be empty, and cannot be sampled from.

For example:

.. code-block:: python

    from codeine import CodingSpace

    space = CodingSpace(
        'MKTLEFQNGSCPRYKKL',
        forbidden_motifs=[
            'ATGAAA',
            'ATGAAG',
        ],
    )

    print(f'Num. valid sequences: {space.n_valid_sequences}')

.. _combinations of the above:

Combining constraints
---------------------

Constraints can be combined freely. For example:

.. code-block:: python

   from codeine import CodingSpace, RestrictionSite
   from codeine.constraints import TandemRepeatConstraint

   aa_seq = 'MKTLEFQNGSCPRYKKL'

   space = CodingSpace(
       aa_seq,
       forbidden_motifs=[
           RestrictionSite.EcoRI,
           RestrictionSite.XhoI,
           'TAGATA',
       ],
       constraints=[
           TandemRepeatConstraint(4, 3),
       ],
       max_homopolymer=5,
       seed=42,
   )

   print(space.n_valid_sequences)
   print(space.sample())