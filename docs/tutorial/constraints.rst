Adding constraints
==================

Constraints reduce the space of valid coding sequences for a given amino acid sequence. Currently **Codeine** handles:

* `Forbidden motifs`_
* `Max homopolymer length`_
* `Restricted codon choices`_

As well as arbitrary `combinations of the above`_. Let's go through these!

.. _Forbidden motifs:

Forbidden motifs
----------------

In biotechnology it is common to exclude DNA sequence motifs that interfere with cloning, synthesis, or downstream applications. **Codeine** can exclude such sequences automatically.

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

The ``max_homopolymer`` constraint limits runs of identical nucleotides.

For example, consider the amino acid sequence ``FFFFF``. In the standard
genetic code, phenylalanine is encoded by either ``TTT`` or ``TTC``. Without
constraints, some coding sequences contain long runs of thymine:

.. code-block:: text

   TTT TTT TTT TTT TTT

Setting ``max_homopolymer=5`` excludes any coding sequence containing six or
more consecutive identical nucleotides.

.. code-block:: python

   from codeine import CodingSpace

   space = CodingSpace(
       'FFFFF',
       max_homopolymer=5,
   )

   print(space.n_valid_sequences)

After applying the constraint, 13 of the 32 possible coding sequences remain valid.

.. _Restricted codon choices:

Fixed codons
------------
``codon_restrictions`` restricts which codons are allowed at specific amino acid positions. They can either be fixed exact codons, or subsets of the set of possible codons for that amino acid. Positions are 1-based.

.. code-block:: python

    from codeine import CodingSpace

    aa_seq = 'SEQVENCE'

    space = CodingSpace(
       aa_seq,
       codon_restrictions={
           1: ['TCG', 'TCA'],
           2: 'GAG',
       },
       seed=42,
    )

Here, position 1 is restricted to ``TCG`` or ``TCA``, and position 2 is fixed to
``GAG``.

.. _combinations of the above:

Combining constraints
---------------------

Constraints can be combined freely.

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
