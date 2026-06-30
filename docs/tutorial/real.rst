Examples
========

Here we consider some basic real world examples using **Codeine**

GFP
---

This example uses Codeine to explore synonymous coding sequences for GFP while
avoiding motifs that may interfere with cloning.

.. code-block:: python

   from codeine import CodingSpace, RestrictionSite

   GFP = '...'

   space = CodingSpace(
       GFP,
       forbidden_motifs=[
           RestrictionSite.EcoRI,
           RestrictionSite.BamHI,
           RestrictionSite.BsaI,
           RestrictionSite.BsmBI,
       ],
       max_homopolymer=5,
       seed=8675309,
   )

   print(space.n_valid_sequences)

   for _ in range(3):
       print(space.sample())

Example output:

.. code-block:: text
   ...

...

Codeine can represent very large coding spaces without enumerating every
sequence.

.. code-block:: python

   from codeine import CodingSpace, RestrictionSite

   SPCAS9 = '...'

   unconstrained = CodingSpace(SPCAS9)

   constrained = CodingSpace(
       SPCAS9,
       forbidden_motifs=[
           RestrictionSite.EcoRI,
           RestrictionSite.BamHI,
           RestrictionSite.BsaI,
           RestrictionSite.BsmBI,
       ],
       max_homopolymer=5,
       seed=42,
   )

   print(unconstrained.n_valid_sequences)
   print(constrained.n_valid_sequences)

   seq = constrained.sample()
   print(seq[:90] + '...')

Example output:

.. code-block:: text

   391000000000000000000000000000000000000000000000000000000000000000000000...
   140000000000000000000000000000000000000000000000000000000000000000000000...
   ATG...