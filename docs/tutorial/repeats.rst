Repeat constraints
==================

Repetitive nucleotide sequences can complicate DNA synthesis, sequencing, and protein production.
As a result, they are often avoided in synthetic biology applications.

**Codeine** provides constraints for excluding:

* `Tandem repeats`_
* `Direct repeats`_
* `Inverted repeats`_
* `Hairpins`_

Let's take a closer look!

.. _repeats-tandem:

Tandem repeats
--------------

Tandem repeats are stretches of DNA that contain the same sequence repeated multiple
times without gaps. Tandem repeats can increase the risk of polymerase slippage, repeat expansion,
and DNA synthesis failures.

The following example uses ``TandemRepeats`` to forbid tandem repeats with a repeat unit of
4 nucleotides occurring three or more times consecutively.

.. code-block:: python

   from codeine import CodingSpace
   from codeine.constraints import TandemRepeats

   space = CodingSpace(
       aa_seq,
       constraints=[
           TandemRepeats(
               repeat_length=4,
               copies=3,
           ),
       ],
   )

   print(space.n_valid_sequences)

Tandem repeats are, by design, exact. Different repeat lengths are therefore handled
by separate constraints, which can be combined as desired:

.. code-block:: python

   space = CodingSpace(
       aa_seq,
       constraints=[
           TandemRepeats(2, 4),  # e.g. ATATATAT
           TandemRepeats(3, 3),  # e.g. GTCGTCGTC
           TandemRepeats(4, 3),  # e.g. AGCTAGCTAGCT
       ],
   )

.. _repeats-direct:

Direct repeats
--------------

A direct repeat is a stretch of nucleotides that is repeated exactly elsewhere in the sequence.
Direct repeats can increase the risk of homologous recombination, leading to deletions and
reduced genetic stability.

The following example uses ``DirectRepeats`` to exclude direct repeats of length 20
nucleotides occurring within a specified distance range.

.. code-block:: python

   from codeine import CodingSpace
   from codeine.constraints import DirectRepeats

   space = CodingSpace(
       aa_seq,
       constraints=[
           DirectRepeats(
               repeat_length=20,
               min_distance=20,
               max_distance=200,
           ),
       ],
   )

   print(space.n_valid_sequences)

Here, the distance is measured between the end of the first repeated sequence and the start of the second.


.. _repeats-inverted:

Inverted repeats
----------------

Inverted repeats are stretches of nucleotides that are exact reverse-complements of other
stretches in the sequence. These can form undesirable secondary structures.

The syntax for ``InvertedRepeats`` is the same as for ``DirectRepeats``:


.. code-block:: python

   from codeine import CodingSpace
   from codeine.constraints import InvertedRepeats

   space = CodingSpace(
       aa_seq,
       constraints=[
           InvertedRepeats(
               repeat_length=12,
               min_distance=20,
               max_distance=200,
           ),
       ],
   )

   print(space.n_valid_sequences)


.. _repeats-hairpins:

Hairpins
--------------

Hairpins are a special case of inverted repeats, where the gap between complementary sequences
is small, forming highly stable "stem and loop" structures.

These can be excluded using with ``Hairpins``:

.. code-block:: python

   from codeine import CodingSpace
   from codeine.constraints import Hairpins

   space = CodingSpace(
       aa_seq,
       constraints=[
           Hairpins(
               stem_length=12,
               min_loop=3,
               max_loop=20,
           ),
       ],
   )

   print(space.n_valid_sequences)

Note that ``Hairpins`` requires exact base pairing: wobble bases are not considered.


Combining repeat constraints
----------------------------

Repeat constraints can be combined freely.

.. code-block:: python

   from codeine.constraints import (
       DirectRepeats,
       Hairpins,
       InvertedRepeats,
       TandemRepeats,
   )

   space = CodingSpace(
       aa_seq,
       constraints=[
           TandemRepeats(
                repeat_length=4,
                copies=3,
           ),
           DirectRepeats(
               repeat_length=15,
               min_distance=30,
               max_distance=300,
           ),
           InvertedRepeats(
               repeat_length=15,
               min_distance=30,
               max_distance=300,
           ),
           Hairpins(
               stem_length=12,
               min_loop=3,
               max_loop=20,
           ),
       ],
   )

Performance
-----------

Compilation time depends on the amino acid sequence and the chosen constraint parameters.
For long or highly repetitive sequences, direct and inverted repeat constraints can
substantially increase compilation time.