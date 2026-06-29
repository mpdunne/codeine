Quickstart
==========

Codeine represents the set of coding sequences that translate to a given amino
acid sequence.

.. code-block:: python

   from codeine import CodingSpace

   space = CodingSpace('SEQVENCE')

   print(space.n_valid_sequences)
   print(space.sample())

Example output:

.. code-block:: text

   1536
   AGCGAACAGGTCGAGAACTGCGAA
