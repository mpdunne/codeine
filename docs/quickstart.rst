Quickstart
==========

To get going with **Codeine**, create a ``CodingSpace``, and start exploring:

.. code-block:: python

   from codeine import CodingSpace

   space = CodingSpace('SEQVENCE')

   print(space.n_valid_sequences)
   print(space.sample())

Example output:

.. code-block:: text

   1536
   AGCGAACAGGTCGAGAACTGCGAA
