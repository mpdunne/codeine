Basic usage
===========

Here we introduce the three basic operations on a coding space: counting,
sampling and enumeration.

Counting
--------

To create a coding space and count sequences:

.. code-block:: python

    from codeine import CodingSpace

    space = CodingSpace('SEQVENCE')

    print(f"{space.n_valid_sequences} valid coding sequences")

Output:

.. code-block:: text

   1536 valid coding sequences


Sampling
--------
Coding sequences for the chosen amino acid sequence can be obtained by randomly samping the space. For example:

.. code-block:: python

    from codeine import CodingSpace

    space = CodingSpace('SEQVENCE', seed=42)

    print(space.sample())

Output:

.. code-block:: text

    TCGGAACAAGTTGAGAACTGCGAA

You can also sample batches of sequences:

.. code-block:: python

    for seq in space.sample(5):
        print(seq)

Output:

.. code-block:: text

    TCGGAACAAGTTGAGAACTGCGAA
    TCAGAACAAGTAGAAAATTGCGAG
    TCCGAGCAGGTTGAGAACTGTGAA
    AGCGAACAAGTTGAGAACTGCGAG
    TCGGAGCAAGTAGAGAACTGCGAG


Enumeration
-----------

Coding spaces can be iterated over or enumerated. Full enumeration is only practical for relatively small spaces. You can do:

.. code-block:: python

    from codeine import CodingSpace

    space = CodingSpace('MILK')

    for seq in space.enumerate():
        print(seq)

Or simply:

.. code-block:: python

    for seq in space:
       print(seq)

Output:

.. code-block:: text

    ATGATTTTAAAA
    ATGATTTTAAAG
    ATGATTTTGAAA
    ATGATTTTGAAG
    ATGATTCTTAAA
    ATGATTCTTAAG
   ...

You can also just grab the first few sequences:

.. code-block:: python

    from codeine import CodingSpace

    space = CodingSpace('SEQVENCE')

    for seq in space[:3]:
       print(seq)

Output:

.. code-block:: text

    ATGATTTTAAAA
    ATGATTTTAAAG
    ATGATTTTGAAA
