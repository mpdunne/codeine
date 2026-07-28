Coding spaces
=============

Here we introduce the three basic operations on a coding space: counting,
sampling and enumeration.

Counting
--------

To create a coding space and count sequences:

.. code-block:: python

    from codeine import CodingSpace

    space = CodingSpace('CYIQNCPLG')

    print(f"{space.n_valid_sequences} valid coding sequences")

Output:

.. code-block:: text

   9216 valid coding sequences

Sampling
--------
Coding sequences for the chosen amino acid sequence can be obtained by randomly sampling the space. For example:

.. code-block:: python

    from codeine import CodingSpace

    space = CodingSpace('CYIQNCPLG', seed=42)

    print(space.sample())

Output:

.. code-block:: text

    TGTTACATACAAAATTGTCCTCTAGGC

In the example above, ``seed`` is provided to ensure reproducible sampling across runs.

You can also sample batches of sequences:

.. code-block:: python

    for seq in space.sample(5):
        print(seq)

Output:

.. code-block:: text

    TGCTACATCCAAAACTGTCCGCTCGGG
    TGTTACATTCAGAACTGCCCTCTGGGA
    TGCTATATCCAGAATTGTCCTCTGGGG
    TGTTATATTCAGAATTGCCCACTCGGA
    TGCTACATACAGAACTGCCCACTCGGT


Enumeration
-----------

Coding spaces can be iterated over or enumerated. Full enumeration is only practical for relatively small spaces. You can do:

.. code-block:: python

    from codeine import CodingSpace

    space = CodingSpace('CYIQNCPLG')

    for seq in space.enumerate():
        print(seq)

Or simply:

.. code-block:: python

    for seq in space:
       print(seq)

Output:

.. code-block:: text

    TGCTACATACAAAACTGCCCACTAGGA
    TGCTACATACAAAACTGCCCACTAGGC
    TGCTACATACAAAACTGCCCACTAGGG
    TGCTACATACAAAACTGCCCACTAGGT
    TGCTACATACAAAACTGCCCACTCGGA
    TGCTACATACAAAACTGCCCACTCGGC
   ...

You can also just grab the first few sequences:

.. code-block:: python

    from codeine import CodingSpace

    space = CodingSpace('CYIQNCPLG')

    for seq in space[:3]:
       print(seq)

Output:

.. code-block:: text

    TGCTACATACAAAACTGCCCACTAGGA
    TGCTACATACAAAACTGCCCACTAGGC
    TGCTACATACAAAACTGCCCACTAGGG


Validation
----------

The ``in`` operator can be used to test whether a coding sequence belongs to a coding space. Membership testing automatically accounts for any specified constraints and does not require enumerating the coding space.

.. code-block:: python

    space = CodingSpace('CYIQNCPLG', fixed_codons={2: 'TAC'})

    print('TGTTACATACAGAACTGCCCTCTTGGG' in space)
    # True

    print('TGTTATATACAGAACTGCCCTCTTGGG' in space)
    # False
