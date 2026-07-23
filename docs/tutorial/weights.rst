Codon weights
=============

``CodonWeights`` control how codons are weighted during sampling.

They do not change which coding sequences are valid: they change how likely
different valid sequences are to be sampled. This can be useful if you want
to bias sequences towards the codon frequencies found in a specific organism.

If no ``CodonWeights`` object is passed to a ``CodingSpace``, then codons are sampled evenly.


Species-specific weights
------------------------

Codeine includes codon usage weights for a small set of common expression hosts
and model organisms.

.. code-block:: python

   from codeine import CodingSpace, CodonWeights

   weights = CodonWeights.ecoli()

   space = CodingSpace('CYIQNCPLG', codon_weights=weights)

The species weights currently available in **Codeine** are:

.. code-block:: text

   CodonWeights.arabidopsis()
   CodonWeights.drosophila()
   CodonWeights.ecoli()
   CodonWeights.human()
   CodonWeights.mouse()
   CodonWeights.yeast()

The species weights are based on codon usage tables from GenScript (`https://www.genscript.com/tools/codon-frequency-table <https://www.genscript.com/tools/codon-frequency-table>`_).

Custom weights
--------------

Custom weights can be supplied as a dictionary.

.. code-block:: python

   from codeine import CodonWeights

   weights = CodonWeights(
       {
           'A': {
               'GCA': 0.5,
               'GCC': 0.2,
               'GCG': 0.2,
               'GCT': 0.1,
           },
           'C': {
               'TGC': 0.6,
               'TGT': 0.4,
           },
           ...
       }
   )

   space = CodingSpace('CYIQNCPLG', codon_weights=weights)

The absolute value of the weights does not matter: only the relative weights
among synonymous codons affect sampling.

For example, these are equivalent:

.. code-block:: python

   {'TGC': 0.8, 'TGT': 0.2}
   {'TGC': 80, 'TGT': 20}

Weights can be specified using RNA or DNA codons, they will be converted internally.


Zero weights
------------

A codon can be assigned a weight of zero. The codon remains part of the valid
coding space, but it will never be selected during weighted sampling.

.. code-block:: python

    weights = CodonWeights(
        {
            ...
            'C': {
                    'TGC': 0.0,
                    'TGT': 1.0,
                },
            ...
        }
    )

At least one codon for each amino acid must have a positive weight.

Weight thresholds
-----------------

Low-frequency codons can be excluded from sampling using ``threshold()``.

.. code-block:: python

    weights = CodonWeights.ecoli().threshold(0.1)

A threshold of ``0.1`` prevents codons representing less than 10% of
the synonymous codon weight for an amino acid from being sampled.

Thresholding does not remove sequences from the coding space and does not
affect counting, enumeration, indexing and membership tests. To remove zero-weight codons from the coding space
entirely, generate an updated translation table using the ``restrict()`` method:

.. code-block:: python

    from codeine import CodingSpace, CodonWeights, TranslationTable

    table = TranslationTable()
    weights = CodonWeights.ecoli().threshold(0.1)

    table, weights = weights.restrict(table)

    space = CodingSpace('CYIQNCPLG', translation_table=table, codon_weights=weights)
