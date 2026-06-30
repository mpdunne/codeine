Codon weights
=============

``CodonWeights`` control how codons are weighted during sampling.

They do not change which coding sequences are valid. They only change how likely
different valid sequences are to be sampled and can be useful if you want to bias sequences towards a codon frequency found in a specific organism.

If no ``CodonWeights`` object is passed to a ``CodingSpace``, then codons are sampled evenly.


Species-specific weights
------------------------

Codeine includes codon usage weights for a small set of common expression hosts
and model organisms.

.. code-block:: python

   from codeine import CodingSpace, CodonWeights

   weights = CodonWeights.ecoli()

   space = CodingSpace('SEQVENCE', codon_weights=weights)

The species weights currently available in **Codeine** are:

.. code-block:: text

   CodonWeights.arabidopsis()
   CodonWeights.drosophila()
   CodonWeights.ecoli()
   CodonWeights.human()
   CodonWeights.mouse()
   CodonWeights.yeast()

The species weights are based on codon count tables from GenScript (`https://www.genscript.com/tools/codon-frequency-table <https://www.genscript.com/tools/codon-frequency-table>`_).

Custom weights
--------------

Custom weights can be supplied as a dictionary.

.. code-block:: python

   from codeine import CodonWeights

   weights = CodonWeights(
       {
           'F': {
               'TTT': 0.8,
               'TTC': 0.2,
           },
           'L': {
               'TTA': 0.1,
               'TTG': 0.2,
               'CTT': 0.2,
               'CTC': 0.2,
               'CTA': 0.1,
               'CTG': 0.2,
           },
           ...
       }
   )

   space = CodingSpace('SEQVENCE', codon_weights=weights)

The absolute value of the weights does not matter: only the relative weights
among synonymous codons affect sampling.

For example, these are equivalent:

.. code-block:: python

   {'TTT': 0.8, 'TTC': 0.2}
   {'TTT': 80, 'TTC': 20}
