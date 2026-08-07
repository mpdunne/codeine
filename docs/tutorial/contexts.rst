Sequence contexts
=================

**Contexts** are regions to the left or right of a coding sequence that we wish to consider when avoiding forbidden motifs, but whose sequences remain fixed.

Sequence contexts are useful in two main scenarios:

* Designing genes for insertion into plasmids or expression vectors;
* Redesigning a small part of an existing coding sequence.

In both cases, forbidden sequences can span the junctions between the coding and context sequences, meaning that the contexts have to be considered
when designing sequences.

Sequence contexts can be supplied to ``CodingSpace`` via the ``context_l`` and
``context_r`` arguments.

.. code-block:: python

    from codeine import CodingSpace, RestrictionSite
    from codeine.constraints import ForbiddenMotifs

    space = CodingSpace(
        'CYIQNCPLG',
        context_l='GCGGCCGCGAAT',
        context_r='ATTCAGAATTCC',
        constraints=[
            ForbiddenMotifs(RestrictionSite.EcoRI)
        ],
    )

    print(space.n_valid_sequences)