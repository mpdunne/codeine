Codeine
=======

**Codeine** is a Python library for exploring constrained
synonymous coding sequence spaces.

Overview
--------

**Codeine** represents the complete space of synonymous
coding sequences satisfying user-defined constraints. This makes it
possible to count, sample and enumerate valid coding sequences exactly.

Traditional codon optimisation tools typically return a single
sequence. By contrast, Codeine represents the full feasible
sequence space, allowing users to explore alternative designs freely.

Constraints are built in to **Codeine**'s sequence generation process,
so that all generated sequences are automatically valid. This means you can generate as
many sequences as needed and score them according to your own objective
functions, without having to repair or filter them.


Features
--------

* Exact counting of sequence spaces
* Uniform and weighted sampling
* Sequence enumeration
* Local sequence redesign
* Mutation spaces
* Restriction site avoidance
* Homopolymer constraints
* Custom translation tables

Example
-------

.. code-block:: python

    from codeine import CodingSpace, RestrictionSite

    space = CodingSpace(
        'CYIQNCPLG',
        forbidden_motifs=[RestrictionSite.EcoRI],
        max_homopolymer=5
    )

    print(space.n_valid_sequences)
    print(space.sample())



.. toctree::
   :maxdepth: 2
   :caption: Contents:
   :hidden:
   :includehidden:

   installation
   quickstart
   tutorial/index
   api
