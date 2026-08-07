Codeine
=======

**Codeine** is a Python library for exploring synonymous coding sequence spaces under biological constraints.


Why Codeine?
------------

Codeine compiles biological constraints into a ``CodingSpace``, an
exact representation of the feasible design space for a given amino
acid sequence.

Sequences sampled or enumerated from a ``CodingSpace`` are
automatically valid under the specified constraints.
This means sequences can be generated freely and scored according to
arbitrary objective functions, without sequence repair or filtering.


Features
--------

* Arbitrarily long amino acid sequences
* Exact counting of sequence spaces
* Sequence enumeration
* Uniform and weighted sampling
* Local sequence redesign
* Mutation spaces
* Fixed codons
* Sequence contexts
* Restriction sites & arbitrary motif avoidance
* Homopolymer and tandem repeat avoidance
* Custom translation tables


Example
-------

.. code-block:: python

    from codeine import CodingSpace, RestrictionSite
    from codeine.constraints import ForbiddenMotifConstraint, HomopolymerConstraint

    space = CodingSpace(
        'CYIQNCPLG',
        constraints=[
            ForbiddenMotifConstraint(RestrictionSite.EcoRI),
            HomopolymerConstraint(5),
        ],
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
