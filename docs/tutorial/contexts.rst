Sequence contexts
=================

**Contexts** are regions to the left or right of a coding sequence that we wish to consider when avoiding forbidden motifs, but whose sequences remain fixed.

Sequence contexts are useful in two main scenarios:

* Designing genes for insertion into plasmids or expression vectors,
  where motifs may span the junction between the vector and the coding sequence.
* Redesigning only part of an existing coding sequence while leaving the
  surrounding sequence unchanged.

Sequence contexts can be supplied to ``CodingSpace`` via the ``context_l`` and
``context_r`` arguments.