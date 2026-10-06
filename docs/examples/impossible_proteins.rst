Impossible coding spaces
========================

A protein can have an enormous number of synonymous coding sequences, while a
set of sequence constraints can still make its valid coding space empty. Here,
we use Codeine to identify two examples where local amino acid sequences make a
set of otherwise reasonable design constraints impossible to satisfy.

Restriction sites can make a coding space empty
-----------------------------------------------

As an example, consider the CRISPR-associated endonuclease Cas1/endonuclease
Cas2 `UniProt G4RJY6 <https://www.uniprot.org/uniprotkb/G4RJY6/entry>`_. Suppose
our cloning workflow requires the final coding sequence to avoid the ``NdeI``
and ``PciI`` restriction sites.

.. code-block:: python

    protein_sequence = (
        'MDEVLLLTGGISITTRALRALLATGATVAVFSPRGEPLGIFMRPVGDATGAKRRCQYKAA'
        'EDGRGLQYAKSWVFKKILGQRDNIKAWRRRLRGYSQYAESLAKALPGAGLHGAMETPRRR'
        'RRGQDGGQAGVRGRPTHPPVPPGAGRRSPGGAPRGQEASLRRDPQRGQSSGALHMYVIVV'
        'YDITENDVRAKVADILRAYGLARIQRSAYVGRLPPALVKELAERLARAVRGANADIAIFK'
        'VDKRTIDTSLRIPPRPPAGHVALH'
    )

First, let's look at the unrestricted coding space:

.. code-block:: python

    from codeine import CodingSpace

    unrestricted_space = CodingSpace(protein_sequence)

    print(len(str(unrestricted_space.count())))

Output:

.. code-block:: text

    148

The unrestricted space therefore contains a 148-digit number of synonymous
coding sequences. Now we exclude the two restriction sites:

.. code-block:: python

    from codeine import RestrictionSite
    from codeine.constraints import ForbiddenMotifs

    space = CodingSpace(
        protein_sequence,
        constraints=[
            ForbiddenMotifs([
                RestrictionSite.NdeI,
                RestrictionSite.PciI,
            ]),
        ],
    )

    print(space.count())

Output:

.. code-block:: text

    0

Despite the size of the unrestricted space, there is no synonymous coding
sequence for this protein that avoids both restriction sites.

Finding the conflict
--------------------

The incompatibility comes from just three consecutive amino acids: ``HMY``.
There are only four synonymous coding sequences for this peptide:

.. code-block:: python

    for sequence in sorted(CodingSpace('HMY')):
        print(sequence)

Output:

.. code-block:: text

    CACATGTAC
    CACATGTAT
    CATATGTAC
    CATATGTAT

Histidine can be encoded by ``CAC`` or ``CAT``, methionine only by ``ATG``,
and tyrosine by ``TAC`` or ``TAT``. The sequences beginning with ``CAC``
contain the ``PciI`` recognition sequence ``ACATGT``, while those beginning
with ``CAT`` contain the ``NdeI`` recognition sequence ``CATATG``.

We can confirm that this local peptide alone has an empty coding space under
the same constraints:

.. code-block:: python

    space = CodingSpace(
        'HMY',
        constraints=[
            ForbiddenMotifs([
                RestrictionSite.NdeI,
                RestrictionSite.PciI,
            ]),
        ],
    )

    print(space.count())

Output:

.. code-block:: text

    0

The same incompatibility also occurs for ``HMF`` and ``HM*``, where ``*``
denotes a stop codon. Phenylalanine codons and all stop codons begin with
``T``, just like tyrosine codons. Therefore, when histidine is encoded by
``CAC``, the first nucleotide after methionine completes the ``PciI`` motif
``ACATGT``; when histidine is encoded by ``CAT``, the ``NdeI`` motif
``CATATG`` is already present.

.. code-block:: python

    for peptide in ['HMY', 'HMF', 'HM*']:
        space = CodingSpace(
            peptide,
            constraints=[
                ForbiddenMotifs([
                    RestrictionSite.NdeI,
                    RestrictionSite.PciI,
                ]),
            ],
        )
        print(peptide, space.count())

Output:

.. code-block:: text

    HMY 0
    HMF 0
    HM* 0

Thus, ``HMY``, ``HMF``, and ``HM*`` cannot be encoded without introducing at
least one of these two restriction sites.

Combining constraints
---------------------

Empty spaces can also appear through interactions between different types of
constraints. Consider the peptide ``MDPP``:

.. code-block:: python

    unrestricted_space = CodingSpace('MDPP')
    print(unrestricted_space.count())

Output:

.. code-block:: text

    32

If we exclude the ``BamHI`` restriction site, half of these sequences remain:

.. code-block:: python

    space = CodingSpace(
        'MDPP',
        constraints=[
            ForbiddenMotifs([
                RestrictionSite.BamHI,
            ]),
        ],
    )

    print(space.count())

Output:

.. code-block:: text

    16

Methionine is encoded only by ``ATG``. When aspartate is encoded by ``GAT``,
the beginning of the first proline codon, ``CCN``, completes the ``BamHI``
recognition sequence ``GGATCC``. The 16 remaining sequences therefore all use
``GAC`` for aspartate.

Now add a maximum homopolymer length of four:

.. code-block:: python

    from codeine.constraints import MaxHomopolymer

    space.add_constraints(
        MaxHomopolymer(4)
    )

    print(space.count())

Output:

.. code-block:: text

    0

The remaining sequences contain the pattern ``GACCCNCCN``. Starting from the
final ``C`` of ``GAC``, the region ``CCCNCC`` contains five cytosines with at
most one interrupting nucleotide. Since ``MaxHomopolymer(4)`` considers a
single interruption by default, all 16 remaining sequences are excluded.

If only fully consecutive homopolymers are considered, the coding space is no
longer empty:

.. code-block:: python

    space = CodingSpace(
        'MDPP',
        constraints=[
            ForbiddenMotifs([
                RestrictionSite.BamHI,
            ]),
            MaxHomopolymer(4, allow_single_interruption=False),
        ],
    )

    print(space.count())

Output:

.. code-block:: text

    12

The same constraint interaction occurs for ``WDPP``. Tryptophan has a single
codon, ``TGG``, which, like the methionine codon ``ATG``, ends in ``G``.

This is not only a theoretical peptide. The serine/threonine protein kinase
UL97 from human cytomegalovirus strain AD169
(`UniProt P16788 <https://www.uniprot.org/uniprotkb/P16788/entry>`_) contains a
``WDPP`` segment, providing a biologically relevant example where this local
sequence can matter for coding-sequence design.

Thus, both ``MDPP`` and ``WDPP`` are incompatible with the combination of
``BamHI`` exclusion and ``MaxHomopolymer(4)`` under the default interrupted
homopolymer behavior.

This example shows how constraints that are individually satisfiable can become
incompatible when combined. Codeine represents the valid coding space directly,
so an empty space can be detected exactly and the local sequence responsible
for the conflict can be investigated explicitly.
