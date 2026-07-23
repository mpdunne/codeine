Translation tables
==================

By default, **Codeine** uses the standard genetic code. It can also be
used to generate sequences under nonstandard genetic codes.

Choosing a translation table
----------------------------

The ``TranslationTable`` object governs how **Codeine** translates and
reverse-translates amino acid sequences. A table can be chosen from the
`NCBI list of genetic code tables <https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi>`_ by specifying
its ID:

.. code-block:: python

   from codeine import TranslationTable

   table1 = TranslationTable(table_id=1)

   print(table1.name)
   print('M', table1.aa_to_codons['M'])

   table2 = TranslationTable(table_id=2)

   print(table2.name)
   print('M', table2.aa_to_codons['M'])

Output:

.. code-block:: text

    Standard
    M ('ATG',)
    Vertebrate Mitochondrial
    M ('ATA', 'ATG')

To use a nonstandard genetic code with **Codeine**, simply pass one to the ``CodingSpace``:

.. code-block:: python

   from codeine import CodingSpace, TranslationTable

   table = TranslationTable(table_id=2)

   space = CodingSpace('CYIQNCPLG', translation_table=table)

Custom tables
-------------

Codeine also supports user-defined custom tables:

.. code-block:: python

    from codeine import TranslationTable

    table = TranslationTable.custom(
        codons_to_aa={
            'AAA': 'A',
            'AAC': 'B',
            'AAG': 'C',
            'AAT': 'D',
            'ACA': 'E',
            'ACC': 'F',
            'ACG': 'G',
            'ACT': '*',
        },
    )

    print(table.codons_to_aa['AAG'])

RNA tables
----------

By default, translation tables use DNA codons. To use RNA codons instead, set
``rna=True``.

.. code-block:: python

   table = TranslationTable(table_id=1, rna=True)

   print(table.translate['M'])

Available tables
----------------

Codeine supports the following `NCBI genetic code tables <https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi>`_:

.. list-table::
   :header-rows: 1

   * - ID
     - Name
   * - 1
     - Standard
   * - 2
     - Vertebrate Mitochondrial
   * - 3
     - Yeast Mitochondrial
   * - 4
     - Mold Mitochondrial
   * - 5
     - Invertebrate Mitochondrial
   * - 6
     - Ciliate Nuclear
   * - 9
     - Echinoderm Mitochondrial
   * - 10
     - Euplotid Nuclear
   * - 11
     - Bacterial
   * - 12
     - Alternative Yeast Nuclear
   * - 13
     - Ascidian Mitochondrial
   * - 14
     - Alternative Flatworm Mitochondrial
   * - 15
     - Blepharisma Macronuclear
   * - 16
     - Chlorophycean Mitochondrial
   * - 21
     - Trematode Mitochondrial
   * - 22
     - Scenedesmus obliquus Mitochondrial
   * - 23
     - Thraustochytrium Mitochondrial
   * - 24
     - Pterobranchia Mitochondrial
   * - 25
     - Candidate Division SR1
   * - 26
     - Pachysolen tannophilus Nuclear
   * - 27
     - Karyorelict Nuclear
   * - 28
     - Condylostoma Nuclear
   * - 29
     - Mesodinium Nuclear
   * - 30
     - Peritrich Nuclear
   * - 31
     - Blastocrithidia Nuclear
   * - 32
     - Balanophoraceae Plastid
   * - 33
     - Cephalodiscidae Mitochondrial