Highly constrained spaces
=========================

Here, we use Codeine to show how strongly a set of common sequence design constraints can reduce a coding sequence space.

As an example, we use SpCas9, an RNA-guided DNA endonuclease from *Streptococcus pyogenes* that is widely used for CRISPR genome editing.

Initialising the coding space
-----------------------------

Let's grab the amino acid sequence for SpCas9, from
`UniProt Q99ZW2 <https://www.uniprot.org/uniprotkb/Q99ZW2/entry>`_.

.. code-block:: python

    spcas9 = (
        'MDKKYSIGLDIGTNSVGWAVITDEYKVPSKKFKVLGNTDRHSIKKNLIGALLFDSGETAEATRLKRTARR'
        'RYTRRKNRICYLQEIFSNEMAKVDDSFFHRLEESFLVEEDKKHERHPIFGNIVDEVAYHEKYPTIYHLRK'
        'KLVDSTDKADLRLIYLALAHMIKFRGHFLIEGDLNPDNSDVDKLFIQLVQTYNQLFEENPINASGVDAKA'
        'ILSARLSKSRRLENLIAQLPGEKKNGLFGNLIALSLGLTPNFKSNFDLAEDAKLQLSKDTYDDDLDNLLA'
        'QIGDQYADLFLAAKNLSDAILLSDILRVNTEITKAPLSASMIKRYDEHHQDLTLLKALVRQQLPEKYKEI'
        'FFDQSKNGYAGYIDGGASQEEFYKFIKPILEKMDGTEELLVKLNREDLLRKQRTFDNGSIPHQIHLGELH'
        'AILRRQEDFYPFLKDNREKIEKILTFRIPYYVGPLARGNSRFAWMTRKSEETITPWNFEEVVDKGASAQS'
        'FIERMTNFDKNLPNEKVLPKHSLLYEYFTVYNELTKVKYVTEGMRKPAFLSGEQKKAIVDLLFKTNRKVT'
        'VKQLKEDYFKKIECFDSVEISGVEDRFNASLGTYHDLLKIIKDKDFLDNEENEDILEDIVLTLTLFEDRE'
        'MIEERLKTYAHLFDDKVMKQLKRRRYTGWGRLSRKLINGIRDKQSGKTILDFLKSDGFANRNFMQLIHDD'
        'SLTFKEDIQKAQVSGQGDSLHEHIANLAGSPAIKKGILQTVKVVDELVKVMGRHKPENIVIEMARENQTT'
        'QKGQKNSRERMKRIEEGIKELGSQILKEHPVENTQLQNEKLYLYYLQNGRDMYVDQELDINRLSDYDVDH'
        'IVPQSFLKDDSIDNKVLTRSDKNRGKSDNVPSEEVVKKMKNYWRQLLNAKLITQRKFDNLTKAERGGLSE'
        'LDKAGFIKRQLVETRQITKHVAQILDSRMNTKYDENDKLIREVKVITLKSKLVSDFRKDFQFYKVREINN'
        'YHHAHDAYLNAVVGTALIKKYPKLESEFVYGDYKVYDVRKMIAKSEQEIGKATAKYFFYSNIMNFFKTEI'
        'TLANGEIRKRPLIETNGETGEIVWDKGRDFATVRKVLSMPQVNIVKKTEVQTGGFSKESILPKRNSDKLI'
        'ARKKDWDPKKYGGFDSPTVAYSVLVVAKVEKGKSKKLKSVKELLGITIMERSSFEKNPIDFLEAKGYKEV'
        'KKDLIIKLPKYSLFELENGRKRMLASAGELQKGNELALPSKYVNFLYLASHYEKLKGSPEDNEQKQLFVE'
        'QHKHYLDEIIEQISEFSKRVILADANLDKVLSAYNKHRDKPIREQAENIIHLFTLTNLGAPAAFKYFDTT'
        'IDRKRYTSTKEVLDATLIHQSITGLYETRIDLSQLGGD'
    )

Next, we initialise a ``CodingSpace`` for this sequence and count the number of possible synonymous coding sequences:

.. code-block:: python

    from codeine import CodingSpace
    from codeine.tools.display import format_count

    space = CodingSpace(spcas9)
    unconstrained_n = space.count()

    print(format_count(unconstrained_n))

Here we see that there are ``3.91 × 10^657`` sequences in this coding space. Which is a lot! What happens when we apply some constraints?

Applying constraints
--------------------

Let's now add several common sequence design constraints and see how each
affects the size of the coding space.

.. code-block:: python

    from codeine import RestrictionSite
    from codeine.constraints import (
        DirectRepeats,
        ForbiddenMotifs,
        InvertedRepeats,
        MaxHomopolymer,
        TandemRepeats,
    )


    def report_remaining(label):
        remaining = 100 * space.count() / unconstrained_n
        print(f'{remaining:.9f}% remaining after applying {label} constraints')


    # Restriction sites
    space.add_constraints(
        ForbiddenMotifs([
            RestrictionSite.BsaI,     # Golden Gate
            RestrictionSite.BsmBI,    # Golden Gate
            RestrictionSite.BbsI,     # Golden Gate
            RestrictionSite.EcoRI,    # Conventional cloning
            RestrictionSite.BamHI,    # Conventional cloning
            RestrictionSite.HindIII,  # Conventional cloning
        ])
    )
    report_remaining('restriction site')

    # Homopolymers
    space.add_constraints(
        MaxHomopolymer(6)
    )
    report_remaining('homopolymer')

    # Tandem repeats
    space.add_constraints([
        TandemRepeats(repeat_length=2, copies=4),
        TandemRepeats(repeat_length=3, copies=3),
        TandemRepeats(repeat_length=4, copies=3),
        TandemRepeats(repeat_length=5, copies=3),
        TandemRepeats(repeat_length=6, copies=3),
    ])
    report_remaining('tandem repeat')

Output:

.. code-block:: text

    0.000981271% remaining after applying restriction site constraints
    0.000014779% remaining after applying homopolymer constraints
    0.000003998% remaining after applying tandem repeat constraints

After applying these constraints, only 0.000003998% of the original coding sequence space remains:
approximately one in every 25 million sequences. The remaining space is still enormous, but valid
sequences are exceedingly sparse within the unconstrained space.

At this level of constraint, generating sequences and subsequently filtering or repairing
invalid designs becomes impractical. Codeine instead represents the valid sequence space directly,
allowing sequences to be counted, sampled and enumerated without subsequent filtering or repair.