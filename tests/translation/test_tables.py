import pickle
import pytest

from codeine.translation.tables import TranslationTable


def test_default_codon_table_works_as_expected():
    tt = TranslationTable()

    assert tt.rna is False

    assert tt.codons_to_aa['ATG'] == 'M'
    assert tt.codons_to_aa['CCC'] == 'P'
    assert tt.codons_to_aa['GAC'] == 'D'

    assert tt.aa_to_codons['M'] == ('ATG',)
    assert set(tt.aa_to_codons['P']) == {'CCT', 'CCC', 'CCA', 'CCG'}
    assert set(tt.aa_to_codons['D']) == {'GAT', 'GAC'}


def test_rna_codon_table_works_as_expected():
    tt = TranslationTable(rna=True)

    assert tt.rna is True

    assert tt.codons_to_aa['AUG'] == 'M'
    assert tt.codons_to_aa['CCC'] == 'P'
    assert tt.codons_to_aa['GAC'] == 'D'

    assert 'ATG' not in tt.codons_to_aa

    assert tt.aa_to_codons['M'] == ('AUG',)
    assert set(tt.aa_to_codons['P']) == {'CCU', 'CCC', 'CCA', 'CCG'}
    assert set(tt.aa_to_codons['D']) == {'GAU', 'GAC'}


def test_unknown_translation_table_raises_value_error():
    with pytest.raises(ValueError, match='Unknown NCBI translation table'):
        TranslationTable(table_id=999999)

    with pytest.raises(ValueError, match='Unknown NCBI translation table'):
        TranslationTable(table_id='Michael')

    with pytest.raises(ValueError, match='Unknown NCBI translation table'):
        TranslationTable(table_id='')


def test_all_forward_codons_appear_in_reverse_table():
    tt = TranslationTable()
    reverse_codons = {
        codon
        for codons in tt.aa_to_codons.values()
        for codon in codons
    }
    assert reverse_codons == set(tt.codons_to_aa)


def test_all_forward_rna_codons_appear_in_reverse_table():
    tt = TranslationTable(rna=True)
    reverse_codons = {
        codon
        for codons in tt.aa_to_codons.values()
        for codon in codons
    }
    assert reverse_codons == set(tt.codons_to_aa)


def test_dna_table_is_dna_only():
    tt = TranslationTable()

    for codon in tt.codons_to_aa:
        assert 'U' not in codon

    for codons in tt.aa_to_codons.values():
        for codon in codons:
            assert 'U' not in codon


def test_rna_table_is_rna_only():
    tt = TranslationTable(rna=True)

    for codon in tt.codons_to_aa:
        assert 'T' not in codon

    for codons in tt.aa_to_codons.values():
        for codon in codons:
            assert 'T' not in codon


def test_codon_tables_are_read_only():
    tt = TranslationTable()

    with pytest.raises(TypeError):
        tt.codons_to_aa['AAA'] = 'X'

    with pytest.raises(TypeError):
        tt.aa_to_codons['X'] = ('AAA',)

    assert isinstance(tt.aa_to_codons['M'], tuple)
    with pytest.raises(AttributeError):
        tt.aa_to_codons['M'].append('XXX')

    with pytest.raises(AttributeError):
        tt.aa_to_codons = {'M': ('AAA',)}

    with pytest.raises(AttributeError):
        tt.codons_to_aa = {'AAA': 'M'}

    with pytest.raises(AttributeError):
        tt.chicken = 'beef'

    with pytest.raises(AttributeError):
        tt.rna = False


def test_normalise_codon():
    tt = TranslationTable(rna=False)
    assert tt.normalise_sequence('aug') == 'ATG'
    assert tt.normalise_sequence('ATG') == 'ATG'
    assert tt.normalise_sequence('ATg') == 'ATG'
    assert tt.normalise_sequence('ccc') == 'CCC'
    assert tt.normalise_sequence('ggg') == 'GGG'

    tt = TranslationTable(rna=True)
    assert tt.normalise_sequence('aug') == 'AUG'
    assert tt.normalise_sequence('ATG') == 'AUG'
    assert tt.normalise_sequence('ATg') == 'AUG'
    assert tt.normalise_sequence('ccc') == 'CCC'
    assert tt.normalise_sequence('ggg') == 'GGG'


def test_getitem_returns_amino_acid():
    tt = TranslationTable()
    assert tt['ATG'] == 'M'
    assert tt['CCC'] == 'P'
    assert tt['GAC'] == 'D'


def test_getitem_bad_key_raises():
    tt = TranslationTable()

    with pytest.raises(KeyError):
        _ = tt['x']

    with pytest.raises(KeyError):
        _ = tt['AT']

    with pytest.raises(KeyError):
        _ = tt['BAT']

    with pytest.raises(KeyError):
        _ = tt['M']

    with pytest.raises(KeyError):
        _ = tt['atgatg']

    with pytest.raises(KeyError):
        _ = tt[0]


def test_getitem_normalises_before_lookup():
    tt = TranslationTable()
    assert tt['aug'] == 'M'


def test_translation_table_pickle():
    table = TranslationTable(table_id=1, rna=False)
    loaded = pickle.loads(pickle.dumps(table))
    assert type(loaded) is TranslationTable
    assert loaded.table_id == table.table_id
    assert loaded.name == table.name
    assert loaded.rna == table.rna
    assert loaded.codons_to_aa == table.codons_to_aa
    assert loaded.aa_to_codons == table.aa_to_codons
    assert loaded['ATG'] == 'M'


def test_translation_table_pickle_rna():
    table = TranslationTable(table_id=1, rna=True)
    loaded = pickle.loads(pickle.dumps(table))
    assert type(loaded) is TranslationTable
    assert loaded.rna is True
    assert loaded['AUG'] == 'M'
    assert 'AUG' in loaded.codons_to_aa
    assert 'ATG' not in loaded.codons_to_aa


def test_translation_table_pickle_preserves_immutability():
    table = TranslationTable()
    loaded = pickle.loads(pickle.dumps(table))
    with pytest.raises(AttributeError):
        loaded.rna = True


def test_normalise_sequence_bad_inputs():
    tt = TranslationTable()
    with pytest.raises(ValueError):
        tt.normalise_sequence('xyz')


def test_non_standard_translation_tables_differ():
    standard = TranslationTable(table_id=1)
    mitochondrial = TranslationTable(table_id=2)

    assert mitochondrial.table_id == 2
    assert standard.codons_to_aa != mitochondrial.codons_to_aa


def test_translate_dna_sequence():
    tt = TranslationTable()
    assert tt.translate('ATGGAA') == 'ME'


def test_translate_lowercase_sequence():
    tt = TranslationTable()
    assert tt.translate('atggaa') == 'ME'


def test_translate_rna_sequence():
    tt = TranslationTable()
    assert tt.translate('AUGGAA') == 'ME'


def test_translate_empty_sequence():
    tt = TranslationTable()
    assert tt.translate('') == ''


def test_translate_rejects_unknown_codon():
    tt = TranslationTable()

    with pytest.raises(ValueError):
        tt.translate('ATGNAA')


def test_translate_rejects_sequences_of_bad_length():
    tt = TranslationTable()

    with pytest.raises(ValueError):
        tt.translate('A')

    with pytest.raises(ValueError):
        tt.translate('AT')

    with pytest.raises(ValueError):
        tt.translate('ATGG')


def test_custom_translation_table_translates_sequence():
    table = TranslationTable.custom({
        'ATG': 'M',
        'TAA': '*',
    })

    assert table.translate('ATGTAA') == 'M*'


def test_custom_translation_table_accepts_rna():
    table = TranslationTable.custom({
        'AUG': 'M',
        'UAA': '*',
    }, rna=True)

    assert table['AUG'] == 'M'
    assert table.translate('AUGUAA') == 'M*'
    assert table.codons_to_aa == {
        'AUG': 'M',
        'UAA': '*',
    }


def test_custom_translation_table_normalises_to_dna_by_default():
    table = TranslationTable.custom({
        'aug': 'M',
        'uaa': '*',
    })

    assert table['ATG'] == 'M'
    assert table['TAA'] == '*'
    assert table.translate('AUGUAA') == 'M*'


def test_custom_translation_table_normalises_to_rna():
    table = TranslationTable.custom({
        'ATG': 'M',
        'TAA': '*',
    }, rna=True)

    assert table['AUG'] == 'M'
    assert table['UAA'] == '*'
    assert table.translate('AUGUAA') == 'M*'


def test_custom_translation_table_sets_metadata():
    table = TranslationTable.custom(
        {'ATG': 'M', 'TAA': '*'},
        name='Nice table',
    )

    assert table.name == 'Nice table'
    assert table.names == ('Nice table',)
    assert table.table_id is None
    assert table.rna is False


def test_custom_translation_table_builds_reverse_mapping():
    table = TranslationTable.custom({
        'TTT': 'F',
        'TTC': 'F',
        'ATG': 'M',
        'TAA': '*',
    })

    assert table.aa_to_codons['F'] == ('TTT', 'TTC')
    assert table.aa_to_codons['M'] == ('ATG',)
    assert table.aa_to_codons['*'] == ('TAA',)


def test_custom_translation_table_sets_stop_codons():
    table = TranslationTable.custom({
        'ATG': 'M',
        'TAA': '*',
        'TAG': '*',
    })

    assert table.stop_codons == ('TAA', 'TAG')


def test_custom_translation_table_allows_nonstandard_amino_acids():
    table = TranslationTable.custom({
        'TGC': 'U',
        'TAG': 'X',
        'TAA': '*',
    })

    assert table['TGC'] == 'U'
    assert table['TAG'] == 'X'
    assert table['TAA'] == '*'


def test_custom_translation_table_allows_partial_tables():
    table = TranslationTable.custom({
        'ATG': 'M',
    })

    assert table.translate('ATG') == 'M'

    with pytest.raises(KeyError):
        table.translate('ATGTAA')


def test_custom_translation_table_is_immutable():
    table = TranslationTable.custom({'ATG': 'M'})

    with pytest.raises(AttributeError):
        table.name = 'Changed'


def test_custom_translation_table_rejects_invalid_dna_codons():
    with pytest.raises(ValueError, match='nucleotide sequence'):
        TranslationTable.custom({
            'ATX': 'M',
        })


def test_custom_translation_table_rejects_invalid_rna_codons():
    with pytest.raises(ValueError, match='nucleotide sequence'):
        TranslationTable.custom({
            'AUX': 'M',
        }, rna=True)


def test_custom_translation_table_rejects_non_triplet_codons():
    with pytest.raises(ValueError, match='length 3'):
        TranslationTable.custom({
            'AT': 'M',
        })


def test_custom_translation_table_rejects_empty_codons():
    with pytest.raises(ValueError, match='length 3'):
        TranslationTable.custom({
            '': 'M',
        })


def test_custom_translation_table_getitem_normalises_input():
    table = TranslationTable.custom({
        'ATG': 'M',
    })

    assert table['aug'] == 'M'
