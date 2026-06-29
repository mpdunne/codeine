import json
from pathlib import Path

from Bio.Data import CodonTable


out_file = Path('codeine/translation/data/tables.json')


def table_to_dict(table_id: int) -> dict:
    table = CodonTable.unambiguous_dna_by_id[table_id]

    codon_to_aa = dict(table.forward_table)

    for codon in table.stop_codons:
        codon_to_aa[codon] = '*'

    return {
        'id': table_id,
        'name': table.names[0],
        'names': list(table.names),
        'start_codons': list(table.start_codons),
        'stop_codons': list(table.stop_codons),
        'codon_to_aa': dict(sorted(codon_to_aa.items())),
    }


def main() -> None:
    tables = {
        str(table_id): table_to_dict(table_id)
        for table_id in sorted(CodonTable.unambiguous_dna_by_id)
    }

    out_file.parent.mkdir(parents=True, exist_ok=True)

    with out_file.open('w') as f:
        json.dump(tables, f, indent=2, sort_keys=True)

    print(f'Wrote {len(tables)} translation tables to {out_file}')


if __name__ == '__main__':
    main()
