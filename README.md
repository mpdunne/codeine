# 🧬 Codeine

**Codeine** is a Python library for constrained reverse translation of protein sequences.

## Installation & documentation

**Codeine** is available on PyPI:

```bash
pip install codeine
```

Full documentation: **https://codeine.readthedocs.io/**.

## Quick start

```python
from codeine import CodingSpace

space = CodingSpace('MKTIIALSYIFCLVF')

print(space.n_valid_sequences)
print(space.sample())
```

## Overview

A protein sequence can typically be encoded by an enormous number of synonymous DNA/RNA coding sequences.

For many biotechnological applications, such as recombinant expression, we must choose a coding sequence while respecting practical constraints, for example:

- avoiding restriction enzyme sites,
- avoiding nucleotide homopolymers,
- avoiding repetitive sequences,
- fixing specific codons,
- mutating relative to a reference sequence.

Identifying valid coding sequences under such constraints quickly becomes challenging, especially for longer proteins.

**Codeine** exactly represents the complete set of coding sequences for a given protein under specified constraints. It enables efficient counting, sampling, enumeration and mutation library design without the need for sequence filtering or repair.

## Examples

Count valid sequences:

```python
from codeine import CodingSpace

space = CodingSpace('MKTLEFQNGSCPRYKKL')

print(space.n_valid_sequences)
```

Sample a single valid sequence:

```python
from codeine import CodingSpace

space = CodingSpace('MKTLEFQNGSCPRYKKL')

seq = space.sample()
print(seq)
```

Sample many:

```python
from codeine import CodingSpace

space = CodingSpace('MKTLEFQNGSCPRYKKL')

for seq in space.sample(n=5):
    print(seq)
```

Apply constraints:

```python
from codeine import CodingSpace, RestrictionSite
from codeine.constraints import ForbiddenMotifs, MaxHomopolymer

space = CodingSpace(
    'MKTLEFQNGSCPRYKKL',
    fixed_codons={
        2: 'AAG',
        16: 'AAG',
    },
    constraints=[
        ForbiddenMotifs([
            RestrictionSite.EcoRI,
            RestrictionSite.BamHI,
            'CTGCAG',
        ]),
        MaxHomopolymer(max_length=4),
    ],
)

print(space.n_valid_sequences)
print(space.sample())
```

Use custom codon weights to change the sampling distribution:

```python
from codeine import CodingSpace, CodonWeights

weights = CodonWeights.ecoli()

space = CodingSpace('MKTLEFQNGSCPRYKKL', codon_weights=weights)

print(space.sample())
```

Use alternative genetic codes and RNA:

```python
from codeine import CodingSpace, TranslationTable

table = TranslationTable(table_id=2, rna=True)

space = CodingSpace('MKTLEFQNGSCPRYKKL', translation_table=table)

print(space.sample())
```

Enumerate all sequences (recommended only for small spaces):

```python
from codeine import CodingSpace

space = CodingSpace('CYIQNCPLG')

for sequence in space:
    print(sequence)
```

Explore mutants of a chosen reference sequence:

```python
from codeine import CodingSpace

space = CodingSpace('MKTLEFQNGSCPRYKKL')

reference = space.sample()

mutants = space.mutants(
    reference,
    free_positions=range(5, 14),
    min_nts=2,
    max_nts=5,
)

print(mutants.sample())
```
