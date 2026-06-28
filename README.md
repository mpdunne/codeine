# 🧬 Codeine

**Codeine** is a Python library for exploring and generating synonymous protein-coding sequences under biological constraints.

## Installation

```bash
pip install .
```

<!--## Installation & documentation

**Codeine** is available on PyPI:

```bash
pip install codeine
```

Full documentation: **https://codeine.readthedocs.io/**.-->

## Quick start

```python
from codeine import CodingSpace

space = CodingSpace('MKTIIALSYIFCLVF')

print(space.n_valid_sequences)
print(space.sample())
```


## Overview

A protein sequence can typically be encoded by an enormous number of synonymous DNA/RNA coding sequences.

For biotechnological applications such as recombinant expression, we typically must choose a coding sequence while respecting one or more practical constraints, for example:
- avoiding restriction enzyme sites,
- avoiding nucleotide homopolymers,
- fixing specific codons,
- mutating relative to a reference sequence.

Identifying valid coding sequences under such constraints quickly becomes challenging, especially for longer protein sequences.

**Codeine** represents the complete space of valid coding sequences for a given protein and experimental setup. It enables highly efficient sampling, enumeration and mutation library design while guaranteeing that the generated sequences satisfy user-specified constraints.

## Examples

Count valid sequences:

```python
from codeine import CodingSpace

space = CodingSpace('MKTLEFQNGSCPRYKKL')

print(space.n_valid_sequences)
```

Sample a single valid sequences:

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

Apply restrictions:

```python
from codeine import CodingSpace, RestrictionSite

space = CodingSpace(
    'MKTLEFQNGSCPRYKKL',
    forbidden_motifs=[
        RestrictionSite.EcoRI,
        RestrictionSite.BamHI,
        'CTGCAG',
    ],
    codon_restrictions={
        2: 'AAG',
        16: 'AAG',
    },
    max_homopolymer=4,
)

print(space.n_valid_sequences)
print(space.sample())
```

Use custom codon weights to change sampling distribution:

```python
from codeine import CodingSpace, CodonWeights

weights = CodonWeights.ecoli()

space = CodingSpace("MKTLEFQNGSCPRYKKL", codon_weights=weights)

seq = space.sample()
print(seq)
```

Nonstandard translation tables and RNA:

```python
from codeine import CodingSpace, TranslationTable

table = TranslationTable(table_id=2, rna=True)

space = CodingSpace("MKTLEFQNGSCPRYKKL", translation_table=table)

seq = space.sample()
print(seq)
```

Enumerate all sequences (only recommended for small spaces):

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

mutant = mutants.sample()
print(mutant)
```