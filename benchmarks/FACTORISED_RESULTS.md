# Factorised compiler: validation and benchmarks

Measured locally on 7 October 2026, Python 3.14.6, macOS 26.6.2 (x86_64).

The optional engine is implemented, with flat remaining the default. The original
AND/OR strategy is preserved: independently reduced factors, component splitting,
cached counts, and a persistent weighted sampling plan. Sequence ordering matches
flat; random streams are reproducible within an engine, not across engines.

## Comparison with the earlier prototype

These are fresh runs of the local `codeine-factorised-weight-fix` source, explicitly
selecting its factorised engine. Each number is the median of three repetitions
using the same interpreter and machine. Construction is included in compile time;
sampling is measured after first-use preparation. These are local measurements,
with visible run-to-run variation, rather than general performance guarantees.

| Profile | Prototype compile (ms) | New compile (ms) | Prototype 10,000 samples (ms) | New 10,000 samples (ms) |
| --- | ---: | ---: | ---: | ---: |
| MIKEYSASSAFRASMIKEYSASSAFRAS / repeat 15 | 542.35 | 616.58 | 297.30 | 240.91 |
| Same, E. coli weights | 520.79 | 587.84 | 263.05 | 206.84 |
| MIKEYMIKEY / repeat 9 | 4.10 | 4.61 | 89.27 | 76.62 |
| MIKEYAAAAAMIKEY-9 | 4.74 | 5.34 | 163.59 | 152.57 |
| MIKEYAAAAAMIKEY-15 | 5.20 | 4.79 | 214.23 | 163.02 |
| RNYKQT-4 | 4.66 | 7.92 | 117.14 | 103.08 |
| IHERQW-4 | 4.51 | 3.37 | 115.09 | 99.26 |

All seven profiles reproduce the prototype counts. The main repeat-15 case has
**542,376,004,976,640** valid sequences, with or without sampling weights. The new
implementation samples faster in these paired runs. Compilation is mixed: the main
case is about 13–14% slower than the prototype, while some small cases improve.

First-sample preparation is a separate cost: on the main case it changed from
31.18 ms to 65.71 ms uniformly, and from 38.98 ms to 60.98 ms with E. coli weights.
Containment for 1,000 checks improved from 633.25 ms to 20.43 ms uniformly.

The prototype's `auto` mode used flat when there were no factors. Therefore its
unconstrained timings cannot be compared as if they were factorised timings.

### Original-chat coverage

The [shared benchmark discussion](https://chatgpt.com/share/6a6f1207-19d0-83eb-94cd-d0221db8bdc5)
explicitly names MIKEYMIKEY/repeat 9, MIKEYAAAAAMIKEY/repeat 9, and the
MIKEY/SASSAFRAS case. It also lists a “longer MIKEY example” with repeat 15 but
omits its exact sequence. The local prototype tests contain MIKEYAAAAAMIKEY/repeat
15, RNYKQT/repeat 4, and IHERQW/repeat 4; these are now included too. The unnamed
chat entry cannot be claimed as an exact reproduction. The
[other shared chat](https://chatgpt.com/share/6ac4ca99-5e4c-83ed-8d36-9423518f2ee2)
does not provide further concrete benchmark definitions.

## Comparison with the current flat compiler

The original 25-case catalogue was run across eight operations, with three
repetitions and a five-second timeout for each worker. Flat completed 138/200 rows;
factorised completed 131/200. There were no ERROR results and all 131 successful
paired counts agreed. Timeouts include setup and all repetitions: they are not
lower bounds on a single operation's latency.

The artificial group was subsequently run with a ten-second limit (15 seconds for
the recovered cases). Flat timed out on the main hard case, weighted and uniform;
factorised completed both. All other artificial cases completed in both engines.

Some representative real workloads from the full sweep:

| Workload | Flat compile (ms) | Factorised compile (ms) | Flat 100 samples (ms) | Factorised 100 samples (ms) |
| --- | ---: | ---: | ---: | ---: |
| baseline/gfp | 6.48 | 16.75 | 8.09 | 60.57 |
| full/sfgfp/gene-synthesis | 259.63 | 87.61 | 11.52 | 10.58 |
| full/caplacizumab/gene-synthesis | 176.70 | 113.02 | 7.14 | 10.22 |
| full/spcas9/basic | 445.40 | 245.34 | 71.70 | 83.48 |

Factorised indexing is materially slower. For example, first-index lookup in
SpCas9/basic took 441 ms versus 0.508 ms for flat; sfGFP/gene-synthesis took
27.9 ms versus 0.136 ms. Indexing and enumeration are functionally complete, but
remain a performance optimisation opportunity. Ordinary sampling can also be
slower than flat. This supports retaining two engines, not replacing the default.

Elastin's repeat cases timed out in both engines. So did the stricter luciferase
repeat/hairpin cases and the pathological catalogue case at this limit.

Mutation benchmarks were corrected to use the same first valid reference sequence
for both engines: seeded samples can differ between engines. Their original full-
sweep timings should not be used as an identical-workload comparison. With the
corrected reference and a 30-second worker limit, all four rows completed in both
engines and their counts agreed. The sfGFP synthesis/distance case compiled in
2.69 s (flat) and 2.70 s (factorised), with 100 prepared samples taking 18.7 ms and
15.6 ms respectively.

## Reproducing and inspecting results

```bash
python benchmarks/run.py --suite full --compiler flat --timeout 30 --output benchmarks/results/flat.json
python benchmarks/run.py --suite full --compiler factorised --timeout 30 --output benchmarks/results/factorised.json
python benchmarks/compare.py benchmarks/results/flat.json --against benchmarks/results/factorised.json
python benchmarks/run.py --case artificial --compiler factorised --timeout 30
```

The initial full sweep used revisions ff22fad (flat) and c115c7b (factorised).
Follow-ups used 9464ec5 and 42acc04. Raw local results are gitignored:

- `flat-full.json` and `factorised-full.json`
- `artificial-flat.json` and `artificial-factorised.json`
- `recovered-flat.json` and `recovered-factorised.json`
- `mutation-flat.json` and `mutation-factorised.json`
- `prototype-comparison.json` (all seven direct prototype comparisons)
- `prototype-check.py` (the local reproduction driver; requires the prototype source)

All reside in `benchmarks/results/`. Schema 2 separates setup from measured
operations; comparisons reject mixing it with historical schema 1 results.

## Validation and limits

- 1,137 broad targeted tests passed (1,706 deselected; one existing collection warning).
- A separate 91-test run covered compiler integration, weighted sampling, mutation
  rebuilding, copies and benchmark comparison checks. These runs overlap; the
  numbers are not a combined distinct-test total.
- The recovered prototype cases were additionally checked against flat for exact
  counts, ordered slices, last indices and sampled membership.
- Lint for syntax/undefined names and `git diff --check` passed.
- Python 3.8 mixed-constraint compilation, sampling, slices, mutation spaces and
  pickle round-trips passed.
- Long-protein tests cover a 1,100-position count-constrained space, including
  sampling and pickling a prepared sampler without changing the recursion limit.

The full exhaustive/statistical test matrix was not run. Factorised pins and
constraint changes rebuild the model; weight changes reuse the counter. Custom
constraints need factor support. Automatic compiler selection remains deferred.
