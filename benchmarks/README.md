# Benchmarks

A small performance suite for Codeine. Tests answer **is it correct?**; these benchmarks answer **did it get materially faster or slower, and where?**

## Running

```bash
python benchmarks/run.py
python benchmarks/run.py --suite full --label before-factorised
python benchmarks/run.py --case direct-repeat --operation compile --timeout 60
```

Each case/operation runs in a fresh subprocess with its own timeout, so one pathological compilation cannot stall the suite. The default timeout is 10 seconds. Statuses are deliberately plain text: `PASS`, `TIMEOUT`, `ERROR`. The default is three repetitions and the terminal reports their median; raw timings are retained in JSON.

Results go into `benchmarks/results/`, which is gitignored. They are local measurements and shouldn't churn the repository. JSON records the Git revision, dirty-tree state, Python/platform, timeout and raw timings.

## Coverage

The suite is broad enough to expose different compiler behaviour without taking a Cartesian product of every protein and constraint. It contains:

- unconstrained baselines for ubiquitin, GFP, mCherry, luciferase, caplacizumab and SpCas9;
- every main constraint family in isolation at two useful stringency levels, using modest proteins where possible and repeat-rich proteins where needed;
- full practical constraint stacks on sfGFP, luciferase, caplacizumab and SpCas9;
- the sfGFP and SpCas9 documentation workloads within those full-stack cases;
- one deliberately awkward repeat-heavy stress case;
- two representative sfGFP mutation spaces.

The individual cases are intentionally not a protein × constraint matrix. GFP covers motif and ordinary homopolymer constraints; luciferase adds a stricter homopolymer case and provides a longer ordinary coding space for inverted-repeat and hairpin constraints; collagen and elastin provide real repetitive proteins for tandem and direct-repeat constraints. SpCas9 is retained as the deliberately long real-world case rather than being used as the default substrate. A small artificial-repeat group also reproduces the original factorisation
experiments: MIKEYMIKEY, MIKEYAAAAAMIKEY, and MIKEYSASSAFRASMIKEYSASSAFRAS,
including a weighted version of the last case.

The full practical stack is exercised on several proteins rather than only one: sfGFP, luciferase, caplacizumab and SpCas9. The SpCas9 case uses its documented constraint stack. The full benchmark suite also includes the heavier elastin direct-repeat case.

Operations cover compilation/counting plus compiled sampling, containment, indexing and slicing.

## Comparing runs

```bash
python benchmarks/compare.py benchmarks/results/before.json --against benchmarks/results/after.json
```

For a more trustworthy comparison, run each revision several times and compare groups. Each side can contain any number of files; matching raw timings are pooled and compared by median. Changes within ±5% are displayed as `SAME`; that is only a readability convention, not a CI threshold.

A compact history is also available:

```bash
python benchmarks/compare.py --history benchmarks/results/*.json
```

## Comparing compiler implementations

If no compiler is specified, benchmarks use Codeine normally and call `compile()` with no compiler argument:

```bash
python benchmarks/run.py --suite full --label before
```

Select the alternative compiler at construction using:

```bash
python benchmarks/run.py --suite full --compiler factorised --label factorised
```

The runner passes the string unchanged to Codeine and records it in the result JSON. Invalid compiler names are left for Codeine to reject.

## Timing definitions (schema 2)

`compile` includes construction and compilation. All other operations start from
an already compiled space: `sample-1` includes first-use sampling preparation,
while `sample-100` and `sample-10000` warm the sampling cache before timing the
batch. `contains` measures 1,000 checks of a valid sequence prepared outside the
timer. `index` and `slice-100` measure the first index and first 100 sequences.

Each worker's timeout includes setup and all repetitions, even when setup is
outside the timed operation. A timeout is therefore a workload timeout, not a
lower bound on that operation's steady-state latency. Counts are recorded as
strings and comparisons flag disagreements. Schema 1 timings included setup in
every operation; comparisons reject mixing the two definitions.

Mutation benchmarks use the first valid sequence as their reference, rather than
an engine-dependent random sample, so both engines receive identical workloads.
Run the artificial repeat cases with ``--case artificial``.
