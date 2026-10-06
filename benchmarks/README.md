# Benchmarks

A broad performance suite for Codeine. Tests answer **is it correct?**; these benchmarks answer **did it get materially faster or slower, and where?**

## Running

```bash
python benchmarks/run.py
python benchmarks/run.py --suite full --label before-factorised
python benchmarks/run.py --case direct-repeat --operation compile --timeout 60
```

Each case/operation runs in a fresh subprocess with its own timeout, so one pathological compilation cannot stall the suite. Statuses are deliberately plain text: `PASS`, `TIMEOUT`, `ERROR`. The default is three repetitions and the terminal reports their median; raw timings are retained in JSON.

Results go into `benchmarks/results/`, which is gitignored. They are local measurements and shouldn't churn the repository. JSON records the Git revision, dirty-tree state, Python/platform, timeout and raw timings.

## Comparing runs

One run against another:

```bash
python benchmarks/compare.py benchmarks/results/before.json --against benchmarks/results/after.json
```

For a more trustworthy comparison, run each revision several times and compare groups. Each side can contain any number of files; matching raw timings are pooled and compared by median:

```bash
python benchmarks/compare.py benchmarks/results/before-*.json \
    --against benchmarks/results/after-*.json
```

This means there is no special "merge" step. A single run is convenient; three to five independent runs are better when a change is close or noisy. Changes within ±5% are displayed as `SAME`; that is only a readability convention, not a CI threshold.

A compact history is also available:

```bash
python benchmarks/compare.py --history benchmarks/results/*.json
```

## Coverage

The full suite includes controlled synthetic length scaling, real proteins from the test corpus, a large Cas9 input, codon weights, motifs, homopolymers, tandem/direct/inverted repeats, hairpins, combined constraints, mutation spaces, and deliberately repetitive/pathological proteins. Operations cover cold compilation/counting plus compiled sampling, containment, indexing and slicing.

`quick` is for frequent development checks. `full` is the comprehensive before/after suite for compiler work. Timing regressions are reported, not enforced in CI.

## Comparing compiler implementations

The cases stay fixed. If no compiler is specified, benchmarks use Codeine normally and call `compile()` with no compiler argument:

```bash
python benchmarks/run.py --suite full --label before
```

When Codeine exposes an alternative compiler, pass its name directly:

```bash
python benchmarks/run.py --suite full --compiler factorised --label factorised
```

The runner does not maintain its own compiler registry and does not interpret the value. It passes the string unchanged to Codeine as `compile(compiler="factorised")` and records `"compiler": "factorised"` in the result JSON. If `--compiler` is omitted, the JSON records `"compiler": null`.

An invalid compiler is Codeine's error to report. The worker does not catch or translate it; the benchmark is reported as `ERROR`, the traceback is retained in the result file, and the runner exits non-zero if any benchmark errors. This avoids duplicating Codeine's validation in the benchmark suite.
