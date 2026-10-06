#!/usr/bin/env python3
"""
Compare saved Codeine benchmark runs.
"""


import argparse
import json
import statistics
from pathlib import Path
from typing import Dict, List, Optional, Tuple

CHANGE_THRESHOLD_PERCENT = 5.0

BenchmarkKey = Tuple[str, str]
PooledResult = Dict[str, list]


def parse_args() -> argparse.Namespace:
    """
    Parse command-line arguments.

    Returns
    -------
    Parsed comparison options.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('files', nargs='+', help='baseline benchmark result files')
    parser.add_argument(
        '--against',
        nargs='+',
        help='candidate benchmark result files',
    )
    parser.add_argument(
        '--history',
        action='store_true',
        help='show a compact run history',
    )
    return parser.parse_args()


def read_result(path: str) -> dict:
    """
    Read a benchmark result document from disk.

    Parameters
    ----------
    path
        Path to a benchmark JSON file.

    Returns
    -------
    The decoded benchmark result document.
    """
    return json.loads(Path(path).read_text())


def empty_result() -> PooledResult:
    """
    Create an empty pooled benchmark result.

    Returns
    -------
    A pooled-result structure representing a missing benchmark.
    """
    return {'timings': [], 'statuses': []}


def load_results(paths: List[str]) -> Tuple[List[dict], Dict[BenchmarkKey, PooledResult]]:
    """
    Load benchmark runs and pool matching measurements.

    Parameters
    ----------
    paths
        Benchmark JSON files to load.

    Returns
    -------
    The decoded documents and results pooled by case and operation.
    """
    documents = []
    pooled = {}

    for path in paths:
        document = read_result(path)
        documents.append(document)

        for result in document['results']:
            key = (result['case'], result['operation'])
            pooled_result = pooled.setdefault(key, empty_result())
            pooled_result['timings'].extend(result.get('timings_s', []))
            pooled_result['statuses'].append(result['status'])

    return documents, pooled


def median_time(result: PooledResult) -> Optional[float]:
    """
    Calculate the median timing for a pooled result.

    Parameters
    ----------
    result
        Pooled benchmark result.

    Returns
    -------
    Median runtime in seconds, or ``None`` when no timing completed.
    """
    timings = result['timings']
    if not timings:
        return None
    return statistics.median(timings)


def result_status(result: PooledResult) -> str:
    """
    Choose the display status for a pooled result.

    Parameters
    ----------
    result
        Pooled benchmark result.

    Returns
    -------
    ``PASS`` when timings exist, otherwise the latest recorded status.
    """
    if result['timings']:
        return 'PASS'
    if result['statuses']:
        return result['statuses'][-1]
    return 'MISSING'


def format_time(seconds: Optional[float]) -> str:
    """
    Format a duration for terminal output.

    Parameters
    ----------
    seconds
        Duration in seconds, or ``None`` for a missing timing.

    Returns
    -------
    A compact human-readable duration.
    """
    if seconds is None:
        return '-'
    if seconds < 0.001:
        return f"{seconds * 1e6:.0f} us"
    if seconds < 1:
        return f"{seconds * 1e3:.1f} ms"
    return f"{seconds:.2f} s"


def classify_change(before: float, after: float) -> Tuple[str, float, float]:
    """
    Classify the change between two successful timings.

    Parameters
    ----------
    before
        Baseline runtime in seconds.
    after
        Candidate runtime in seconds.

    Returns
    -------
    The classification, percentage change, and candidate/baseline ratio.
    """
    ratio = after / before
    change_percent = (ratio - 1) * 100

    if change_percent < -CHANGE_THRESHOLD_PERCENT:
        status = 'FASTER'
    elif change_percent > CHANGE_THRESHOLD_PERCENT:
        status = 'SLOWER'
    else:
        status = 'SAME'

    return status, change_percent, ratio


def print_history(paths: List[str]) -> None:
    """
    Print a compact summary of saved benchmark runs.

    Parameters
    ----------
    paths
        Benchmark JSON files to include in the history table.
    """
    print(f"{'RUN':30} {'REV':10} {'PASS':>6} {'TIMEOUT':>8} {'ERROR':>7}")
    print('-' * 67)

    for path in paths:
        document = read_result(path)
        results = document['results']
        label = document.get('label') or Path(path).stem
        revision = document.get('revision', '?')

        passed = count_status(results, 'PASS')
        timed_out = count_status(results, 'TIMEOUT')
        errors = count_status(results, 'ERROR')

        print(
            f"{label[:30]:30} {revision:10} "
            f"{passed:6} {timed_out:8} {errors:7}"
        )


def count_status(results: List[dict], status: str) -> int:
    """
    Count benchmark results with a given status.

    Parameters
    ----------
    results
        Benchmark result records.
    status
        Status string to count.

    Returns
    -------
    Number of matching result records.
    """
    return sum(result['status'] == status for result in results)


def comparison_status(
    before_result: PooledResult,
    after_result: PooledResult,
) -> Tuple[str, str, Optional[float]]:
    """
    Describe the change between two pooled results.

    Parameters
    ----------
    before_result
        Pooled baseline result.
    after_result
        Pooled candidate result.

    Returns
    -------
    Display status, formatted percentage change, and timing ratio. The ratio
    is ``None`` when either side has no successful timing.
    """
    before = median_time(before_result)
    after = median_time(after_result)

    if before is not None and after is not None:
        status, change_percent, ratio = classify_change(before, after)
        return status, f"{change_percent:+.1f}%", ratio

    before_status = result_status(before_result)
    after_status = result_status(after_result)

    if before_status == after_status:
        return after_status, '-', None

    return f"{before_status}->{after_status}", '-', None


def print_comparison_row(
    key: BenchmarkKey,
    before_result: PooledResult,
    after_result: PooledResult,
    status: str,
    change: str,
) -> None:
    """
    Print one row of the comparison table.

    Parameters
    ----------
    key
        Case and operation identifying the benchmark.
    before_result
        Pooled baseline result.
    after_result
        Pooled candidate result.
    status
        Display classification for the comparison.
    change
        Formatted percentage change.
    """
    before = median_time(before_result)
    after = median_time(after_result)
    name = f"{key[0]} :: {key[1]}"

    print(
        f"{name[:54]:54} "
        f"{format_time(before):>9} "
        f"{format_time(after):>9} "
        f"{change:>9}  "
        f"{status}"
    )


def compare_results(baseline_paths: List[str], candidate_paths: List[str]) -> None:
    """
    Compare baseline and candidate benchmark runs.

    Parameters
    ----------
    baseline_paths
        Benchmark JSON files forming the baseline.
    candidate_paths
        Benchmark JSON files forming the candidate set.
    """
    baseline_documents, baseline = load_results(baseline_paths)
    candidate_documents, candidate = load_results(candidate_paths)

    print(
        f"Baseline: {len(baseline_documents)} run(s)  ->  "
        f"Candidate: {len(candidate_documents)} run(s)"
    )
    print(f"{'BENCHMARK':54} {'BEFORE':>9} {'AFTER':>9} {'CHANGE':>9}  STATUS")
    print('-' * 96)

    counts = {'FASTER': 0, 'SAME': 0, 'SLOWER': 0}
    transitions = 0
    ratios = []

    keys = sorted(set(baseline) | set(candidate))
    for key in keys:
        before_result = baseline.get(key, empty_result())
        after_result = candidate.get(key, empty_result())
        status, change, ratio = comparison_status(before_result, after_result)

        if status in counts:
            counts[status] += 1
        elif '->' in status:
            transitions += 1

        if ratio is not None:
            ratios.append(ratio)

        print_comparison_row(
            key,
            before_result,
            after_result,
            status,
            change,
        )

    print(
        f"\nSummary: {counts['FASTER']} faster, "
        f"{counts['SAME']} unchanged (±{CHANGE_THRESHOLD_PERCENT:g}%), "
        f"{counts['SLOWER']} slower, {transitions} status changes"
    )

    if ratios:
        geometric_mean = statistics.geometric_mean(ratios)
        overall_change = (geometric_mean - 1) * 100
        print(
            'Overall matched timing ratio: '
            f"{geometric_mean:.3f}x candidate/baseline "
            f"({overall_change:+.1f}%)"
        )


def main() -> None:
    """
    Run the benchmark-comparison command-line interface.
    """
    args = parse_args()

    if args.history:
        print_history(args.files)
        return

    if not args.against:
        raise SystemExit('compare.py: error: use --against to supply candidate run(s)')

    compare_results(args.files, args.against)


if __name__ == '__main__':
    main()
