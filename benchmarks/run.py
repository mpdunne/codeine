#!/usr/bin/env python3
"""
Run Codeine's performance benchmark suite.
"""


import argparse
import datetime as dt
import json
import platform
import statistics
import subprocess
import sys
from pathlib import Path
from typing import List, Optional, Sequence

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from benchmarks.cases import Case, all_cases

OPERATIONS = (
    'compile',
    'count',
    'sample-1',
    'sample-100',
    'contains',
    'index',
    'slice-100',
)


def git(*args: str) -> Optional[str]:
    """
    Run a Git query in the repository root.

    Parameters
    ----------
    *args
        Arguments passed to the ``git`` executable.

    Returns
    -------
    Stripped command output, or ``None`` when Git is unavailable or the
    repository query fails.
    """
    try:
        return subprocess.check_output(
            ['git', *args],
            cwd=ROOT,
            text=True,
            stderr=subprocess.DEVNULL,
        ).strip()
    except (OSError, subprocess.CalledProcessError):
        return None


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


def positive_int(value: str) -> int:
    """
    Parse a positive integer for ``argparse``.

    Parameters
    ----------
    value
        Command-line value to parse.

    Returns
    -------
    Parsed positive integer.

    Raises
    ------
    argparse.ArgumentTypeError
        If the value is less than one.
    """
    parsed = int(value)
    if parsed < 1:
        raise argparse.ArgumentTypeError('must be at least 1')
    return parsed


def positive_float(value: str) -> float:
    """
    Parse a positive floating-point value for ``argparse``.

    Parameters
    ----------
    value
        Command-line value to parse.

    Returns
    -------
    Parsed positive float.

    Raises
    ------
    argparse.ArgumentTypeError
        If the value is not greater than zero.
    """
    parsed = float(value)
    if parsed <= 0:
        raise argparse.ArgumentTypeError('must be greater than 0')
    return parsed


def parse_args() -> argparse.Namespace:
    """
    Parse command-line arguments.

    Returns
    -------
    Parsed benchmark-runner options.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--suite', choices=('quick', 'full'), default='quick')
    parser.add_argument(
        '--case',
        action='append',
        default=[],
        help='run cases whose names contain this text; may be repeated',
    )
    parser.add_argument(
        '--operation',
        action='append',
        choices=OPERATIONS,
        default=[],
        help='operation to benchmark; may be repeated',
    )
    parser.add_argument('--repeat', type=positive_int, default=3)
    parser.add_argument('--timeout', type=positive_float, default=30)
    parser.add_argument('--label', help='human-readable label stored with the run')
    parser.add_argument(
        '--compiler',
        help=(
            'compiler name passed unchanged to Codeine; if omitted, '
            'compile() is called normally'
        ),
    )
    parser.add_argument('--output', type=Path, help='result JSON path')
    return parser.parse_args()


def select_cases(suite: str, queries: Sequence[str]) -> List[Case]:
    """
    Select benchmark cases for a run.

    Parameters
    ----------
    suite
        Suite name that each selected case must belong to.
    queries
        Optional name fragments used to filter cases.

    Returns
    -------
    Matching cases in stable catalogue order.
    """
    selected = []

    for case in all_cases():
        if suite == 'quick' and not case.quick:
            continue
        if queries and not any(query in case.name for query in queries):
            continue
        selected.append(case)

    return selected


def worker_command(
    case: Case,
    operation: str,
    repeat: int,
    compiler: Optional[str],
) -> List[str]:
    """
    Build the subprocess command for one benchmark.

    Parameters
    ----------
    case
        Benchmark case to run.
    operation
        Operation to measure.
    repeat
        Number of timing repetitions.
    compiler
        Optional compiler name to pass through to Codeine.

    Returns
    -------
    Command-line arguments for the isolated benchmark worker.
    """
    command = [
        sys.executable,
        '-m',
        'benchmarks.worker',
        case.name,
        operation,
        '--repeat',
        str(repeat),
    ]

    if compiler is not None:
        command.extend(['--compiler', compiler])

    return command


def benchmark_result(
    case: Case,
    operation: str,
    status: str,
    timings: Optional[List[float]] = None,
    error: Optional[str] = None,
) -> dict:
    """
    Build a serialisable benchmark result record.

    Parameters
    ----------
    case
        Benchmark case that was run.
    operation
        Operation that was measured.
    status
        Final benchmark status.
    timings
        Successful timing measurements in seconds.
    error
        Captured error output, if the benchmark failed.

    Returns
    -------
    Result record suitable for inclusion in the output JSON document.
    """
    timings = timings or []
    median = statistics.median(timings) if timings else None

    return {
        'case': case.name,
        'group': case.group,
        'operation': operation,
        'status': status,
        'timings_s': timings,
        'median_s': median,
        'error': error,
    }


def run_benchmark(
    case: Case,
    operation: str,
    repeat: int,
    timeout: float,
    compiler: Optional[str],
) -> dict:
    """
    Run one benchmark in an isolated subprocess.

    Parameters
    ----------
    case
        Benchmark case to run.
    operation
        Operation to measure.
    repeat
        Number of timing repetitions.
    timeout
        Maximum subprocess runtime in seconds.
    compiler
        Optional compiler name to pass through to Codeine.

    Returns
    -------
    A serialisable PASS, TIMEOUT, or ERROR result record.
    """
    command = worker_command(case, operation, repeat, compiler)

    try:
        completed = subprocess.run(
            command,
            cwd=ROOT,
            text=True,
            capture_output=True,
            timeout=timeout,
        )
    except subprocess.TimeoutExpired:
        return benchmark_result(case, operation, 'TIMEOUT')

    if completed.returncode:
        output = completed.stderr or completed.stdout
        error = output[-4000:].strip()
        return benchmark_result(case, operation, 'ERROR', error=error)

    timings = json.loads(completed.stdout)['timings_s']
    return benchmark_result(case, operation, 'PASS', timings=timings)


def default_output_path(revision: str, dirty: bool) -> Path:
    """
    Build the default path for a benchmark result file.

    Parameters
    ----------
    revision
        Git revision label for the current checkout.
    dirty
        Whether the checkout contains uncommitted changes.

    Returns
    -------
    Timestamped JSON path below ``benchmarks/results``.
    """
    timestamp = dt.datetime.now().astimezone().strftime('%Y%m%d-%H%M%S')
    suffix = '_dirty' if dirty else ''
    filename = f"{timestamp}_{revision}{suffix}.json"
    return ROOT / 'benchmarks' / 'results' / filename


def print_header(
    suite: str,
    revision: str,
    dirty: bool,
    compiler: Optional[str],
) -> None:
    """
    Print benchmark-run metadata and the table header.

    Parameters
    ----------
    suite
        Name of the benchmark suite being run.
    revision
        Git revision label for the current checkout.
    dirty
        Whether the checkout contains uncommitted changes.
    compiler
        Optional compiler selected for the run.
    """
    dirty_label = ' (dirty)' if dirty else ''
    print(f"Codeine benchmarks [{suite}]  revision {revision}{dirty_label}")

    if compiler is not None:
        print(f"Compiler: {compiler}")

    print()
    print(f"{'BENCHMARK':58} {'MEDIAN':>10}  STATUS")
    print('-' * 82)


def print_result(result: dict, timeout: float) -> None:
    """
    Print one completed benchmark result.

    Parameters
    ----------
    result
        Serialised benchmark result record.
    timeout
        Configured timeout in seconds, used to display timeouts.
    """
    name = f"{result['case']} :: {result['operation']}"

    if result['status'] == 'TIMEOUT':
        elapsed = f">{timeout:g} s"
    else:
        elapsed = format_time(result['median_s'])

    print(f"{name:58} {elapsed:>10}  {result['status']}")


def make_document(
    args: argparse.Namespace,
    revision: str,
    dirty: bool,
    results: list,
) -> dict:
    """
    Build the JSON document for a benchmark run.

    Parameters
    ----------
    args
        Parsed benchmark-runner options.
    revision
        Git revision label for the current checkout.
    dirty
        Whether the checkout contains uncommitted changes.
    results
        Completed benchmark result records.

    Returns
    -------
    Complete serialisable benchmark-run document.
    """
    return {
        'schema_version': 1,
        'timestamp': dt.datetime.now().astimezone().isoformat(),
        'label': args.label,
        'suite': args.suite,
        'revision': revision,
        'dirty': dirty,
        'python': platform.python_version(),
        'platform': platform.platform(),
        'repeat': args.repeat,
        'timeout_s': args.timeout,
        'compiler': args.compiler,
        'results': results,
    }


def main() -> int:
    """
    Run the benchmark command-line interface.

    Returns
    -------
    zero when every benchmark passes, otherwise one.
    """
    args = parse_args()
    cases = select_cases(args.suite, args.case)
    operations = args.operation or list(OPERATIONS)

    if not cases:
        print('No benchmark cases matched.', file=sys.stderr)
        return 2

    revision = git('rev-parse', '--short', 'HEAD') or 'unknown'
    dirty_output = git('status', '--porcelain')
    dirty = bool(dirty_output) if dirty_output is not None else False

    output = args.output or default_output_path(revision, dirty)
    output.parent.mkdir(parents=True, exist_ok=True)

    print_header(args.suite, revision, dirty, args.compiler)

    results = []
    for case in cases:
        for operation in operations:
            result = run_benchmark(
                case=case,
                operation=operation,
                repeat=args.repeat,
                timeout=args.timeout,
                compiler=args.compiler,
            )
            results.append(result)
            print_result(result, args.timeout)

    document = make_document(args, revision, dirty, results)
    output.write_text(json.dumps(document, indent=2) + '\n')

    passed = sum(result['status'] == 'PASS' for result in results)
    errors = sum(result['status'] == 'ERROR' for result in results)

    print(f"\nCompleted {passed}/{len(results)} benchmarks. Results: {output}")
    return 1 if errors else 0


if __name__ == '__main__':
    raise SystemExit(main())
