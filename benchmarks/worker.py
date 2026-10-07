"""
Run one Codeine benchmark in an isolated process.
"""


import argparse
import json
import time

from benchmarks.cases import Case, all_cases


def parse_args() -> argparse.Namespace:
    """
    Parse worker command-line arguments.

    Returns
    -------
    Parsed case, operation, repeat count, and optional compiler name.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('case', help='benchmark case name')
    parser.add_argument('operation', help='operation to benchmark')
    parser.add_argument('--repeat', type=int, required=True)
    parser.add_argument('--compiler')
    return parser.parse_args()


def run_operation(space, operation: str) -> None:
    """
    Run a measured operation on an already compiled space.
    """
    if operation == 'count':
        space.count()
    elif operation == 'sample-1':
        space.sample()
    elif operation == 'sample-100':
        space.sample(100)
    elif operation == 'sample-10000':
        space.sample(10000)
    elif operation == 'contains':
        space.contains(space[0])
    elif operation == 'index':
        space[0]
    elif operation == 'slice-100':
        space[:100]
    else:
        raise ValueError(f'Unknown benchmark operation: {operation}')


def find_case(name: str) -> Case:
    """
    Find a benchmark case by exact name.

    Parameters
    ----------
    name
        Stable benchmark case name.

    Returns
    -------
    The matching benchmark case.

    Raises
    ------
    ValueError
        If no case has the requested name.
    """
    for case in all_cases():
        if case.name == name:
            return case
    raise ValueError(f"Unknown benchmark case: {name}")


def main() -> None:
    """
    Run one benchmark case and operation repeatedly in this process.
    """
    args = parse_args()
    case = find_case(args.case)
    timings = []

    counts = []

    for _ in range(args.repeat):
        start = time.perf_counter()
        space = case.build(compiler=args.compiler or 'flat')
        space.compile()
        compilation_time = time.perf_counter() - start
        counts.append(str(space.count()))

        if args.operation == 'compile':
            timings.append(compilation_time)
            continue

        # First-sample latency is separate from steady batch sampling.
        if args.operation in ('sample-100', 'sample-10000'):
            space.sample()

        sequence = space[0] if args.operation == 'contains' else None
        start = time.perf_counter()

        if args.operation == 'contains':
            for _ in range(1000):
                space.contains(sequence)
        else:
            run_operation(space, args.operation)

        timings.append(time.perf_counter() - start)

    print(json.dumps({'timings_s': timings, 'counts': counts}))


if __name__ == '__main__':
    main()
