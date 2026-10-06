"""
Run one Codeine benchmark in an isolated process.
"""


import argparse
import json
import time
from typing import Optional

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


def compile_space(space, compiler: Optional[str]) -> None:
    """
    Compile a coding space.

    Parameters
    ----------
    space
        Coding space to compile.
    compiler
        Optional compiler name passed unchanged to Codeine. When
        omitted, ``compile()`` is called without a compiler argument.
    """
    if compiler is None:
        space.compile()
    else:
        space.compile(compiler=compiler)


def run_operation(space, operation: str, compiler: Optional[str]) -> None:
    """
    Run one operation against a freshly built coding space.

    Parameters
    ----------
    space
        Coding space on which to run the operation.
    operation
        Benchmark operation name.
    compiler
        Optional compiler name passed through during compilation.

    Raises
    ------
    ValueError
        If ``operation`` is not recognised.
    """
    compile_space(space, compiler)

    if operation == 'compile':
        return
    if operation == 'count':
        space.count()
        return
    if operation == 'sample-1':
        space.sample()
        return
    if operation == 'sample-100':
        space.sample(100)
        return
    if operation == 'contains':
        sample = space.sample()
        space.contains(sample)
        return
    if operation == 'index':
        space[0]
        return
    if operation == 'slice-100':
        space[:100]
        return

    raise ValueError(f"Unknown benchmark operation: {operation}")


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

    for _ in range(args.repeat):
        start = time.perf_counter()
        space = case.build()
        run_operation(space, args.operation, args.compiler)
        timings.append(time.perf_counter() - start)

    print(json.dumps({'timings_s': timings}))


if __name__ == '__main__':
    main()
