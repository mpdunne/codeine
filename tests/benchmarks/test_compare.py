import json

import pytest

from benchmarks.compare import compare_results, comparison_status


def test_comparison_flags_count_mismatch():
    before = {'timings': [1.0], 'statuses': ['PASS'], 'counts': ['42']}
    after = {'timings': [0.5], 'statuses': ['PASS'], 'counts': ['43']}

    assert comparison_status(before, after)[0] == 'COUNT MISMATCH'


def test_comparison_rejects_different_timing_definitions(tmp_path):
    paths = []

    for version in [1, 2]:
        path = tmp_path / f'{version}.json'
        path.write_text(json.dumps({'schema_version': version, 'results': []}))
        paths.append(str(path))

    with pytest.raises(ValueError, match='timing definitions differ'):
        compare_results(paths[:1], paths[1:])
