"""Opt-in runtime profiling test for native regression benchmarking."""

from __future__ import annotations

import os

import pytest

from scripts.benchmarks.profile_regression_runtime import REPORT_CSV, REPORT_MD, run_regression_runtime_profile


def _should_run_perf() -> bool:
    return os.environ.get("ePCSAFT_RUN_PERF", "").strip().lower() in {"1", "true", "yes", "on"}


def test_runtime_profile_regression() -> None:
    if not _should_run_perf():
        pytest.skip("Set ePCSAFT_RUN_PERF=1 to run regression runtime profiling.")

    rows = run_regression_runtime_profile()
    assert rows, "Regression runtime profile produced no rows."
    assert REPORT_CSV.exists(), f"Expected regression runtime profile CSV was not written: {REPORT_CSV}"
    assert REPORT_MD.exists(), f"Expected regression runtime profile Markdown was not written: {REPORT_MD}"
    assert any(str(row.get("backend")) == "public_default" for row in rows)
    assert any(str(row.get("backend")) == "least_squares_native" for row in rows)
