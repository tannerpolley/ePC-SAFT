"""Opt-in runtime profiling test for MIAC fit generation.

Run with:

    set ePCSAFT_RUN_PERF=1
    python -m pytest tests/profile/test_miac_profile.py -s
"""

from __future__ import annotations

import os

import pytest

from analyses.data_validation.miac_fits.scripts.profile_miac_runtime import REPORT_CSV, REPORT_MD, run_miac_runtime_profile


def _should_run_perf() -> bool:
    return os.environ.get("ePCSAFT_RUN_PERF", "").strip().lower() in {"1", "true", "yes", "on"}


def test_runtime_profile_miac() -> None:
    if not _should_run_perf():
        pytest.skip("Set ePCSAFT_RUN_PERF=1 to run MIAC runtime profiling.")

    rows = run_miac_runtime_profile()
    assert rows, "MIAC runtime profile produced no rows."
    assert REPORT_CSV.exists(), f"Expected MIAC runtime profile CSV was not written: {REPORT_CSV}"
    assert REPORT_MD.exists(), f"Expected MIAC runtime profile Markdown was not written: {REPORT_MD}"
    assert any(str(row.get("category")) == "payload" for row in rows), "Expected at least one payload benchmark row."
    assert any(str(row.get("category")) == "point" for row in rows), "Expected at least one per-point benchmark row."
