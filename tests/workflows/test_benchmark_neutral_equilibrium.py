from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path


def test_neutral_benchmark_module_exposes_required_cases():
    from epcsaft.benchmarks.neutral_equilibrium import CASE_BUILDERS

    assert tuple(CASE_BUILDERS) == (
        "neutral_state",
        "tp_flash",
        "bubble_p",
        "dew_p",
        "lle_seeded",
    )


def test_neutral_benchmark_module_applies_baseline_speedup(tmp_path):
    from epcsaft.benchmarks.neutral_equilibrium import run_neutral_equilibrium_benchmarks

    baseline_path = tmp_path / "baseline.json"
    baseline_path.write_text(
        json.dumps(
            {
                "cases": [
                    {
                        "case": "tp_flash",
                        "median_ns": 1000,
                    }
                ]
            }
        ),
        encoding="utf-8",
    )

    payload = run_neutral_equilibrium_benchmarks(
        warmup=1,
        repeat=2,
        case="tp_flash",
        baseline_json=baseline_path,
    )

    assert payload["cases"][0]["case"] == "tp_flash"
    assert payload["cases"][0]["baseline_median_ns"] == 1000
    assert payload["cases"][0]["speedup_vs_baseline"] > 0.0


def test_neutral_benchmark_script_runs_tp_flash_case_and_writes_json(tmp_path):
    output_path = tmp_path / "neutral_equilibrium.json"

    result = subprocess.run(
        [
            sys.executable,
            "scripts/benchmark_neutral_equilibrium.py",
            "--case",
            "tp_flash",
            "--warmup",
            "1",
            "--repeat",
            "2",
            "--json",
            str(output_path),
        ],
        capture_output=True,
        text=True,
        check=False,
    )

    assert result.returncode == 0, result.stdout + result.stderr
    assert "tp_flash" in result.stdout
    assert output_path.exists()

    payload = json.loads(output_path.read_text(encoding="utf-8"))
    assert payload["warmup"] == 1
    assert payload["repeat"] == 2
    assert payload["cases"][0]["case"] == "tp_flash"
    assert payload["cases"][0]["failures"] == 0
    assert payload["cases"][0]["fallback_used"] is False
    assert "median_ns" in payload["cases"][0]
    assert "fingerprint" in payload["cases"][0]
