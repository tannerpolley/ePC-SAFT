from __future__ import annotations

import json
import subprocess
import sys


def test_reactive_regression_benchmark_cases_and_order():
    from epcsaft.benchmarks.reactive_regression import CASE_BUILDERS

    assert tuple(CASE_BUILDERS) == (
        "reactive_speciation_batch_tiny",
        "reactive_bubble_batch_tiny",
        "reactive_regression_objective_tiny",
        "reactive_regression_parameter_perturbation",
    )


def test_reactive_regression_benchmark_schema_for_one_case():
    from epcsaft.benchmarks.reactive_regression import run_reactive_regression_benchmarks

    payload = run_reactive_regression_benchmarks(warmup=1, repeat=1, case="reactive_regression_objective_tiny")
    case_payload = payload["cases"][0]
    required = (
        "case",
        "row_count",
        "parameter_count",
        "warmup",
        "repeat",
        "median_ns",
        "mean_ns",
        "p10_ns",
        "p90_ns",
        "min_ns",
        "max_ns",
        "success_count",
        "failure_count",
        "residual_count",
        "fallback_used",
        "cache_hits",
        "cache_misses",
        "speciation_solves",
        "bubble_solves",
        "density_solves",
        "activity_calls",
        "fugacity_calls",
        "fingerprint",
        "diagnostics_keys",
    )
    for field in required:
        assert field in case_payload
    assert case_payload["case"] == "reactive_regression_objective_tiny"
    assert case_payload["success_count"] >= 1


def test_reactive_regression_benchmark_baseline_merge(tmp_path):
    from epcsaft.benchmarks.reactive_regression import run_reactive_regression_benchmarks

    baseline_path = tmp_path / "baseline.json"
    baseline_path.write_text(
        json.dumps({"cases": [{"case": "reactive_bubble_batch_tiny", "median_ns": 1000}]}, indent=2),
        encoding="utf-8",
    )

    payload = run_reactive_regression_benchmarks(
        warmup=1,
        repeat=1,
        case="reactive_bubble_batch_tiny",
        baseline_json=baseline_path,
    )
    case_payload = payload["cases"][0]
    assert case_payload["case"] == "reactive_bubble_batch_tiny"
    assert case_payload["baseline_median_ns"] == 1000
    assert case_payload["speedup_vs_baseline"] > 0.0


def test_reactive_regression_benchmark_script_executes_and_writes_json(tmp_path):
    output_path = tmp_path / "reactive_regression.json"
    result = subprocess.run(
        [
            sys.executable,
            "scripts/benchmark_reactive_regression.py",
            "--case",
            "reactive_speciation_batch_tiny",
            "--warmup",
            "1",
            "--repeat",
            "1",
            "--json",
            str(output_path),
        ],
        capture_output=True,
        text=True,
        check=False,
    )

    assert result.returncode == 0, result.stdout + result.stderr
    assert "reactive_speciation_batch_tiny" in result.stdout
    assert output_path.exists()

    payload = json.loads(output_path.read_text(encoding="utf-8"))
    assert payload["warmup"] == 1
    assert payload["repeat"] == 1
    assert payload["cases"][0]["case"] == "reactive_speciation_batch_tiny"
    assert payload["cases"][0]["cache_hits"] >= 0
    assert payload["cases"][0]["success_count"] >= 1
