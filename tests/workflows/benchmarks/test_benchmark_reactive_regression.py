from __future__ import annotations

import json
import subprocess
import sys


def test_reactive_regression_benchmark_cases_and_order():
    from scripts.benchmarks.helpers.reactive_regression import CASE_BUILDERS, DEFAULT_CASES

    assert tuple(CASE_BUILDERS) == (
        "reactive_speciation_batch_tiny",
        "reactive_bubble_batch_tiny",
        "reactive_regression_objective_tiny",
        "reactive_regression_parameter_perturbation",
        "reactive_regression_pressure_speciation_35_row_surrogate",
        "mea_trace_carbonate_35_row_public_surrogate",
    )
    assert DEFAULT_CASES == (
        "reactive_speciation_batch_tiny",
        "reactive_bubble_batch_tiny",
        "reactive_regression_objective_tiny",
        "reactive_regression_parameter_perturbation",
    )


def test_reactive_regression_benchmark_schema_for_one_case():
    from scripts.benchmarks.helpers.reactive_regression import run_reactive_regression_benchmarks

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
        "measured_success_repeat_count",
        "failure_messages",
        "baseline_repeat",
        "baseline_warmup",
        "residual_count",
        "cache_hits",
        "cache_misses",
        "context_cache_hits",
        "context_cache_misses",
        "objective_seed_hits",
        "objective_seed_misses",
        "native_reference_state_cache_hits",
        "native_reference_state_cache_misses",
        "density_warm_start_hits",
        "speciation_solves",
        "bubble_solves",
        "density_solves",
        "activity_calls",
        "fugacity_calls",
        "fingerprint",
        "diagnostics_keys",
        "target_family_counts",
    )
    for field in required:
        assert field in case_payload
    assert case_payload["case"] == "reactive_regression_objective_tiny"
    assert case_payload["success_count"] >= 1
    assert case_payload["measured_success_repeat_count"] == 1
    assert case_payload["failure_messages"] == []
    assert case_payload["baseline_repeat"] == 1
    assert case_payload["baseline_warmup"] == 1
    assert case_payload["context_cache_hits"] >= 0
    assert case_payload["native_reference_state_cache_hits"] is None


def test_reactive_regression_benchmark_baseline_merge(tmp_path):
    from scripts.benchmarks.helpers.reactive_regression import run_reactive_regression_benchmarks

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


def test_reactive_regression_benchmark_excludes_failed_repeats_from_timing():
    from scripts.benchmarks.helpers import reactive_regression as bench

    calls = {"count": 0}

    def flaky_runner():
        calls["count"] += 1
        if calls["count"] == 1:
            raise RuntimeError("synthetic benchmark failure")
        return bench.BenchmarkObservation(
            fingerprint={"case": "synthetic"},
            diagnostics={"diagnostics_keys": []},
            row_count=1,
            parameter_count=1,
            success_count=1,
            failure_count=0,
            residual_count=1,
            cache_hits=0,
            cache_misses=0,
            speciation_solves=1,
            bubble_solves=0,
            density_solves=0,
            activity_calls=1,
            fugacity_calls=0,
        )

    prepared = bench.PreparedBenchmarkCase(
        case="synthetic_failure_timing",
        description="synthetic failure timing case",
        runner=flaky_runner,
    )
    payload = bench._benchmark_case(prepared, warmup=0, repeat=2)

    assert payload["measured_success_repeat_count"] == 1
    assert payload["failure_count"] == 1
    assert payload["failure_messages"] == ["RuntimeError: synthetic benchmark failure"]
    assert payload["min_ns"] > 1


def test_reactive_regression_benchmark_raises_when_all_repeats_fail():
    from scripts.benchmarks.helpers import reactive_regression as bench

    prepared = bench.PreparedBenchmarkCase(
        case="synthetic_all_failed",
        description="synthetic all failed case",
        runner=lambda: (_ for _ in ()).throw(RuntimeError("all failed")),
    )

    import pytest

    with pytest.raises(RuntimeError, match="failed for every measured repeat"):
        bench._benchmark_case(prepared, warmup=0, repeat=2)


def test_reactive_regression_benchmark_has_35_row_public_surrogate():
    from scripts.benchmarks.helpers.reactive_regression import run_reactive_regression_benchmarks

    payload = run_reactive_regression_benchmarks(
        warmup=0,
        repeat=1,
        case="mea_trace_carbonate_35_row_public_surrogate",
    )
    case_payload = payload["cases"][0]

    assert case_payload["row_count"] == 35
    assert case_payload["measured_success_repeat_count"] == 1
    assert case_payload["speciation_solves"] >= 35
    assert case_payload["bubble_solves"] == 0
    assert "surrogate" in case_payload["case"].lower()
    assert "synthetic rows" in case_payload["fingerprint"]["surrogate_note"].lower()


def test_reactive_regression_benchmark_has_35_row_pressure_speciation_surrogate():
    from scripts.benchmarks.helpers.reactive_regression import run_reactive_regression_benchmarks

    payload = run_reactive_regression_benchmarks(
        warmup=0,
        repeat=1,
        case="reactive_regression_pressure_speciation_35_row_surrogate",
    )
    case_payload = payload["cases"][0]

    assert case_payload["row_count"] == 35
    assert case_payload["measured_success_repeat_count"] == 1
    assert case_payload["speciation_solves"] >= 33
    assert case_payload["bubble_solves"] >= 2
    assert case_payload["target_family_counts"]["partial_pressure"] >= 2
    assert case_payload["target_family_counts"]["speciation"] >= 35
    assert case_payload["target_family_counts"]["activity"] >= 33
    assert case_payload["fingerprint"]["target_family_counts"]["partial_pressure"] >= 2
    assert "pressure" in case_payload["fingerprint"]["surrogate_note"].lower()


def test_reactive_regression_benchmark_script_executes_and_writes_json(tmp_path):
    output_path = tmp_path / "reactive_regression.json"
    result = subprocess.run(
        [
            sys.executable,
            "scripts/benchmarks/benchmark_reactive_regression.py",
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
