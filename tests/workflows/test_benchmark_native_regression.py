from __future__ import annotations

import json
import subprocess
import sys


def test_native_regression_benchmark_cases_and_order():
    from epcsaft.benchmarks.native_regression import CASE_BUILDERS, DEFAULT_CASES

    assert tuple(CASE_BUILDERS) == (
        "native_neutral_density_tiny",
        "native_binary_kij_tiny",
        "native_reactive_born_ssmds_tiny",
        "native_mea_pressure_speciation_35_row_surrogate",
    )
    assert DEFAULT_CASES == tuple(CASE_BUILDERS)


def test_native_regression_benchmark_schema_and_reactive_coverage():
    from epcsaft.benchmarks.native_regression import run_native_regression_benchmarks

    payload = run_native_regression_benchmarks(
        warmup=0,
        repeat=1,
        case="native_reactive_born_ssmds_tiny",
    )
    case_payload = payload["cases"][0]

    assert case_payload["case"] == "native_reactive_born_ssmds_tiny"
    assert case_payload["success"] is True
    assert case_payload["status"] == "converged"
    assert case_payload["fixed_shape_residuals"] is True
    assert case_payload["production_finite_difference_allowed"] is False
    assert {"pressure", "speciation", "activity"} <= set(case_payload["target_families"])
    assert {"born_radius", "solvation_factor"} <= set(case_payload["parameter_kinds"])


def test_native_regression_benchmark_has_35_row_public_surrogate():
    from epcsaft.benchmarks.native_regression import run_native_regression_benchmarks

    payload = run_native_regression_benchmarks(
        warmup=0,
        repeat=1,
        case="native_mea_pressure_speciation_35_row_surrogate",
    )
    case_payload = payload["cases"][0]

    assert case_payload["row_count"] == 35
    assert case_payload["residual_count"] == 105
    assert {"pressure", "speciation"} <= set(case_payload["target_families"])
    assert {"born_radius", "solvation_factor", "reaction_equilibrium_constant"} <= set(case_payload["parameter_kinds"])


def test_native_regression_benchmark_baseline_merge(tmp_path):
    from epcsaft.benchmarks.native_regression import run_native_regression_benchmarks

    baseline_path = tmp_path / "baseline.json"
    baseline_path.write_text(
        json.dumps({"cases": [{"case": "native_binary_kij_tiny", "median_ns": 1000}]}, indent=2),
        encoding="utf-8",
    )

    payload = run_native_regression_benchmarks(
        warmup=0,
        repeat=1,
        case="native_binary_kij_tiny",
        baseline_json=baseline_path,
    )

    assert payload["cases"][0]["baseline_median_ns"] == 1000
    assert payload["cases"][0]["speedup_vs_baseline"] > 0.0


def test_native_regression_benchmark_script_executes_and_writes_json(tmp_path):
    output_path = tmp_path / "native_regression.json"
    result = subprocess.run(
        [
            sys.executable,
            "scripts/benchmark_native_regression.py",
            "--case",
            "native_neutral_density_tiny",
            "--warmup",
            "0",
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
    assert "native_neutral_density_tiny" in result.stdout
    assert output_path.exists()


def test_native_ceres_thermo_regression_benchmark_schema():
    from epcsaft.benchmarks.native_ceres_thermo_regression import run_native_ceres_thermo_regression_benchmark

    payload = run_native_ceres_thermo_regression_benchmark(warmup=0, repeat=1)

    assert payload["case"] == "reactive_speciation_logk_implicit"
    assert payload["optimizer_backend"] in {"backend_unavailable", "ceres"}
    assert payload["derivative_backend"] in {"implicit", "analytic_implicit"}
    assert payload["python_objective_used"] is False
    assert payload["finite_difference_used"] is False
    assert payload["initial_cost"] >= 0.0
    assert payload["final_cost"] >= 0.0


def test_native_ceres_thermo_regression_benchmark_script_executes_and_writes_json(tmp_path):
    output_path = tmp_path / "native_ceres_thermo_regression.json"
    result = subprocess.run(
        [
            sys.executable,
            "scripts/benchmark_native_ceres_thermo_regression.py",
            "--warmup",
            "0",
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
    assert "reactive_speciation_logk_implicit" in result.stdout
    payload = json.loads(output_path.read_text(encoding="utf-8"))
    assert payload["python_objective_used"] is False
    assert payload["finite_difference_used"] is False
