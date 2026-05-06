from __future__ import annotations

import json
import importlib
import os
from pathlib import Path

import numpy as np
import pytest

from epcsaft.equilibrium_core.confidence import benchmark_case_to_native_inputs
from epcsaft.equilibrium_core.confidence import load_benchmark_suite
from epcsaft.equilibrium_core.confidence import run_confidence_suite
from epcsaft.equilibrium_core.confidence import run_smoke_cases


def test_khudaida_benchmark_fixture_loads_charge_neutral_cases() -> None:
    suite = load_benchmark_suite("khudaida_2026")

    assert suite.name == "khudaida_2026"
    assert suite.species == ("H2O", "Ethanol", "Butanol", "Na+", "Cl-")
    assert len(suite.cases) >= 35
    assert suite.thresholds["solver_residual_norm"] == pytest.approx(1.0e-6)

    for case in suite.cases:
        assert case.feed_formula.shape == (4,)
        assert case.experimental_organic_formula.shape == (4,)
        assert case.experimental_aqueous_formula.shape == (4,)
        assert np.isclose(np.sum(case.feed_formula), 1.0)
        assert np.isclose(np.sum(case.experimental_organic_formula), 1.0)
        assert np.isclose(np.sum(case.experimental_aqueous_formula), 1.0)
        native = benchmark_case_to_native_inputs(case)
        assert native.feed.shape == (5,)
        assert abs(float(native.feed[3] - native.feed[4])) <= 1.0e-12
        assert abs(float(native.initial_aq[3] - native.initial_aq[4])) <= 1.0e-12
        assert abs(float(native.initial_org[3] - native.initial_org[4])) <= 1.0e-12


def test_khudaida_smoke_cases_return_results_or_diagnostic_failures() -> None:
    suite = load_benchmark_suite("khudaida_2026")
    predictions = run_smoke_cases(suite, case_keys=("0.05:293.15:1",))

    assert len(predictions) == 1
    for prediction in predictions:
        diagnostics = prediction.diagnostics
        assert diagnostics["phase_equilibrium_model"] == "electrolyte_lle_v5_native_charge_constrained_solve"
        assert diagnostics["equilibrium_route"] == "electrolyte_lle"
        assert diagnostics["solver_language"] == "c++"
        assert "ceres" not in json.dumps(diagnostics).lower()
        assert prediction.status in {"accepted", "diagnostic_failure"}
        if prediction.status == "accepted":
            assert prediction.metrics is not None
            assert np.isfinite(prediction.metrics.grand_aad)
        else:
            assert diagnostics["acceptance_gate"] == "predictive_solve_failed"
            assert diagnostics["best_failure_reason"]


def test_confidence_suite_smoke_mode_writes_bounded_report(tmp_path: Path) -> None:
    report = run_confidence_suite("khudaida_2026", mode="smoke", output_root=tmp_path)
    summary = json.loads(report.summary_path.read_text(encoding="utf-8"))

    assert summary["mode"] == "smoke"
    assert summary["case_count"] == 2
    assert summary["oracle_rows"] == 1
    assert summary["stress_rows"] == 0
    assert summary["sensitivity_rows"] == 1


def test_khudaida_paper_validation_recompute_uses_native_lle() -> None:
    common = importlib.import_module("analyses.2026_khudaida.scripts._common")
    exp_row = common._experimental_rows(0.05, 293.15)[0]
    feed_row = common._digitized_feed_rows_for_figure(2, 293.15, 0.05)[0]
    result = common._solve_formula_feed(
        293.15,
        feed_row["feed_formula"],
        [exp_row["organic_formula"], exp_row["aqueous_formula"]],
    )

    assert result is not None
    assert result["source"] == "epcsaft_native_v5"
    assert result["converged"] is True
    assert result["residual_norm"] <= 1.0e-6
    assert np.all(np.isfinite(result["organic_formula"]))
    assert np.all(np.isfinite(result["aqueous_formula"]))


@pytest.mark.skipif(
    os.environ.get("EPCSAFT_EQUILIBRIUM_CONFIDENCE", "").lower() not in {"1", "true", "yes", "on"},
    reason="Set EPCSAFT_EQUILIBRIUM_CONFIDENCE=1 to run opt-in electrolyte LLE confidence checks.",
)
def test_opt_in_confidence_report_generates_full_outputs(tmp_path: Path) -> None:
    report = run_confidence_suite("khudaida_2026", mode="full", output_root=tmp_path)

    assert report.summary_path.exists()
    assert report.benchmark_csv.exists()
    assert report.continuation_csv.exists()
    assert report.oracle_csv.exists()
    assert report.sensitivity_csv.exists()
    assert report.residual_gate_plot.exists()
    assert report.error_plot.exists()
    assert report.continuation_plot.exists()
    assert report.sensitivity_plot.exists()

    summary = json.loads(report.summary_path.read_text(encoding="utf-8"))
    assert summary["suite"] == "khudaida_2026"
    assert summary["mode"] == "full"
    assert summary["case_count"] >= 35
    assert summary["accepted_count"] + summary["diagnostic_failure_count"] == summary["case_count"]
    assert summary["ceres_reported"] is False
