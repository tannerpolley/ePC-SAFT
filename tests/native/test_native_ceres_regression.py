from __future__ import annotations

import pytest

from epcsaft.native_regression import solve_native_regression_residual_records


def test_native_regression_backend_reports_fixed_shape_objective_cost() -> None:
    result = solve_native_regression_residual_records(
        [
            {
                "row_id": "r1",
                "family": "activity",
                "target": "MEA",
                "predicted": 1.25,
                "observed": 1.00,
                "scale": 2.0,
            }
        ],
        [
            {
                "name": "MEA.d_born",
                "path": "d_born[MEA+]",
                "kind": "born_radius",
                "initial": 3.0,
                "lower": 2.0,
                "upper": 5.0,
            }
        ],
        options={"optimizer_backend": "auto", "derivative_backend": "analytic"},
    )

    assert result["success"] is True
    assert result["status"] == "converged"
    assert result["final_cost"] == result["initial_cost"]
    assert result["objective_result"]["residuals"] == [0.5]
    assert result["objective_result"]["residual_schema"][0]["family"] == "activity"
    assert result["objective_result"]["row_diagnostics"][0]["status"] == "converged"


def test_native_regression_backend_uses_analytic_sensitivities_to_update_parameters() -> None:
    result = solve_native_regression_residual_records(
        [
            {
                "row_id": "r1",
                "family": "activity",
                "target": "MEA",
                "predicted": 1.20,
                "observed": 1.00,
                "scale": 1.0,
                "sensitivities": {"MEA.d_born": 0.5},
            }
        ],
        [
            {
                "name": "MEA.d_born",
                "path": "d_born[MEA+]",
                "kind": "born_radius",
                "initial": 3.0,
                "lower": 2.0,
                "upper": 5.0,
            }
        ],
        options={"optimizer_backend": "auto", "derivative_backend": "analytic"},
    )

    assert result["success"] is True
    assert result["status"] == "converged"
    assert result["iterations"] == 1
    assert result["function_evaluations"] == 2
    assert result["parameters"][0] == pytest.approx(2.6)
    assert result["final_cost"] < result["initial_cost"]
    assert result["objective_result"]["residual_norm"] < 1.0e-9


def test_native_regression_backend_marks_all_rows_failed() -> None:
    result = solve_native_regression_residual_records(
        [
            {
                "row_id": "r1",
                "family": "pressure",
                "target": "P",
                "predicted": 0.0,
                "observed": 0.0,
                "success": False,
                "recoverable_failure": True,
            }
        ],
        [
            {
                "name": "CO2.k_ij",
                "path": "k_ij[CO2,H2O]",
                "kind": "binary_interaction",
                "initial": 0.0,
                "lower": -1.0,
                "upper": 1.0,
            }
        ],
        options={"optimizer_backend": "auto", "derivative_backend": "analytic"},
    )

    assert result["success"] is False
    assert result["status"] == "all_rows_failed"
    assert result["objective_result"]["failure_count"] == 1
    assert result["objective_result"]["row_diagnostics"][0]["penalty_applied"] is True
