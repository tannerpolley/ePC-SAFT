from __future__ import annotations

import epcsaft
from epcsaft.native_regression import solve_native_regression_residual_records


def _records() -> list[dict[str, object]]:
    return [
        {
            "row_id": "row_1",
            "family": "speciation",
            "target": "CO2",
            "predicted": 0.45,
            "observed": 0.40,
            "scale": 10.0,
        },
        {
            "row_id": "row_2",
            "family": "pressure",
            "target": "P",
            "predicted": 101400.0,
            "observed": 101325.0,
            "scale": 1.0e-5,
        },
    ]


def _parameters() -> list[dict[str, object]]:
    return [
        {
            "name": "CO2.k_ij",
            "path": "k_ij[CO2,H2O]",
            "kind": "binary_interaction",
            "initial": 0.01,
            "lower": -0.5,
            "upper": 0.5,
            "scale": 1.0,
        }
    ]


def test_native_regression_solve_rejects_unsupported_derivative_derivatives() -> None:
    result = solve_native_regression_residual_records(
        _records(),
        _parameters(),
        options={"derivative_backend": "unsupported_derivative"},
    )

    assert result["success"] is False
    assert result["status"] == "invalid_input"
    assert "unsupported_derivative is debug-only" in result["message"]


def test_native_regression_solve_accepts_analytic_derivative_policy() -> None:
    result = epcsaft.solve_native_regression_residual_records(
        _records(),
        _parameters(),
        options={"derivative_backend": "analytic", "optimizer_backend": "auto"},
    )

    assert result["success"] is True
    assert result["status"] == "converged"
    assert result["derivative_backend"] == "analytic"
    assert result["optimizer_backend"] == "analytic_linear_native"
    assert result["parameter_names"] == ["CO2.k_ij"]
    assert result["objective_result"]["fixed_shape_residuals"] is True
    assert result["objective_result"]["failure_count"] == 0


def test_native_regression_solve_reports_bounds_inconsistent() -> None:
    bad_parameters = _parameters()
    bad_parameters[0] = {**bad_parameters[0], "lower": 1.0, "upper": -1.0}

    result = solve_native_regression_residual_records(
        _records(),
        bad_parameters,
        options={"derivative_backend": "analytic"},
    )

    assert result["success"] is False
    assert result["status"] == "bounds_inconsistent"



