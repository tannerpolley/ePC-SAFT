from __future__ import annotations

import pytest

import epcsaft
from epcsaft import _core
from epcsaft.native_regression import (
    CANONICAL_NATIVE_REGRESSION_STATUSES,
    evaluate_native_regression_residual_records,
    native_regression_contract_schema,
)


def test_native_regression_contract_schema_is_exposed_from_core() -> None:
    schema = _core._native_regression_contract_schema()

    assert schema["statuses"] == [
        "converged",
        "max_iterations",
        "line_search_failed",
        "singular_jacobian",
        "all_rows_failed",
        "nonfinite_objective",
        "bounds_inconsistent",
        "invalid_input",
    ]
    assert schema["fixed_shape_residuals"] is True
    assert schema["production_finite_difference_allowed"] is False
    assert {"pressure", "speciation", "activity"}.issubset(schema["supported_target_families"])
    assert {"binary_interaction", "born_radius"}.issubset(schema["supported_parameter_kinds"])
    assert schema["row_diagnostic_fields"] == [
        "row_id",
        "success",
        "status",
        "message",
        "residual_start",
        "residual_count",
        "penalty_applied",
        "solve_backend",
        "derivative_backend",
    ]


def test_python_native_regression_contract_wrapper_matches_public_exports() -> None:
    schema = native_regression_contract_schema()

    assert tuple(schema["statuses"]) == CANONICAL_NATIVE_REGRESSION_STATUSES
    assert epcsaft.CANONICAL_NATIVE_REGRESSION_STATUSES == CANONICAL_NATIVE_REGRESSION_STATUSES
    assert epcsaft.native_regression_contract_schema()["statuses"] == list(CANONICAL_NATIVE_REGRESSION_STATUSES)


def test_runtime_capabilities_reference_native_status_target() -> None:
    contract = epcsaft.capabilities()["regression"]["reactive_electrolyte_batch_context"]["fit_status_contract"]

    assert tuple(contract["canonical_statuses_target"]) == CANONICAL_NATIVE_REGRESSION_STATUSES
    assert "bounded_incomplete" not in contract["canonical_statuses_target"]


def test_native_residual_record_evaluation_preserves_fixed_shape_and_penalties() -> None:
    result = evaluate_native_regression_residual_records(
        [
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
                "predicted": 1.0,
                "observed": 2.0,
                "success": False,
                "recoverable_failure": True,
                "failure_message": "bubble solve failed",
            },
        ],
        penalty_residual=123.0,
    )

    assert result["residuals"] == pytest.approx([0.5, 123.0])
    assert result["success_count"] == 1
    assert result["failure_count"] == 1
    assert result["fixed_shape_residuals"] is True
    assert [entry["residual_index"] for entry in result["residual_schema"]] == [0, 1]
    assert [entry["family"] for entry in result["residual_schema"]] == ["speciation", "pressure"]
    row_2 = {row["row_id"]: row for row in result["row_diagnostics"]}["row_2"]
    assert row_2["success"] is False
    assert row_2["status"] == "line_search_failed"
    assert row_2["penalty_applied"] is True
    assert row_2["message"] == "bubble solve failed"


def test_native_residual_record_evaluation_rejects_invalid_scale() -> None:
    try:
        evaluate_native_regression_residual_records(
            [
                {
                    "row_id": "row_1",
                    "family": "speciation",
                    "target": "CO2",
                    "predicted": 0.45,
                    "observed": 0.40,
                    "scale": 0.0,
                }
            ]
        )
    except ValueError as exc:
        assert "scale" in str(exc)
    else:  # pragma: no cover
        raise AssertionError("invalid scale was accepted")
