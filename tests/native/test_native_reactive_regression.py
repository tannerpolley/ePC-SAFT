from __future__ import annotations

from epcsaft.native_regression import native_regression_contract_schema, solve_native_regression_residual_records


def test_native_reactive_regression_contract_covers_pressure_speciation_and_born_kij() -> None:
    schema = native_regression_contract_schema()

    assert {"pressure", "speciation"}.issubset(schema["supported_target_families"])
    assert "born_radius" in schema["supported_parameter_kinds"]
    assert "binary_interaction" in schema["supported_parameter_kinds"]
    assert schema["production_finite_difference_allowed"] is False


def test_native_reactive_regression_minimum_pressure_speciation_slice_runs() -> None:
    result = solve_native_regression_residual_records(
        [
            {
                "row_id": "speciation_row",
                "row_kind": "reactive_speciation",
                "family": "speciation",
                "target": "CO2(aq)",
                "predicted": 0.12,
                "observed": 0.10,
                "scale": 20.0,
            },
            {
                "row_id": "bubble_row",
                "row_kind": "reactive_electrolyte_bubble",
                "family": "pressure",
                "target": "bubble_pressure",
                "predicted": 102000.0,
                "observed": 101325.0,
                "scale": 1.0e-5,
            },
        ],
        [
            {
                "name": "MEA.CO2.k_ij",
                "path": "k_ij[MEA,CO2]",
                "kind": "binary_interaction",
                "initial": 0.02,
                "lower": -0.2,
                "upper": 0.2,
            },
            {
                "name": "MEAH+.d_born",
                "path": "d_born[MEAH+]",
                "kind": "born_radius",
                "initial": 3.2,
                "lower": 2.0,
                "upper": 6.0,
            },
        ],
        options={"optimizer_backend": "auto", "derivative_backend": "analytic"},
    )

    assert result["success"] is True
    assert result["status"] == "converged"
    assert result["objective_result"]["success_count"] == 2
    assert [entry["family"] for entry in result["objective_result"]["residual_schema"]] == ["speciation", "pressure"]
    assert result["parameter_names"] == ["MEA.CO2.k_ij", "MEAH+.d_born"]
