from __future__ import annotations

import math

import epcsaft


def test_native_thermo_regression_evaluates_reactive_bubble_pressure_row_in_cpp() -> None:
    species = ["H2O", "Na+", "Cl-"]
    mix = epcsaft.ePCSAFTMixture.from_dataset("2026_Khudaida", species, [0.98, 0.01, 0.01], 298.15)

    result = epcsaft.evaluate_native_thermo_regression_rows(
        mix,
        {
            "species": species,
            "rows": [
                {
                    "row_id": "bubble_1",
                    "row_mode": "reactive_electrolyte_bubble",
                    "T": 298.15,
                    "x_liq": [0.98, 0.01, 0.01],
                    "vapor_species": ["H2O"],
                    "options": {"initial_pressure": 101325.0, "return_best_effort": True},
                    "targets": [{"family": "pressure", "target": "P", "observed": 101325.0, "scale": 1.0e-5}],
                }
            ],
        },
    )

    assert result["success_count"] == 1
    assert result["failure_count"] == 0
    assert result["fixed_shape_residuals"] is True
    assert len(result["residuals"]) == 1
    assert math.isfinite(result["residuals"][0])
    assert result["row_diagnostics"][0]["solve_backend"] == "native_electrolyte_bubble"
    assert result["row_diagnostics"][0]["derivative_backend"] == "not_differentiated"
