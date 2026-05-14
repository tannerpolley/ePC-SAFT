from __future__ import annotations

import epcsaft


def test_ceres_liquid_electrolyte_regression_reports_backend_unavailable() -> None:
    result = epcsaft.fit_liquid_electrolyte_parameters(
        species=("H2O", "Na+", "Cl-"),
        data_rows=[
            {
                "T": 298.15,
                "P": 101325.0,
                "rho": 55000.0,
                "epsilon_r_exp": 78.3,
                "mean_ionic_activity": 0.75,
            }
        ],
        parameter_set="2026_Khudaida",
        parameters_to_fit=("d_born", "f_solv"),
        initial_guess={"d_born": 3.4, "f_solv": 1.0},
        bounds={"d_born": (1.0, 8.0), "f_solv": (0.1, 3.0)},
        solver_options={"optimizer_backend": "ceres"},
    )

    assert result.success is False
    assert result.backend == "backend_unavailable"
    assert result.optimizer_backend == "ceres"
    assert result.derivative_backend == "backend_unavailable"
    assert result.jacobian_backend == "backend_unavailable"
    assert result.jacobian_fallback_used is False
    assert result.python_objective_used is False
    assert "backend_unavailable" in result.backend_unavailable_reason
    assert result.problem.mode == "liquid_electrolyte"
    assert result.problem.fit_targets == ("d_born", "f_solv")
    assert result.row_diagnostics
    assert all(row["supported"] is False for row in result.row_diagnostics)
