from __future__ import annotations

import pytest

import epcsaft
from tests.helpers.regression_cases import _methane_like_records, _minimal_neutral_metadata


def test_ceres_pure_neutral_regression_owns_optimizer_loop() -> None:
    ceres = epcsaft.runtime_build_info()["optional_dependencies"]["ceres"]
    if not ceres["compiled"]:
        pytest.skip("Ceres support is not enabled in this native build.")
    initial_guess = {"m": 1.08, "s": 3.55, "e": 155.0}

    result = epcsaft.fit_pure_neutral(
        _methane_like_records(),
        "Methane",
        assoc_scheme="",
        fixed_parameters=_minimal_neutral_metadata(16.043e-3),
        initial_guess=initial_guess,
        bounds={
            "m": (0.5, 3.5),
            "s": (2.0, 5.0),
            "e": (50.0, 400.0),
        },
        optimizer_backend="ceres",
    )

    assert result.success, result.message
    assert result.optimizer_backend == "ceres"
    assert result.backend == "ceres"
    assert result.derivative_backend != "finite_difference"
    assert result.jacobian_backend != "finite_difference"
    assert result.python_objective_used is False
    assert result.objective_final < result.objective_initial
    assert any(abs(result.fitted_values[name] - initial_guess[name]) > 1.0e-8 for name in ("m", "s", "e"))
    assert result.parameter_map == pytest.approx(result.fitted_values)
    assert result.row_diagnostics
    assert result.n_residual_evaluations >= 1
    assert result.n_jacobian_evaluations >= 1
