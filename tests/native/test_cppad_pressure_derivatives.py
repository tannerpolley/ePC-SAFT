from __future__ import annotations

import numpy as np

from tests.helpers.native_cases import _neutral_state


def _state():
    mix, _species, _pressure, density, temperature, composition = _neutral_state()
    return mix.state(T=temperature, x=composition, rho=density)


def test_pressure_density_derivative_result_reports_backend_contract() -> None:
    state = _state()

    result = state.pressure_density_derivative_result()

    assert set(("supported", "backend", "message", "value", "jacobian", "shape")).issubset(result)
    assert result["backend"] in {"cppad", "backend_unavailable"}
    assert result["shape"] == [1, 1]
    assert "finite_difference" not in str(result).lower()
    if result["supported"]:
        assert np.asarray(result["jacobian"], dtype=float).shape == (1, 1)


def test_pressure_unsupported_derivative_results_are_explicit() -> None:
    state = _state()

    composition = state.pressure_composition_derivative_result()
    parameters = state.pressure_parameter_derivative_result()

    assert composition["supported"] is False
    assert composition["backend"] == "backend_unavailable"
    assert composition["shape"] == [1, state.x.size]
    assert parameters["supported"] is False
    assert parameters["backend"] == "backend_unavailable"
    assert "finite_difference" not in str(composition).lower()
