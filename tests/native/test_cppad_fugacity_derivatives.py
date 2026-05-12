from __future__ import annotations

import numpy as np

from tests.helpers.native_cases import _neutral_state


def _state():
    mix, _species, _pressure, density, temperature, composition = _neutral_state()
    return mix.state(T=temperature, x=composition, rho=density)


def test_ares_composition_derivative_result_uses_explicit_backend_labels() -> None:
    state = _state()

    result = state.ares_composition_derivative_result()

    assert result["supported"] is True
    assert result["backend"] in {"analytic", "cppad", "legacy_eigen_forward"}
    assert result["shape"] == [1, state.x.size]
    assert np.asarray(result["jacobian"], dtype=float).shape == (1, state.x.size)
    assert "autodiff" not in str(result["backend"]).lower()
    assert "finite_difference" not in str(result).lower()


def test_ln_fugacity_composition_derivative_result_truthfully_unavailable() -> None:
    state = _state()

    result = state.ln_fugacity_composition_derivative_result()

    assert result["supported"] is False
    assert result["backend"] == "backend_unavailable"
    assert result["shape"] == [state.x.size, state.x.size]
    assert "second composition derivatives" in result["message"]
    assert "finite_difference" not in str(result).lower()


def test_ln_fugacity_parameter_derivative_result_truthfully_unavailable_without_born_path() -> None:
    state = _state()

    result = state.ln_fugacity_parameter_derivative_result()

    assert result["supported"] is False
    assert result["backend"] == "backend_unavailable"
    assert result["shape"] == [state.x.size, 0]
