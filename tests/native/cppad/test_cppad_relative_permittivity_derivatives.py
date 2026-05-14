from __future__ import annotations

import numpy as np

from tests.helpers.native_cases import _ionic_state


def _state():
    mix, _species, _pressure, density, temperature, composition = _ionic_state()
    return mix.state(T=temperature, x=composition, rho=density)


def test_relative_permittivity_composition_derivative_result_is_analytic() -> None:
    state = _state()

    result = state.relative_permittivity_composition_derivative_result()

    assert result["supported"] is True
    assert result["backend"] in {"analytic", "legacy_eigen_forward"}
    assert result["shape"] == [1, state.x.size]
    assert np.asarray(result["jacobian"], dtype=float).shape == (1, state.x.size)
    assert "finite_difference" not in str(result).lower()


def test_relative_permittivity_parameter_derivative_result_for_linear_rule() -> None:
    state = _state()

    result = state.relative_permittivity_parameter_derivative_result()

    assert result["supported"] is True
    assert result["backend"] == "analytic"
    assert result["shape"] == [1, state.x.size]
    np.testing.assert_allclose(np.asarray(result["jacobian"], dtype=float), state.x.reshape(1, -1))
