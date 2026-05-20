from __future__ import annotations

import numpy as np
import pytest

from epcsaft import ePCSAFTMixture
from epcsaft._types import InputError
from tests.helpers.native_cases import _neutral_state


def _state():
    mix, _species, _pressure, density, temperature, composition = _neutral_state()
    return mix.state(T=temperature, x=composition, rho=density)


def _neutral_binary_state():
    species = ["A", "B"]
    params = {
        "m": np.asarray([1.0000, 1.6069]),
        "s": np.asarray([3.7039, 3.5206]),
        "e": np.asarray([150.03, 191.42]),
        "k_ij": np.asarray([[0.0, 3.0e-4], [3.0e-4, 0.0]]),
        "l_ij": np.asarray([[0.0, 2.0e-4], [2.0e-4, 0.0]]),
    }
    mix = ePCSAFTMixture.from_params(params, species=species)
    return mix.state(T=233.15, x=np.asarray([0.35, 0.65]), rho=1000.0)


def _pure_neutral_state():
    params = {
        "m": np.asarray([1.6069]),
        "s": np.asarray([3.5206]),
        "e": np.asarray([191.42]),
    }
    mix = ePCSAFTMixture.from_params(params, species=["A"])
    return mix.state(T=233.15, x=np.asarray([1.0]), rho=1000.0)


def _centered_pressure_slope(state, step):
    mix = state._mixture
    plus = mix.state(T=state._T, x=state._x, rho=state.molar_density() + step, phase=state._phase).pressure()
    minus = mix.state(T=state._T, x=state._x, rho=state.molar_density() - step, phase=state._phase).pressure()
    return (plus - minus) / (2.0 * step)


def test_pressure_density_derivative_result_reports_backend_contract() -> None:
    state = _state()

    result = state.pressure_density_derivative_result()

    assert set(("supported", "backend", "message", "value", "jacobian", "shape")).issubset(result)
    assert result["backend"] == "cppad"
    assert result["shape"] == [1, 1]
    if result["supported"]:
        assert np.asarray(result["jacobian"], dtype=float).shape == (1, 1)


def test_pressure_density_derivative_result_includes_compressibility_density_dependence() -> None:
    state = _neutral_binary_state()

    result = state.pressure_density_derivative_result()

    reported = float(np.asarray(result["jacobian"], dtype=float)[0, 0])
    step = state.molar_density() * 1.0e-4
    expected = _centered_pressure_slope(state, step)
    ideal_slope = 8.31446261815324 * state._T
    assert not np.isclose(expected, ideal_slope, rtol=1.0e-3, atol=1.0e-6)
    np.testing.assert_allclose(reported, expected, rtol=5.0e-4, atol=1.0e-2)


def test_pressure_unsupported_derivative_routes_raise() -> None:
    state = _state()

    with pytest.raises(InputError, match="unsupported"):
        state.pressure_composition_derivative_result()

    with pytest.raises(InputError, match="unsupported"):
        state.pressure_parameter_derivative_result()


def test_pressure_parameter_derivative_result_supports_neutral_binary_kij() -> None:
    state = _neutral_binary_state()

    result = state.pressure_parameter_derivative_result()

    assert set(("supported", "backend", "message", "value", "jacobian", "shape")).issubset(result)
    assert result["supported"] is True
    assert result["backend"] == "cppad"
    assert result["shape"] == [1, 2]
    assert result["parameter_order"] == ("k_ij:A:B", "l_ij:A:B")
    assert np.asarray(result["jacobian"], dtype=float).shape == (1, 2)


def test_pressure_parameter_derivative_result_supports_pure_neutral_m_sigma_epsilon() -> None:
    state = _pure_neutral_state()

    result = state.pressure_parameter_derivative_result()

    assert result["supported"] is True
    assert result["backend"] == "cppad"
    assert result["shape"] == [1, 3]
    assert result["parameter_order"] == ("m:A", "sigma:A", "epsilon:A")
    assert np.asarray(result["jacobian"], dtype=float).shape == (1, 3)
