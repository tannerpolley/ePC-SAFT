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


def test_pressure_density_derivative_result_reports_backend_contract() -> None:
    state = _state()

    result = state.pressure_density_derivative_result()

    assert set(("supported", "backend", "message", "value", "jacobian", "shape")).issubset(result)
    assert result["backend"] in {"cppad", "not_available"}
    assert result["shape"] == [1, 1]
    assert "numerical_derivative" not in str(result).lower()
    if result["supported"]:
        assert np.asarray(result["jacobian"], dtype=float).shape == (1, 1)


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
    assert "numerical_derivative" not in str(result).lower()
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
    assert "numerical_derivative" not in str(result).lower()
