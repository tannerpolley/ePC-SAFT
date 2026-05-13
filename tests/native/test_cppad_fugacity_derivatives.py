from __future__ import annotations

import numpy as np

from epcsaft import ePCSAFTMixture
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


def test_ln_fugacity_parameter_derivative_result_supports_neutral_binary_kij() -> None:
    state = _neutral_binary_state()

    result = state.ln_fugacity_parameter_derivative_result()

    assert result["supported"] is True
    assert result["backend"] == "cppad"
    assert "finite_difference" not in str(result).lower()
    assert result["shape"] == [state.x.size, 2]
    assert result["parameter_order"] == ("k_ij:A:B", "l_ij:A:B")
    assert result["component_order"] == ("A", "B")
    assert np.asarray(result["jacobian"], dtype=float).shape == (state.x.size, 2)


def test_ln_fugacity_parameter_derivative_result_supports_pure_neutral_m_sigma_epsilon() -> None:
    state = _pure_neutral_state()

    result = state.ln_fugacity_parameter_derivative_result()

    assert result["supported"] is True
    assert result["backend"] == "cppad"
    assert result["shape"] == [1, 3]
    assert result["parameter_order"] == ("m:A", "sigma:A", "epsilon:A")
    assert result["component_order"] == ("A",)
    assert np.asarray(result["jacobian"], dtype=float).shape == (1, 3)
    assert "finite_difference" not in str(result).lower()


def test_chemical_potential_parameter_derivative_result_supports_neutral_binary_kij() -> None:
    state = _neutral_binary_state()

    result = state.chemical_potential_parameter_derivative_result()

    assert result["supported"] is True
    assert result["backend"] == "cppad"
    assert "finite_difference" not in str(result).lower()
    assert result["shape"] == [state.x.size, 2]
    assert result["parameter_order"] == ("k_ij:A:B", "l_ij:A:B")
    assert result["component_order"] == ("A", "B")
    assert result["value_basis"] == "residual_chemical_potential"
    assert np.asarray(result["jacobian"], dtype=float).shape == (state.x.size, 2)


def test_chemical_potential_parameter_derivative_result_supports_pure_neutral_m_sigma_epsilon() -> None:
    state = _pure_neutral_state()

    result = state.chemical_potential_parameter_derivative_result()

    assert result["supported"] is True
    assert result["backend"] == "cppad"
    assert result["shape"] == [1, 3]
    assert result["parameter_order"] == ("m:A", "sigma:A", "epsilon:A")
    assert result["component_order"] == ("A",)
    assert result["value_basis"] == "residual_chemical_potential"
    assert np.asarray(result["jacobian"], dtype=float).shape == (1, 3)
    assert "finite_difference" not in str(result).lower()
