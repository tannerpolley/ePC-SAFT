from __future__ import annotations

import numpy as np
import pytest

import epcsaft
from epcsaft import _core


def _neutral_binary_mixture() -> epcsaft.ePCSAFTMixture:
    params = {
        "m": np.asarray([1.0, 1.6069]),
        "s": np.asarray([3.7039, 3.5206]),
        "e": np.asarray([150.03, 191.42]),
        "k_ij": np.asarray([[0.0, 3.0e-4], [3.0e-4, 0.0]]),
    }
    return epcsaft.ePCSAFTMixture.from_params(params, species=["Methane", "Ethane"])


def test_eos_phase_block_builds_amount_volume_variables_and_pressure_gate() -> None:
    mix = _neutral_binary_mixture()
    amounts = np.asarray([0.8, 1.2], dtype=float)
    temperature = 300.0
    density = 120.0
    volume = float(amounts.sum() / density)
    composition = amounts / amounts.sum()
    state = mix.state(T=temperature, rho=density, x=composition, phase="vapor")
    pressure = state.pressure()

    payload = _core._native_eos_phase_block(mix._native, temperature, pressure, amounts.tolist(), volume)

    assert payload["block"] == "eos_phase"
    assert payload["variable_names"] == ["n_0", "n_1", "volume"]
    assert payload["constraint_names"] == ["pressure_consistency"]
    assert payload["derivative_backend"] == "analytic"
    assert payload["composition"] == pytest.approx(composition.tolist(), abs=1.0e-14)
    assert payload["density"] == pytest.approx(density, rel=1.0e-14)
    assert payload["eos_pressure"] == pytest.approx(pressure, rel=1.0e-14)
    assert abs(payload["pressure_consistency_residual"]) <= 1.0e-8
    assert abs(payload["gradient"][-1]) <= 1.0e-10
    assert payload["objective_terms"]["pressure_work"] == pytest.approx(
        pressure * volume / payload["gas_constant_temperature"]
    )


def test_eos_phase_block_gradient_matches_chemical_potential_and_pressure_identities() -> None:
    mix = _neutral_binary_mixture()
    amounts = np.asarray([0.8, 1.2], dtype=float)
    temperature = 300.0
    density = 120.0
    volume = float(amounts.sum() / density)
    composition = amounts / amounts.sum()
    state = mix.state(T=temperature, rho=density, x=composition, phase="vapor")
    eos_pressure = state.pressure()
    target_pressure = eos_pressure + 2500.0

    payload = _core._native_eos_phase_block(mix._native, temperature, target_pressure, amounts.tolist(), volume)

    residual_mu = np.asarray(state.residual_chemical_potential(), dtype=float)
    expected_amount_gradient = np.log(amounts / volume) + residual_mu
    expected_volume_gradient = (target_pressure - eos_pressure) / payload["gas_constant_temperature"]

    assert payload["gradient"][:-1] == pytest.approx(expected_amount_gradient.tolist(), rel=1.0e-12, abs=1.0e-12)
    assert payload["gradient"][-1] == pytest.approx(expected_volume_gradient, rel=1.0e-12, abs=1.0e-12)
    assert payload["pressure_consistency_residual"] == pytest.approx(eos_pressure - target_pressure, abs=1.0e-8)
    assert np.isfinite(payload["objective"])


def test_eos_phase_block_pressure_jacobian_uses_exact_cppad_homogeneity_identity() -> None:
    mix = _neutral_binary_mixture()
    amounts = np.asarray([0.8, 1.2], dtype=float)
    temperature = 300.0
    density = 120.0
    volume = float(amounts.sum() / density)
    composition = amounts / amounts.sum()
    state = mix.state(T=temperature, rho=density, x=composition, phase="vapor")
    target_pressure = state.pressure() + 2500.0

    payload = _core._native_eos_phase_block(mix._native, temperature, target_pressure, amounts.tolist(), volume)

    jacobian = np.asarray(payload["pressure_jacobian"], dtype=float)
    assert payload["pressure_jacobian_backend"] == "cppad"
    assert payload["pressure_jacobian_shape"] == (1, 3)
    assert jacobian.size == 3
    assert jacobian[-1] == pytest.approx(
        -payload["pressure_density_derivative"] * payload["density"] / payload["volume"],
        rel=1.0e-12,
        abs=1.0e-8,
    )
    assert float(np.dot(amounts, jacobian[:-1]) + volume * jacobian[-1]) == pytest.approx(0.0, abs=1.0e-8)


def test_eos_phase_block_reports_pressure_constraint_jacobian_from_exact_curvature_identity() -> None:
    mix = _neutral_binary_mixture()
    amounts = np.asarray([0.8, 1.2], dtype=float)
    temperature = 300.0
    density = 120.0
    volume = float(amounts.sum() / density)
    composition = amounts / amounts.sum()
    state = mix.state(T=temperature, rho=density, x=composition, phase="vapor")

    payload = _core._native_eos_phase_block(mix._native, temperature, state.pressure(), amounts.tolist(), volume)

    assert payload["objective_curvature_backend"] == "cppad"
    assert payload["constraint_jacobian_backend"] == "cppad"
    assert payload["objective_curvature_shape"] == (3, 3)
    assert payload["constraint_jacobian_shape"] == (1, 3)
    objective_curvature = np.asarray(payload["objective_curvature_row_major"], dtype=float).reshape((3, 3))
    constraint_jacobian = np.asarray(payload["constraint_jacobian_row_major"], dtype=float).reshape((1, 3))
    expected_pressure_jacobian = -payload["gas_constant_temperature"] * objective_curvature[-1, :]

    assert np.all(np.isfinite(objective_curvature))
    assert np.all(np.isfinite(constraint_jacobian))
    assert constraint_jacobian[0, :] == pytest.approx(expected_pressure_jacobian, rel=1.0e-11, abs=1.0e-8)
