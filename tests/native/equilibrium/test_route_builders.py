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


def test_neutral_two_phase_eos_nlp_contract_uses_phase_system_blocks() -> None:
    mix = _neutral_binary_mixture()
    temperature = 300.0
    phase_amounts = [
        np.asarray([0.7, 0.3], dtype=float),
        np.asarray([0.1, 0.9], dtype=float),
    ]
    volumes = [float(phase_amounts[0].sum() / 80.0), float(phase_amounts[1].sum() / 140.0)]
    feed_amounts = phase_amounts[0] + phase_amounts[1]
    target_pressure = mix.state(
        T=temperature,
        rho=phase_amounts[0].sum() / volumes[0],
        x=phase_amounts[0] / phase_amounts[0].sum(),
        phase="liquid",
    ).pressure()

    payload = _core._native_neutral_two_phase_eos_nlp_contract(
        mix._native,
        temperature,
        target_pressure,
        [phase.tolist() for phase in phase_amounts],
        volumes,
        feed_amounts.tolist(),
    )
    phase_system = _core._native_eos_phase_system(
        mix._native,
        temperature,
        target_pressure,
        [phase.tolist() for phase in phase_amounts],
        volumes,
        feed_amounts.tolist(),
    )

    assert payload["problem_name"] == "neutral_two_phase_eos"
    assert payload["derivative_backend"] == "analytic_cppad"
    assert payload["phase_count"] == 2
    assert payload["species_count"] == 2
    assert payload["variable_count"] == 6
    assert payload["constraint_count"] == 4
    assert payload["jacobian_nonzero_count"] == 24
    assert payload["initial_point"] == pytest.approx([0.7, 0.3, volumes[0], 0.1, 0.9, volumes[1]])
    assert payload["objective_at_initial"] == pytest.approx(phase_system["objective"])
    assert payload["gradient_at_initial"] == pytest.approx(phase_system["gradient"], rel=1.0e-12, abs=1.0e-12)
    assert payload["constraints_at_initial"] == pytest.approx(
        phase_system["constraints"],
        rel=1.0e-12,
        abs=1.0e-8,
    )
    assert payload["jacobian_values_at_initial"] == pytest.approx(
        phase_system["constraint_jacobian_row_major"],
        rel=1.0e-12,
        abs=1.0e-8,
    )
    assert len(payload["variable_lower_bounds"]) == payload["variable_count"]
    assert len(payload["variable_upper_bounds"]) == payload["variable_count"]
    assert len(payload["constraint_lower_bounds"]) == payload["constraint_count"]
    assert len(payload["constraint_upper_bounds"]) == payload["constraint_count"]
    assert np.all(np.asarray(payload["variable_lower_bounds"], dtype=float) > 0.0)
    assert payload["constraint_lower_bounds"] == pytest.approx([0.0, 0.0, 0.0, 0.0])
    assert payload["constraint_upper_bounds"] == pytest.approx([0.0, 0.0, 0.0, 0.0])


def test_neutral_two_phase_eos_ipopt_solve_is_native_and_dependency_gated() -> None:
    mix = _neutral_binary_mixture()
    temperature = 300.0
    phase_amounts = [
        np.asarray([0.4, 0.6], dtype=float),
        np.asarray([0.4, 0.6], dtype=float),
    ]
    density = 100.0
    volumes = [float(phase.sum() / density) for phase in phase_amounts]
    feed_amounts = phase_amounts[0] + phase_amounts[1]
    target_pressure = mix.state(
        T=temperature,
        rho=density,
        x=phase_amounts[0] / phase_amounts[0].sum(),
        phase="liquid",
    ).pressure()

    payload = _core._native_neutral_two_phase_eos_ipopt_solve(
        mix._native,
        temperature,
        target_pressure,
        [phase.tolist() for phase in phase_amounts],
        volumes,
        feed_amounts.tolist(),
        30,
        1.0e-8,
    )

    assert payload["backend"] == "ipopt"
    assert payload["problem_name"] == "neutral_two_phase_eos"
    assert payload["derivative_backend"] == "analytic_cppad"
    assert payload["exact_gradient_required"] is True
    assert payload["exact_jacobian_required"] is True
    if not payload["compiled"]:
        assert payload["ran"] is False
        assert payload["accepted"] is False
        assert payload["status"] == "requires_ipopt_build"
        return

    assert payload["ran"] is True
    assert payload["accepted"] is True
    assert payload["material_balance_norm"] <= 1.0e-7
    assert payload["pressure_consistency_norm"] <= 1.0e-3
    assert np.asarray(payload["variables"], dtype=float).shape == (6,)


def test_neutral_two_phase_eos_postsolve_rejects_collapsed_phases() -> None:
    mix = _neutral_binary_mixture()
    temperature = 300.0
    phase_amounts = [
        np.asarray([0.4, 0.6], dtype=float),
        np.asarray([0.4, 0.6], dtype=float),
    ]
    density = 120.0
    volumes = [float(phase.sum() / density) for phase in phase_amounts]
    feed_amounts = phase_amounts[0] + phase_amounts[1]
    target_pressure = mix.state(
        T=temperature,
        rho=density,
        x=phase_amounts[0] / phase_amounts[0].sum(),
        phase="liquid",
    ).pressure()

    payload = _core._native_neutral_two_phase_eos_postsolve(
        mix._native,
        temperature,
        target_pressure,
        [phase.tolist() for phase in phase_amounts],
        volumes,
        feed_amounts.tolist(),
        1.0e-8,
        1.0e-6,
        1.0e-6,
        1.0e-3,
    )

    assert payload["accepted"] is False
    assert payload["rejection_reason"] == "phase_distance"
    assert payload["material_balance_norm"] <= 1.0e-12
    assert payload["pressure_consistency_norm"] <= 1.0e-6
    assert payload["chemical_potential_consistency_norm"] == pytest.approx(0.0, abs=1.0e-14)
    assert payload["phase_distance"] == pytest.approx(0.0, abs=1.0e-14)


def test_neutral_two_phase_eos_postsolve_reports_pressure_gate() -> None:
    mix = _neutral_binary_mixture()
    temperature = 300.0
    phase_amounts = [
        np.asarray([0.7, 0.3], dtype=float),
        np.asarray([0.1, 0.9], dtype=float),
    ]
    volumes = [float(phase_amounts[0].sum() / 80.0), float(phase_amounts[1].sum() / 140.0)]
    feed_amounts = phase_amounts[0] + phase_amounts[1]
    target_pressure = mix.state(
        T=temperature,
        rho=phase_amounts[0].sum() / volumes[0],
        x=phase_amounts[0] / phase_amounts[0].sum(),
        phase="liquid",
    ).pressure()

    payload = _core._native_neutral_two_phase_eos_postsolve(
        mix._native,
        temperature,
        target_pressure,
        [phase.tolist() for phase in phase_amounts],
        volumes,
        feed_amounts.tolist(),
        1.0e-8,
        1.0e-6,
        1.0e-6,
        1.0e-3,
    )

    assert payload["accepted"] is False
    assert payload["rejection_reason"] == "pressure_consistency"
    assert payload["material_balance_norm"] <= 1.0e-12
    assert payload["pressure_consistency_norm"] > 1.0e-6
    assert "chemical_potential_consistency_norm" in payload
    assert payload["phase_distance"] > 1.0e-3


def test_neutral_two_phase_eos_postsolve_reports_chemical_potential_gate() -> None:
    mix = _neutral_binary_mixture()
    temperature = 300.0
    phase_amounts = [
        np.asarray([0.7, 0.3], dtype=float),
        np.asarray([0.1, 0.9], dtype=float),
    ]
    density = 120.0
    volumes = [float(phase.sum() / density) for phase in phase_amounts]
    feed_amounts = phase_amounts[0] + phase_amounts[1]
    target_pressure = mix.state(
        T=temperature,
        rho=density,
        x=phase_amounts[0] / phase_amounts[0].sum(),
        phase="liquid",
    ).pressure()

    payload = _core._native_neutral_two_phase_eos_postsolve(
        mix._native,
        temperature,
        target_pressure,
        [phase.tolist() for phase in phase_amounts],
        volumes,
        feed_amounts.tolist(),
        1.0e-8,
        1.0e12,
        1.0e-9,
        1.0e-3,
    )

    assert payload["accepted"] is False
    assert payload["rejection_reason"] == "chemical_potential_consistency"
    assert payload["material_balance_norm"] <= 1.0e-12
    assert payload["pressure_consistency_norm"] <= 1.0e12
    assert payload["chemical_potential_consistency_norm"] > 1.0e-9
    assert payload["phase_distance"] > 1.0e-3
