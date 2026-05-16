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


def test_neutral_tp_flash_route_contract_builds_native_initial_point_from_feed() -> None:
    mix = _neutral_binary_mixture()
    temperature = 300.0
    target_pressure = 1.0e5
    feed_amounts = np.asarray([0.3, 0.7], dtype=float)

    payload = _core._native_neutral_tp_flash_eos_nlp_contract(
        mix._native,
        temperature,
        target_pressure,
        feed_amounts.tolist(),
    )

    initial = np.asarray(payload["initial_point"], dtype=float).reshape(2, 3)
    phase_amounts = initial[:, :2]
    volumes = initial[:, 2]
    phase_compositions = phase_amounts / np.sum(phase_amounts, axis=1, keepdims=True)

    assert payload["problem_name"] == "neutral_two_phase_eos"
    assert payload["derivative_backend"] == "analytic_cppad"
    assert payload["phase_count"] == 2
    assert payload["species_count"] == 2
    assert np.all(phase_amounts > 0.0)
    assert np.all(volumes > 0.0)
    assert np.sum(phase_amounts, axis=0) == pytest.approx(feed_amounts)
    assert phase_compositions[0] != pytest.approx(phase_compositions[1])
    assert np.asarray(payload["constraints_at_initial"], dtype=float)[:2] == pytest.approx([0.0, 0.0])


def test_neutral_lle_route_contract_builds_native_initial_point_from_feed() -> None:
    mix = _neutral_binary_mixture()
    temperature = 298.15
    target_pressure = 1.013e5
    feed_amounts = np.asarray([0.45, 0.55], dtype=float)

    payload = _core._native_neutral_lle_eos_nlp_contract(
        mix._native,
        temperature,
        target_pressure,
        feed_amounts.tolist(),
    )

    initial = np.asarray(payload["initial_point"], dtype=float).reshape(2, 3)
    phase_amounts = initial[:, :2]
    volumes = initial[:, 2]
    phase_compositions = phase_amounts / np.sum(phase_amounts, axis=1, keepdims=True)

    assert payload["problem_name"] == "neutral_two_phase_eos"
    assert payload["derivative_backend"] == "analytic_cppad"
    assert np.all(phase_amounts > 0.0)
    assert np.all(volumes > 0.0)
    assert np.sum(phase_amounts, axis=0) == pytest.approx(feed_amounts)
    assert phase_compositions[0] != pytest.approx(phase_compositions[1])
    assert np.asarray(payload["constraints_at_initial"], dtype=float)[:2] == pytest.approx([0.0, 0.0])


@pytest.mark.parametrize(
    ("binding_name", "problem_name"),
    [
        ("_native_neutral_bubble_p_eos_nlp_contract", "neutral_bubble_p_eos"),
        ("_native_neutral_dew_p_eos_nlp_contract", "neutral_dew_p_eos"),
    ],
)
def test_neutral_fixed_temperature_pressure_route_contract_pins_specified_phase(
    binding_name: str,
    problem_name: str,
) -> None:
    mix = _neutral_binary_mixture()
    temperature = 300.0
    composition = np.asarray([0.35, 0.65], dtype=float)

    payload = getattr(_core, binding_name)(
        mix._native,
        temperature,
        composition.tolist(),
    )

    initial = np.asarray(payload["initial_point"], dtype=float)
    jacobian = np.asarray(payload["jacobian_values_at_initial"], dtype=float).reshape(
        payload["constraint_count"],
        payload["variable_count"],
    )
    local_variable_count = composition.size + 1
    pressure_col = payload["variable_count"] - 1
    first_amounts = initial[: composition.size]
    second_amounts = initial[local_variable_count : local_variable_count + composition.size]
    fixed_amounts = first_amounts if "bubble" in problem_name else second_amounts

    assert payload["problem_name"] == problem_name
    assert payload["derivative_backend"] == "analytic_cppad"
    assert payload["phase_count"] == 2
    assert payload["species_count"] == composition.size
    assert payload["variable_count"] == 2 * local_variable_count + 1
    assert payload["constraint_count"] == 2 * composition.size + 4
    assert payload["jacobian_nonzero_count"] == payload["variable_count"] * payload["constraint_count"]
    assert np.all(initial > 0.0)
    assert fixed_amounts / fixed_amounts.sum() == pytest.approx(composition)
    assert payload["constraints_at_initial"][: composition.size + 1] == pytest.approx([0.0, 0.0, 0.0])
    pressure_row_start = composition.size + 1
    assert jacobian[pressure_row_start, pressure_col] == pytest.approx(-1.0)
    assert jacobian[pressure_row_start + 1, pressure_col] == pytest.approx(-1.0)
    assert payload["constraint_lower_bounds"][-1] > 0.0
    assert payload["constraint_upper_bounds"][-1] > payload["constraint_lower_bounds"][-1]
    assert payload["constraints_at_initial"][-1] >= payload["constraint_lower_bounds"][-1]
    assert jacobian[-1, local_variable_count - 1] == pytest.approx(-1.0)
    assert jacobian[-1, 2 * local_variable_count - 1] == pytest.approx(1.0)


@pytest.mark.parametrize(
    ("binding_name", "problem_name"),
    [
        ("_native_neutral_bubble_p_eos_route_result", "neutral_bubble_p_eos"),
        ("_native_neutral_dew_p_eos_route_result", "neutral_dew_p_eos"),
    ],
)
def test_neutral_fixed_temperature_pressure_route_result_uses_ipopt_adapter_gate(
    binding_name: str,
    problem_name: str,
) -> None:
    mix = _neutral_binary_mixture()
    payload = getattr(_core, binding_name)(
        mix._native,
        300.0,
        [0.35, 0.65],
        30,
        1.0e-8,
        1.0e-7,
        1.0e-5,
        1.0e-7,
        1.0e-4,
    )

    assert payload["backend"] == "ipopt"
    assert payload["problem_name"] == problem_name
    assert payload["derivative_backend"] == "analytic_cppad"
    assert payload["exact_gradient_required"] is True
    assert payload["exact_jacobian_required"] is True
    if not payload["compiled"]:
        assert payload["ran"] is False
        assert payload["solver_accepted"] is False
        assert payload["accepted"] is False
        assert payload["status"] == "requires_ipopt_build"
        assert payload["phase_amounts"] == []
        assert payload["phase_volumes"] == []
        assert "fixed_composition_norm" in payload["postsolve"]
        assert "phase_amount_total_norm" in payload["postsolve"]
        return

    assert payload["ran"] is True
    if not payload["solver_accepted"]:
        assert payload["accepted"] is False
        assert payload["status"] == "solver_rejected"
        return

    assert np.asarray(payload["variables"], dtype=float).shape == (7,)
    assert np.asarray(payload["phase_amounts"], dtype=float).shape == (2, 2)
    assert np.asarray(payload["phase_volumes"], dtype=float).shape == (2,)
    assert payload["postsolve"]["derivative_backend"] == "analytic_cppad"
    assert payload["status"] in {"accepted", "solver_rejected", "postsolve_rejected"}


def test_neutral_two_phase_eos_route_result_translates_solver_and_postsolve() -> None:
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

    payload = _core._native_neutral_two_phase_eos_route_result(
        mix._native,
        temperature,
        target_pressure,
        [phase.tolist() for phase in phase_amounts],
        volumes,
        feed_amounts.tolist(),
        30,
        1.0e-8,
        1.0e-7,
        1.0e-5,
        1.0e-7,
        1.0e-4,
    )

    assert payload["backend"] == "ipopt"
    assert payload["problem_name"] == "neutral_two_phase_eos"
    assert payload["derivative_backend"] == "analytic_cppad"
    assert payload["exact_gradient_required"] is True
    assert payload["exact_jacobian_required"] is True
    if not payload["compiled"]:
        assert payload["ran"] is False
        assert payload["solver_accepted"] is False
        assert payload["accepted"] is False
        assert payload["status"] == "requires_ipopt_build"
        assert payload["phase_amounts"] == []
        assert payload["phase_volumes"] == []
        assert payload["postsolve"]["accepted"] is False
        return

    assert payload["ran"] is True
    assert np.asarray(payload["variables"], dtype=float).shape == (6,)
    assert np.asarray(payload["phase_amounts"], dtype=float).shape == (2, 2)
    assert np.asarray(payload["phase_volumes"], dtype=float).shape == (2,)
    assert payload["postsolve"]["derivative_backend"] == "analytic_cppad"
    assert payload["accepted"] == (payload["solver_accepted"] and payload["postsolve"]["accepted"])
    if payload["accepted"]:
        assert payload["status"] == "accepted"
        assert payload["postsolve"]["rejection_reason"] == "accepted"
    else:
        assert payload["status"] in {"solver_rejected", "postsolve_rejected"}


def test_neutral_lle_route_result_uses_ipopt_adapter_gate() -> None:
    mix = _neutral_binary_mixture()
    payload = _core._native_neutral_lle_eos_route_result(
        mix._native,
        298.15,
        1.013e5,
        [0.45, 0.55],
        30,
        1.0e-8,
        1.0e-7,
        1.0e-5,
        1.0e-7,
        1.0e-4,
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
    assert np.asarray(payload["phase_amounts"], dtype=float).shape == (2, 2)
    assert np.asarray(payload["phase_volumes"], dtype=float).shape == (2,)
    assert payload["status"] in {"accepted", "solver_rejected", "postsolve_rejected"}


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
