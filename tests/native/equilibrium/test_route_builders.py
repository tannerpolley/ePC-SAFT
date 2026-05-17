from __future__ import annotations

import numpy as np
import pytest

import epcsaft
from epcsaft import _core
from tests.helpers.runtime_cases import _ionic_params


def _neutral_binary_mixture() -> epcsaft.ePCSAFTMixture:
    params = {
        "m": np.asarray([1.0, 1.6069]),
        "s": np.asarray([3.7039, 3.5206]),
        "e": np.asarray([150.03, 191.42]),
        "k_ij": np.asarray([[0.0, 3.0e-4], [3.0e-4, 0.0]]),
    }
    return epcsaft.ePCSAFTMixture.from_params(params, species=["Methane", "Ethane"])


def _ionic_mixture() -> epcsaft.ePCSAFTMixture:
    params = _ionic_params()
    params["assoc_scheme"] = [None, None, None]
    params["e_assoc"] = np.zeros(3)
    params["vol_a"] = np.zeros(3)
    return epcsaft.ePCSAFTMixture.from_params(params, species=["water", "Na+", "Cl-"])


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
    assert payload["constraint_lower_bounds"][-1] == pytest.approx(1.0e-8)
    assert payload["constraint_upper_bounds"][-1] > 1.0e6
    assert payload["constraints_at_initial"][-1] >= payload["constraint_lower_bounds"][-1]


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
    assert payload["constraint_lower_bounds"][-1] == pytest.approx(1.0e-8)
    assert payload["constraint_upper_bounds"][-1] > 1.0e6
    assert payload["constraints_at_initial"][-1] >= payload["constraint_lower_bounds"][-1]


def test_electrolyte_lle_route_contract_adds_phase_charge_rows() -> None:
    mix = _ionic_mixture()
    temperature = 298.15
    target_pressure = 1.013e5
    feed_amounts = np.asarray([0.9998, 1.0e-4, 1.0e-4], dtype=float)
    charges = np.asarray([0.0, 1.0, -1.0], dtype=float)

    payload = _core._native_electrolyte_lle_eos_nlp_contract(
        mix._native,
        temperature,
        target_pressure,
        feed_amounts.tolist(),
    )

    initial = np.asarray(payload["initial_point"], dtype=float).reshape(2, 4)
    phase_amounts = initial[:, :3]
    volumes = initial[:, 3]
    phase_system = _core._native_eos_phase_system(
        mix._native,
        temperature,
        target_pressure,
        phase_amounts.tolist(),
        volumes.tolist(),
        feed_amounts.tolist(),
        charges.tolist(),
    )

    assert payload["problem_name"] == "electrolyte_lle_eos"
    assert payload["derivative_backend"] == "analytic_cppad"
    assert payload["phase_count"] == 2
    assert payload["species_count"] == 3
    assert payload["variable_count"] == 8
    assert payload["constraint_count"] == 8
    assert payload["jacobian_nonzero_count"] == 64
    assert np.all(phase_amounts > 0.0)
    assert np.all(volumes > 0.0)
    assert np.sum(phase_amounts, axis=0) == pytest.approx(feed_amounts)
    assert phase_amounts @ charges == pytest.approx([0.0, 0.0], abs=1.0e-14)
    assert payload["constraints_at_initial"][:3] == pytest.approx([0.0, 0.0, 0.0], abs=1.0e-14)
    assert payload["constraints_at_initial"][5:7] == pytest.approx([0.0, 0.0], abs=1.0e-14)
    assert payload["constraint_lower_bounds"][-1] == pytest.approx(1.0e-8)
    assert payload["constraint_upper_bounds"][-1] > 1.0e6
    assert payload["constraints_at_initial"][-1] >= payload["constraint_lower_bounds"][-1]
    assert phase_system["constraint_names"][-2:] == ["phase_0.charge_balance", "phase_1.charge_balance"]
    assert payload["constraints_at_initial"][:-1] == pytest.approx(
        phase_system["constraints"],
        rel=1.0e-12,
        abs=1.0e-8,
    )
    payload_jacobian = np.asarray(payload["jacobian_values_at_initial"], dtype=float).reshape(
        payload["constraint_count"],
        payload["variable_count"],
    )
    assert payload_jacobian[:-1].reshape(-1).tolist() == pytest.approx(
        phase_system["constraint_jacobian_row_major"],
        rel=1.0e-12,
        abs=1.0e-8,
    )
    assert np.count_nonzero(np.abs(payload_jacobian[-1]) > 0.0) > 0


def test_electrolyte_lle_route_result_uses_ipopt_adapter_gate_and_charge_rows() -> None:
    mix = _ionic_mixture()
    payload = _core._native_electrolyte_lle_eos_route_result(
        mix._native,
        298.15,
        1.013e5,
        [0.9998, 1.0e-4, 1.0e-4],
        30,
        1.0e-8,
        0.0,
        1.0e-7,
        1.0e-5,
        1.0e-8,
        1.0e-7,
        1.0e-4,
    )

    assert payload["backend"] == "ipopt"
    assert payload["problem_name"] == "electrolyte_lle_eos"
    assert payload["derivative_backend"] == "analytic_cppad"
    assert payload["exact_gradient_required"] is True
    assert payload["exact_jacobian_required"] is True
    if not payload["compiled"]:
        assert payload["ran"] is False
        assert payload["solver_accepted"] is False
        assert payload["accepted"] is False
        assert payload["status"] == "ipopt_dependency_required"
        assert payload["phase_amounts"] == []
        assert payload["phase_volumes"] == []
        assert payload["postsolve"]["accepted"] is False
        return

    assert payload["ran"] is True
    if not payload["solver_accepted"]:
        assert payload["accepted"] is False
        assert payload["status"] == "solver_rejected"
        return

    assert np.asarray(payload["variables"], dtype=float).shape == (8,)
    assert np.asarray(payload["constraints"], dtype=float).shape == (7,)
    assert np.asarray(payload["phase_amounts"], dtype=float).shape == (2, 3)
    assert np.asarray(payload["phase_volumes"], dtype=float).shape == (2,)
    assert payload["postsolve"]["derivative_backend"] == "analytic_cppad"
    assert payload["postsolve"]["charge_balance_norm"] <= 1.0e-8 or payload["status"] == "postsolve_rejected"


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
    liquid_volume_col = local_variable_count - 1
    vapor_volume_col = 2 * local_variable_count - 1
    lower_bounds = np.asarray(payload["variable_lower_bounds"], dtype=float)
    upper_bounds = np.asarray(payload["variable_upper_bounds"], dtype=float)

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
    assert lower_bounds[liquid_volume_col] <= initial[liquid_volume_col] <= upper_bounds[liquid_volume_col]
    assert lower_bounds[vapor_volume_col] <= initial[vapor_volume_col] <= upper_bounds[vapor_volume_col]
    assert upper_bounds[liquid_volume_col] < lower_bounds[vapor_volume_col]


@pytest.mark.parametrize(
    ("binding_name", "problem_name"),
    [
        ("_native_neutral_bubble_t_eos_nlp_contract", "neutral_bubble_t_eos"),
        ("_native_neutral_dew_t_eos_nlp_contract", "neutral_dew_t_eos"),
    ],
)
def test_neutral_fixed_pressure_temperature_route_contract_pins_specified_phase(
    binding_name: str,
    problem_name: str,
) -> None:
    mix = _neutral_binary_mixture()
    pressure = 1.0e5
    composition = np.asarray([0.35, 0.65], dtype=float)

    payload = getattr(_core, binding_name)(
        mix._native,
        pressure,
        composition.tolist(),
    )

    initial = np.asarray(payload["initial_point"], dtype=float)
    jacobian = np.asarray(payload["jacobian_values_at_initial"], dtype=float).reshape(
        payload["constraint_count"],
        payload["variable_count"],
    )
    local_variable_count = composition.size + 1
    temperature_col = payload["variable_count"] - 1
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
    assert initial[temperature_col] == pytest.approx(300.0)
    assert fixed_amounts / fixed_amounts.sum() == pytest.approx(composition)
    assert payload["constraints_at_initial"][: composition.size + 1] == pytest.approx([0.0, 0.0, 0.0])
    pressure_row_start = composition.size + 1
    assert np.isfinite(jacobian[pressure_row_start, temperature_col])
    assert np.isfinite(jacobian[pressure_row_start + 1, temperature_col])
    assert payload["constraint_lower_bounds"][-1] > 0.0
    assert payload["constraint_upper_bounds"][-1] > payload["constraint_lower_bounds"][-1]
    assert payload["constraints_at_initial"][-1] >= payload["constraint_lower_bounds"][-1]
    assert jacobian[-1, local_variable_count - 1] == pytest.approx(-1.0)
    assert jacobian[-1, 2 * local_variable_count - 1] == pytest.approx(1.0)


def test_electrolyte_bubble_pressure_contract_adds_phase_charge_rows() -> None:
    mix = _ionic_mixture()
    temperature = 298.15
    liquid_composition = np.asarray([0.9998, 1.0e-4, 1.0e-4], dtype=float)
    charges = np.asarray([0.0, 1.0, -1.0], dtype=float)

    payload = _core._native_electrolyte_bubble_p_eos_nlp_contract(
        mix._native,
        temperature,
        liquid_composition.tolist(),
    )

    initial = np.asarray(payload["initial_point"], dtype=float)
    jacobian = np.asarray(payload["jacobian_values_at_initial"], dtype=float).reshape(
        payload["constraint_count"],
        payload["variable_count"],
    )
    local_variable_count = liquid_composition.size + 1
    liquid_amounts = initial[: liquid_composition.size]
    vapor_amounts = initial[local_variable_count : local_variable_count + liquid_composition.size]
    charge_row_start = liquid_composition.size - 1 + 2

    assert payload["problem_name"] == "electrolyte_bubble_p_eos"
    assert payload["derivative_backend"] == "analytic_cppad"
    assert payload["phase_count"] == 2
    assert payload["species_count"] == 3
    assert payload["variable_count"] == 9
    assert payload["constraint_count"] == 12
    assert payload["jacobian_nonzero_count"] == 108
    assert np.all(initial > 0.0)
    assert liquid_amounts / liquid_amounts.sum() == pytest.approx(liquid_composition)
    assert liquid_amounts @ charges == pytest.approx(0.0, abs=1.0e-14)
    assert vapor_amounts @ charges == pytest.approx(0.0, abs=1.0e-14)
    assert payload["constraints_at_initial"][: charge_row_start + 2] == pytest.approx(
        [0.0] * (charge_row_start + 2),
        abs=1.0e-14,
    )
    assert jacobian[charge_row_start, : liquid_composition.size] == pytest.approx(charges)
    assert jacobian[charge_row_start + 1, local_variable_count : local_variable_count + liquid_composition.size] == (
        pytest.approx(charges)
    )


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
        0.0,
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
        assert payload["status"] == "ipopt_dependency_required"
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


@pytest.mark.parametrize(
    ("binding_name", "problem_name"),
    [
        ("_native_neutral_bubble_t_eos_route_result", "neutral_bubble_t_eos"),
        ("_native_neutral_dew_t_eos_route_result", "neutral_dew_t_eos"),
    ],
)
def test_neutral_fixed_pressure_temperature_route_result_uses_ipopt_adapter_gate(
    binding_name: str,
    problem_name: str,
) -> None:
    mix = _neutral_binary_mixture()
    payload = getattr(_core, binding_name)(
        mix._native,
        1.0e5,
        [0.35, 0.65],
        30,
        1.0e-8,
        0.0,
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
        assert payload["status"] == "ipopt_dependency_required"
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
        assert payload["status"] == "ipopt_dependency_required"
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
        0.0,
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
        assert payload["status"] == "ipopt_dependency_required"
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
    assert payload["ln_fugacity_consistency_norm"] == pytest.approx(0.0, abs=1.0e-14)
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
    assert "ln_fugacity_consistency_norm" in payload
    assert np.isfinite(payload["ln_fugacity_consistency_norm"])
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
    assert payload["ln_fugacity_consistency_norm"] > 1.0e-9
    assert payload["phase_distance"] > 1.0e-3
