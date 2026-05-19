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


def _ascani_electrolyte_mixture() -> tuple[epcsaft.ePCSAFTMixture, list[float]]:
    aq = np.asarray([0.798324680201737, 0.016320352824141723, 0.09267748348706063, 0.09267748348706063])
    org = np.asarray([0.37006036048879404, 0.6214918588210971, 0.004223890345054407, 0.004223890345054407])
    beta_org = 0.613766575013417
    feed = ((1.0 - beta_org) * aq + beta_org * org).tolist()
    mix = epcsaft.ePCSAFTMixture.from_dataset("2022_Ascani", ["H2O", "Butanol", "Na+", "Cl-"], feed, 298.15)
    return mix, feed


def _dense_jacobian_from_sparse_contract(payload: dict) -> np.ndarray:
    dense = np.zeros((payload["constraint_count"], payload["variable_count"]), dtype=float)
    rows = np.asarray(payload["jacobian_rows"], dtype=int)
    cols = np.asarray(payload["jacobian_cols"], dtype=int)
    values = np.asarray(payload["jacobian_values_at_initial"], dtype=float)
    assert rows.shape == cols.shape == values.shape == (payload["jacobian_nonzero_count"],)
    for row, col, value in zip(rows, cols, values, strict=True):
        dense[row, col] += value
    return dense


def test_neutral_tp_flash_route_uses_exact_hessian_when_requested() -> None:
    mix = _neutral_binary_mixture()
    payload = _core._native_neutral_tp_flash_eos_route_result(
        mix._native,
        300.0,
        1.0e5,
        [0.4, 0.6],
        30,
        1.0e-8,
        0.0,
        "exact",
        20,
        1.0e-7,
        1.0e-5,
        1.0e-7,
        1.0e-4,
        None,
    )

    if not payload["compiled"]:
        assert payload["status"] == "ipopt_dependency_required"
        return

    assert payload["ran"] is True
    assert payload["solver_accepted"] is True
    assert payload["hessian_approximation"] == "exact"
    assert payload["exact_hessian_available"] is True
    assert payload["hessian_backend"] == "cppad_phase_system"
    assert payload["eval_h_calls"] > 0


def test_neutral_lle_route_uses_exact_hessian_when_requested() -> None:
    mix = _neutral_binary_mixture()
    payload = _core._native_neutral_lle_eos_route_result(
        mix._native,
        298.15,
        1.013e5,
        [0.45, 0.55],
        30,
        1.0e-8,
        0.0,
        "exact",
        20,
        1.0e-7,
        1.0e-5,
        1.0e-7,
        1.0e-4,
        None,
    )

    if not payload["compiled"]:
        assert payload["status"] == "ipopt_dependency_required"
        return

    assert payload["ran"] is True
    assert payload["solver_accepted"] is True
    assert payload["hessian_approximation"] == "exact"
    assert payload["exact_hessian_available"] is True
    assert payload["hessian_backend"] == "cppad_phase_system"
    assert payload["eval_h_calls"] > 0


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
    assert payload["jacobian_nonzero_count"] == 10
    assert payload["initial_point"] == pytest.approx([0.7, 0.3, volumes[0], 0.1, 0.9, volumes[1]])
    assert payload["objective_at_initial"] == pytest.approx(phase_system["objective"])
    assert payload["gradient_at_initial"] == pytest.approx(phase_system["gradient"], rel=1.0e-12, abs=1.0e-12)
    assert payload["constraints_at_initial"] == pytest.approx(
        phase_system["constraints"],
        rel=1.0e-12,
        abs=1.0e-8,
    )
    np.testing.assert_allclose(
        _dense_jacobian_from_sparse_contract(payload),
        np.asarray(phase_system["constraint_jacobian_row_major"], dtype=float).reshape(4, 6),
        rtol=1.0e-12,
        atol=1.0e-8,
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


def test_neutral_stability_tpd_contract_builds_exact_native_nlp() -> None:
    mix = _neutral_binary_mixture()
    feed = np.asarray([0.3, 0.7], dtype=float)

    payload = _core._native_neutral_stability_tpd_nlp_contract(
        mix._native,
        300.0,
        1.0e5,
        feed.tolist(),
        "vap",
        "vap",
    )

    initial = np.asarray(payload["initial_point"], dtype=float)
    gradient = np.asarray(payload["gradient_at_initial"], dtype=float)
    jacobian = np.asarray(payload["jacobian_values_at_initial"], dtype=float)

    assert payload["problem_name"] == "neutral_stability_tpd"
    assert payload["derivative_backend"] == "cppad_implicit"
    assert payload["parent_phase"] == "vap"
    assert payload["trial_phase"] == "vap"
    assert payload["species_count"] == 2
    assert payload["variable_count"] == 2
    assert payload["constraint_count"] == 1
    assert payload["jacobian_nonzero_count"] == 2
    assert payload["feed_composition"] == pytest.approx(feed)
    assert len(payload["parent_reduced_potential"]) == 2
    assert np.all(initial > 0.0)
    assert initial.sum() == pytest.approx(1.0)
    assert initial != pytest.approx(feed)
    assert payload["constraints_at_initial"] == pytest.approx([0.0])
    assert np.all(np.isfinite(gradient))
    assert jacobian == pytest.approx([1.0, 1.0])


def test_neutral_stability_tpd_route_result_uses_ipopt_adapter_gate() -> None:
    mix = _neutral_binary_mixture()
    payload = _core._native_neutral_stability_tpd_route_result(
        mix._native,
        300.0,
        1.0e5,
        [0.3, 0.7],
        "vap",
        "vap",
        30,
        1.0e-8,
        0.0,
        "limited-memory",
        20,
        1.0e-8,
        [],
        None,
        linear_solver="mumps",
        acceptable_tolerance=9.0e-7,
        constraint_violation_tolerance=8.0e-8,
        dual_infeasibility_tolerance=7.0e-8,
        complementarity_tolerance=6.0e-8,
    )

    assert payload["backend"] == "ipopt"
    assert payload["problem_name"] == "neutral_stability_tpd"
    assert payload["derivative_backend"] == "cppad_implicit"
    assert payload["parent_phase"] == "vap"
    assert payload["trial_phase"] == "vap"
    assert payload["exact_gradient_required"] is True
    assert payload["exact_jacobian_required"] is True
    if not payload["compiled"]:
        assert payload["ran"] is False
        assert payload["solver_accepted"] is False
        assert payload["accepted"] is False
        assert payload["status"] == "ipopt_dependency_required"
        assert payload["trial_composition"] == []
        return

    assert payload["ran"] is True
    if not payload["solver_accepted"]:
        assert payload["accepted"] is False
        assert payload["status"] == "solver_rejected"
        return
    assert payload["accepted"] is True
    assert payload["status"] == "accepted"
    assert np.asarray(payload["trial_composition"], dtype=float).sum() == pytest.approx(1.0, abs=1.0e-7)


def test_neutral_stability_tpd_route_uses_exact_hessian_when_requested() -> None:
    mix = _neutral_binary_mixture()
    payload = _core._native_neutral_stability_tpd_route_result(
        mix._native,
        300.0,
        1.0e5,
        [0.3, 0.7],
        "vap",
        "vap",
        30,
        1.0e-8,
        0.0,
        "exact",
        20,
        1.0e-8,
        [],
        None,
    )

    if not payload["compiled"]:
        assert payload["status"] == "ipopt_dependency_required"
        return

    assert payload["ran"] is True
    assert payload["solver_accepted"] is True
    assert payload["accepted"] is True
    assert payload["status"] == "accepted"
    assert payload["hessian_approximation"] == "exact"
    assert payload["exact_hessian_available"] is True
    assert payload["hessian_backend"] != "limited-memory"
    assert payload["eval_h_calls"] > 0


def test_electrolyte_stability_tpd_contract_adds_charge_constraint() -> None:
    mix = _ionic_mixture()
    feed = np.asarray([0.9998, 1.0e-4, 1.0e-4], dtype=float)
    charges = np.asarray([0.0, 1.0, -1.0], dtype=float)

    payload = _core._native_electrolyte_stability_tpd_nlp_contract(
        mix._native,
        298.15,
        1.013e5,
        feed.tolist(),
    )

    initial = np.asarray(payload["initial_point"], dtype=float)
    gradient = np.asarray(payload["gradient_at_initial"], dtype=float)
    jacobian = np.asarray(payload["jacobian_values_at_initial"], dtype=float).reshape(2, 3)

    assert payload["problem_name"] == "electrolyte_stability_tpd"
    assert payload["derivative_backend"] == "cppad_implicit"
    assert payload["parent_phase"] == "liq"
    assert payload["trial_phase"] == "liq"
    assert payload["species_count"] == 3
    assert payload["variable_count"] == 3
    assert payload["constraint_count"] == 2
    assert payload["jacobian_nonzero_count"] == 6
    assert payload["feed_composition"] == pytest.approx(feed)
    assert len(payload["parent_reduced_potential"]) == 3
    assert initial == pytest.approx(feed)
    assert np.dot(initial, charges) == pytest.approx(0.0, abs=1.0e-14)
    assert payload["constraints_at_initial"] == pytest.approx([0.0, 0.0], abs=1.0e-14)
    assert np.all(np.isfinite(gradient))
    assert jacobian[0] == pytest.approx([1.0, 1.0, 1.0])
    assert jacobian[1] == pytest.approx(charges)


def test_electrolyte_stability_tpd_route_result_uses_ipopt_adapter_gate() -> None:
    mix = _ionic_mixture()
    payload = _core._native_electrolyte_stability_tpd_route_result(
        mix._native,
        298.15,
        1.013e5,
        [0.9998, 1.0e-4, 1.0e-4],
        30,
        1.0e-8,
        0.0,
        "limited-memory",
        20,
        1.0e-8,
        [],
        None,
        linear_solver="mumps",
        acceptable_tolerance=9.0e-7,
        constraint_violation_tolerance=8.0e-8,
        dual_infeasibility_tolerance=7.0e-8,
        complementarity_tolerance=6.0e-8,
    )

    assert payload["backend"] == "ipopt"
    assert payload["problem_name"] == "electrolyte_stability_tpd"
    assert payload["derivative_backend"] == "cppad_implicit"
    assert payload["parent_phase"] == "liq"
    assert payload["trial_phase"] == "liq"
    assert payload["seed_name"] == "canonical_charge_neutral_feed"
    assert payload["exact_gradient_required"] is True
    assert payload["exact_jacobian_required"] is True
    if not payload["compiled"]:
        assert payload["ran"] is False
        assert payload["solver_accepted"] is False
        assert payload["accepted"] is False
        assert payload["status"] == "ipopt_dependency_required"
        assert payload["trial_composition"] == []
        return

    assert payload["ran"] is True
    if not payload["solver_accepted"]:
        assert payload["accepted"] is False
        assert payload["status"] == "solver_rejected"
        return
    assert payload["accepted"] is True
    assert payload["status"] == "accepted"
    trial = np.asarray(payload["trial_composition"], dtype=float)
    assert trial.sum() == pytest.approx(1.0, abs=1.0e-7)
    assert np.dot(trial, np.asarray([0.0, 1.0, -1.0])) == pytest.approx(0.0, abs=1.0e-7)


def test_electrolyte_stability_tpd_route_uses_exact_hessian_when_requested() -> None:
    mix = _ionic_mixture()
    payload = _core._native_electrolyte_stability_tpd_route_result(
        mix._native,
        298.15,
        1.013e5,
        [0.9998, 1.0e-4, 1.0e-4],
        30,
        1.0e-8,
        0.0,
        "exact",
        20,
        1.0e-8,
        [],
        None,
    )

    if not payload["compiled"]:
        assert payload["status"] == "ipopt_dependency_required"
        return

    assert payload["ran"] is True
    assert payload["solver_accepted"] is True
    assert payload["accepted"] is True
    assert payload["status"] == "accepted"
    assert payload["hessian_approximation"] == "exact"
    assert payload["exact_hessian_available"] is True
    assert payload["hessian_backend"] != "limited-memory"
    assert payload["eval_h_calls"] > 0


def test_electrolyte_stability_exact_hessian_dilute_salt_route_keeps_callback_finite() -> None:
    formula = np.asarray(
        [0.6732574103166201, 0.0354880934611478, 0.2323336121370503, 0.05892088408518178],
        dtype=float,
    )
    neutrals = formula[:3] / float(np.sum(formula[:3]))
    salt = 1.0e-6
    formula = np.asarray([*(neutrals * (1.0 - salt)), salt], dtype=float)
    feed = np.asarray([formula[0], formula[1], formula[2], formula[3], formula[3]], dtype=float)
    feed = feed / float(np.sum(feed))
    mix = epcsaft.ePCSAFTMixture.from_dataset(
        "2026_Khudaida",
        ["H2O", "Ethanol", "Butanol", "Na+", "Cl-"],
        feed,
        303.15,
    )

    payload = _core._native_electrolyte_stability_tpd_route_result(
        mix._native,
        303.15,
        1.0e5,
        feed.tolist(),
        30,
        1.0e-8,
        0.0,
        "exact",
        10,
        1.0e-8,
        feed.tolist(),
        None,
    )

    if not payload["compiled"]:
        assert payload["status"] == "ipopt_dependency_required"
        return

    assert payload["ran"] is True
    assert payload["hessian_approximation"] == "exact"
    assert payload["exact_hessian_available"] is True
    assert payload["eval_h_calls"] > 0
    assert payload["solver_status"] != "invalid_number_detected"
    assert payload["last_callback_exception"] == ""
    assert payload["last_callback_failure"] == ""


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


def test_electrolyte_lle_route_contract_uses_liquid_root_transformed_variables() -> None:
    mix, feed_amounts = _ascani_electrolyte_mixture()
    temperature = 298.15
    target_pressure = 1.0e5

    payload = _core._native_electrolyte_lle_eos_nlp_contract(
        mix._native,
        temperature,
        target_pressure,
        feed_amounts,
    )

    assert payload["problem_name"] == "electrolyte_lle_eos"
    assert payload["derivative_backend"] == "cppad_implicit"
    assert payload["density_backend"] == "liquid_pressure_root"
    assert payload["phase_count"] == 2
    assert payload["species_count"] == 4
    assert payload["variable_model"] == "ascani_transformed_salt_pairs"
    assert payload["variable_count"] == 3
    assert payload["constraint_count"] == 7
    assert payload["jacobian_nonzero_count"] == 21
    assert len(payload["initial_point"]) == 3
    assert len(payload["variable_lower_bounds"]) == 3
    assert len(payload["variable_upper_bounds"]) == 3
    assert np.allclose(payload["constraint_lower_bounds"][:3], 0.0)
    assert np.allclose(payload["constraint_upper_bounds"][:3], 0.0)
    assert payload["constraint_lower_bounds"][3] >= 0.1
    assert payload["constraint_upper_bounds"][3] > 1.0e6
    assert payload["constraints_at_initial"][3] >= payload["constraint_lower_bounds"][3]
    assert np.all(np.asarray(payload["constraints_at_initial"][4:], dtype=float) > 0.0)
    payload_jacobian = np.asarray(payload["jacobian_values_at_initial"], dtype=float).reshape(
        payload["constraint_count"],
        payload["variable_count"],
    )
    assert np.all(np.isfinite(payload_jacobian))
    assert np.count_nonzero(np.abs(payload_jacobian[0]) > 0.0) > 0
    assert np.count_nonzero(np.abs(payload_jacobian[3]) > 0.0) > 0


def test_electrolyte_lle_route_result_uses_ipopt_adapter_gate_and_charge_rows() -> None:
    mix, feed = _ascani_electrolyte_mixture()
    payload = _core._native_electrolyte_lle_eos_route_result(
        mix._native,
        298.15,
        1.0e5,
        feed,
        500,
        1.0e-8,
        0.0,
        "limited-memory",
        20,
        1.0e-8,
        1.0e-3,
        1.0e-7,
        1.0e-6,
        0.1,
        None,
        linear_solver="mumps",
        acceptable_tolerance=9.0e-7,
        constraint_violation_tolerance=8.0e-8,
        dual_infeasibility_tolerance=7.0e-8,
        complementarity_tolerance=6.0e-8,
    )

    assert payload["backend"] == "ipopt"
    assert payload["problem_name"] == "electrolyte_lle_eos"
    assert payload["derivative_backend"] == "cppad_implicit"
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
    assert payload["solver_accepted"] is True
    assert payload["accepted"] is True
    assert payload["status"] == "accepted"
    assert np.asarray(payload["variables"], dtype=float).shape == (3,)
    assert np.asarray(payload["constraints"], dtype=float).shape == (7,)
    assert np.asarray(payload["phase_amounts"], dtype=float).shape == (2, 4)
    assert np.asarray(payload["phase_volumes"], dtype=float).shape == (2,)
    assert payload["postsolve"]["derivative_backend"] == "cppad_implicit"
    assert payload["postsolve"]["charge_balance_norm"] <= 1.0e-8
    assert payload["postsolve"]["material_balance_norm"] <= 1.0e-8
    assert payload["postsolve"]["ln_fugacity_consistency_norm"] <= 1.0e-6
    assert payload["postsolve"]["phase_distance"] >= 0.1

    phase_compositions = np.asarray(payload["postsolve"]["phase_compositions"], dtype=float)
    phase_amounts = np.asarray(payload["phase_amounts"], dtype=float)
    phase_volumes = np.asarray(payload["phase_volumes"], dtype=float)
    route_densities = phase_amounts.sum(axis=1) / phase_volumes
    for composition, route_density in zip(phase_compositions, route_densities, strict=True):
        liquid_density = mix.state(T=298.15, P=1.0e5, x=composition, phase="liq").density()
        vapor_density = mix.state(T=298.15, P=1.0e5, x=composition, phase="vap").density()
        assert route_density == pytest.approx(liquid_density, rel=1.0e-10)
        assert route_density / vapor_density > 100.0


def test_electrolyte_lle_route_uses_exact_hessian_when_requested() -> None:
    mix, feed = _ascani_electrolyte_mixture()
    payload = _core._native_electrolyte_lle_eos_route_result(
        mix._native,
        298.15,
        1.0e5,
        feed,
        500,
        1.0e-8,
        0.0,
        "exact",
        20,
        1.0e-8,
        1.0e-3,
        1.0e-7,
        1.0e-6,
        0.1,
        None,
    )

    if not payload["compiled"]:
        assert payload["status"] == "ipopt_dependency_required"
        return

    assert payload["ran"] is True
    assert payload["solver_accepted"] is True
    assert payload["accepted"] is True
    assert payload["status"] == "accepted"
    assert payload["hessian_approximation"] == "exact"
    assert payload["exact_hessian_available"] is True
    assert payload["hessian_backend"] != "limited-memory"
    assert payload["eval_h_calls"] > 0


def test_electrolyte_lle_exact_hessian_dilute_salt_route_keeps_callback_finite() -> None:
    formula = np.asarray(
        [0.6732574103166201, 0.0354880934611478, 0.2323336121370503, 0.05892088408518178],
        dtype=float,
    )
    neutrals = formula[:3] / float(np.sum(formula[:3]))
    salt = 1.0e-6
    formula = np.asarray([*(neutrals * (1.0 - salt)), salt], dtype=float)
    feed = np.asarray([formula[0], formula[1], formula[2], formula[3], formula[3]], dtype=float)
    feed = feed / float(np.sum(feed))
    mix = epcsaft.ePCSAFTMixture.from_dataset(
        "2026_Khudaida",
        ["H2O", "Ethanol", "Butanol", "Na+", "Cl-"],
        feed,
        303.15,
    )

    payload = _core._native_electrolyte_lle_eos_route_result(
        mix._native,
        303.15,
        1.0e5,
        feed.tolist(),
        20,
        1.0e-8,
        0.0,
        "exact",
        20,
        1.0e-7,
        1.0e-3,
        1.0e-7,
        1.0e-6,
        0.1,
        None,
    )

    if not payload["compiled"]:
        assert payload["status"] == "ipopt_dependency_required"
        return

    assert payload["ran"] is True
    assert payload["hessian_approximation"] == "exact"
    assert payload["exact_hessian_available"] is True
    assert payload["eval_h_calls"] > 0
    assert payload["solver_status"] != "invalid_number_detected"
    assert payload["last_callback_exception"] == ""


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
    jacobian = _dense_jacobian_from_sparse_contract(payload)
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
    assert payload["jacobian_nonzero_count"] == 28
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
    jacobian = _dense_jacobian_from_sparse_contract(payload)
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
    assert payload["jacobian_nonzero_count"] == 30
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
    jacobian = _dense_jacobian_from_sparse_contract(payload)
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
    assert payload["jacobian_nonzero_count"] == 52
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
        "limited-memory",
        20,
        1.0e-7,
        1.0e-5,
        1.0e-7,
        1.0e-4,
        None,
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
        ("_native_neutral_bubble_p_eos_route_result", "neutral_bubble_p_eos"),
        ("_native_neutral_dew_p_eos_route_result", "neutral_dew_p_eos"),
    ],
)
def test_neutral_fixed_temperature_pressure_route_uses_exact_hessian_when_requested(
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
        "exact",
        20,
        1.0e-7,
        1.0e-5,
        1.0e-7,
        1.0e-4,
        None,
    )

    if not payload["compiled"]:
        assert payload["status"] == "ipopt_dependency_required"
        return

    assert payload["ran"] is True
    assert payload["problem_name"] == problem_name
    assert payload["hessian_approximation"] == "exact"
    assert payload["exact_hessian_available"] is True
    assert payload["hessian_backend"] != "limited-memory"
    assert payload["eval_h_calls"] > 0
    assert payload["solver_status"] != "invalid_number_detected"
    assert payload["last_callback_exception"] == ""


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
        "limited-memory",
        20,
        1.0e-7,
        1.0e-5,
        1.0e-7,
        1.0e-4,
        None,
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
def test_neutral_fixed_pressure_temperature_route_uses_exact_hessian_when_requested(
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
        "exact",
        20,
        1.0e-7,
        1.0e-5,
        1.0e-7,
        1.0e-4,
        None,
    )

    if not payload["compiled"]:
        assert payload["status"] == "ipopt_dependency_required"
        return

    assert payload["ran"] is True
    assert payload["problem_name"] == problem_name
    assert payload["hessian_approximation"] == "exact"
    assert payload["exact_hessian_available"] is True
    assert payload["hessian_backend"] != "limited-memory"
    assert payload["eval_h_calls"] > 0
    assert payload["solver_status"] != "invalid_number_detected"
    assert payload["last_callback_exception"] == ""


def test_electrolyte_bubble_pressure_route_uses_exact_hessian_when_requested() -> None:
    mix = _ionic_mixture()
    payload = _core._native_electrolyte_bubble_p_eos_route_result(
        mix._native,
        298.15,
        [0.9998, 1.0e-4, 1.0e-4],
        30,
        1.0e-8,
        "exact",
        20,
        1.0e-7,
        1.0e-5,
        1.0e-7,
        1.0e-4,
        1.0e-7,
        None,
    )

    if not payload["compiled"]:
        assert payload["status"] == "ipopt_dependency_required"
        return

    assert payload["ran"] is True
    assert payload["problem_name"] == "electrolyte_bubble_p_eos"
    assert payload["hessian_approximation"] == "exact"
    assert payload["exact_hessian_available"] is True
    assert payload["hessian_backend"] != "limited-memory"
    assert payload["eval_h_calls"] > 0
    assert payload["solver_status"] != "invalid_number_detected"
    assert payload["last_callback_exception"] == ""


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
        "limited-memory",
        20,
        1.0e-7,
        1.0e-5,
        1.0e-7,
        1.0e-4,
        None,
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
        "limited-memory",
        20,
        1.0e-7,
        1.0e-5,
        1.0e-7,
        1.0e-4,
        None,
        linear_solver="mumps",
        acceptable_tolerance=9.0e-7,
        constraint_violation_tolerance=8.0e-8,
        dual_infeasibility_tolerance=7.0e-8,
        complementarity_tolerance=6.0e-8,
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


def test_reactive_two_phase_eos_contract_uses_conserved_balances_and_standard_potentials() -> None:
    mix = _neutral_binary_mixture()
    temperature = 300.0
    phase_amounts = [
        np.asarray([0.1, 0.4], dtype=float),
        np.asarray([0.2, 0.3], dtype=float),
    ]
    volumes = [float(phase_amounts[0].sum() / 80.0), float(phase_amounts[1].sum() / 120.0)]
    species_totals = phase_amounts[0] + phase_amounts[1]
    target_pressure = mix.state(
        T=temperature,
        rho=phase_amounts[0].sum() / volumes[0],
        x=phase_amounts[0] / phase_amounts[0].sum(),
        phase="liquid",
    ).pressure()

    payload = _core._native_reactive_two_phase_eos_nlp_contract(
        mix._native,
        temperature,
        target_pressure,
        [phase.tolist() for phase in phase_amounts],
        volumes,
        1,
        [1.0, 1.0],
        [float(species_totals.sum())],
        1,
        [-1.0, 1.0],
        [float(np.log(3.0))],
    )
    phase_system = _core._native_eos_phase_system(
        mix._native,
        temperature,
        target_pressure,
        [phase.tolist() for phase in phase_amounts],
        volumes,
        species_totals.tolist(),
    )

    standard_mu = np.asarray(payload["standard_mu_rt"], dtype=float)
    base_gradient = np.asarray(phase_system["gradient"], dtype=float)
    expected_gradient = base_gradient.copy()
    expected_gradient[:2] += standard_mu
    expected_gradient[3:5] += standard_mu
    expected_objective = phase_system["objective"] + float(standard_mu @ species_totals)
    contract_jacobian = _dense_jacobian_from_sparse_contract(payload)
    phase_system_jacobian = np.asarray(phase_system["constraint_jacobian_row_major"], dtype=float).reshape(4, 6)

    assert payload["problem_name"] == "reactive_two_phase_eos"
    assert payload["derivative_backend"] == "analytic_cppad"
    assert payload["phase_count"] == 2
    assert payload["species_count"] == 2
    assert payload["balance_row_count"] == 1
    assert payload["reaction_count"] == 1
    assert payload["variable_count"] == 6
    assert payload["constraint_count"] == 3
    assert payload["jacobian_nonzero_count"] == 10
    assert payload["objective_at_initial"] == pytest.approx(expected_objective, rel=1.0e-12, abs=1.0e-10)
    assert payload["gradient_at_initial"] == pytest.approx(expected_gradient, rel=1.0e-12, abs=1.0e-10)
    assert payload["constraints_at_initial"][0] == pytest.approx(0.0, abs=1.0e-12)
    assert payload["constraints_at_initial"][1:] == pytest.approx(
        phase_system["constraints"][2:],
        rel=1.0e-12,
        abs=1.0e-8,
    )
    assert standard_mu[1] - standard_mu[0] == pytest.approx(-np.log(3.0))
    assert contract_jacobian[0] == pytest.approx([1.0, 1.0, 0.0, 1.0, 1.0, 0.0])
    assert contract_jacobian[1:] == pytest.approx(phase_system_jacobian[2:], rel=1.0e-12, abs=1.0e-8)


def test_reactive_two_phase_eos_postsolve_checks_reaction_stationarity() -> None:
    mix = _neutral_binary_mixture()
    temperature = 300.0
    phase_amounts = [
        np.asarray([0.1, 0.4], dtype=float),
        np.asarray([0.2, 0.3], dtype=float),
    ]
    volumes = [float(phase_amounts[0].sum() / 80.0), float(phase_amounts[1].sum() / 120.0)]
    species_totals = phase_amounts[0] + phase_amounts[1]
    target_pressure = mix.state(
        T=temperature,
        rho=phase_amounts[0].sum() / volumes[0],
        x=phase_amounts[0] / phase_amounts[0].sum(),
        phase="liquid",
    ).pressure()

    payload = _core._native_reactive_two_phase_eos_postsolve(
        mix._native,
        temperature,
        target_pressure,
        [phase.tolist() for phase in phase_amounts],
        volumes,
        1,
        [1.0, 1.0],
        [float(species_totals.sum())],
        1,
        [-1.0, 1.0],
        [float(np.log(3.0))],
        1.0e-12,
        1.0e12,
        1.0e-12,
        1.0e-3,
    )

    assert payload["derivative_backend"] == "analytic_cppad"
    assert payload["phase_count"] == 2
    assert payload["species_count"] == 2
    assert payload["balance_row_count"] == 1
    assert payload["reaction_count"] == 1
    assert payload["conserved_balance_norm"] == pytest.approx(0.0, abs=1.0e-12)
    assert np.isfinite(payload["pressure_consistency_norm"])
    assert payload["reaction_stationarity_norm"] > 1.0e-12
    assert payload["phase_distance"] > 1.0e-3
    assert payload["accepted"] is False
    assert payload["rejection_reason"] == "reaction_stationarity"
    assert len(payload["constraints"]) == 3
    assert len(payload["reaction_stationarity_residuals"]) == 2


def test_reactive_two_phase_eos_postsolve_accepts_candidate_under_declared_tolerances() -> None:
    mix = _neutral_binary_mixture()
    temperature = 300.0
    phase_amounts = [
        np.asarray([0.1, 0.4], dtype=float),
        np.asarray([0.2, 0.3], dtype=float),
    ]
    volumes = [float(phase_amounts[0].sum() / 80.0), float(phase_amounts[1].sum() / 120.0)]
    species_totals = phase_amounts[0] + phase_amounts[1]
    target_pressure = mix.state(
        T=temperature,
        rho=phase_amounts[0].sum() / volumes[0],
        x=phase_amounts[0] / phase_amounts[0].sum(),
        phase="liquid",
    ).pressure()

    payload = _core._native_reactive_two_phase_eos_postsolve(
        mix._native,
        temperature,
        target_pressure,
        [phase.tolist() for phase in phase_amounts],
        volumes,
        1,
        [1.0, 1.0],
        [float(species_totals.sum())],
        1,
        [-1.0, 1.0],
        [float(np.log(3.0))],
        1.0e-12,
        1.0e12,
        1.0e12,
        1.0e-3,
    )

    assert payload["derivative_backend"] == "analytic_cppad"
    assert payload["accepted"] is True
    assert payload["rejection_reason"] == "accepted"
    assert payload["conserved_balance_norm"] == pytest.approx(0.0, abs=1.0e-12)
    assert payload["phase_distance"] > 1.0e-3


def test_reactive_electrolyte_two_phase_eos_postsolve_rejects_phase_charge_imbalance() -> None:
    species = ["A", "B", "C+", "D-"]
    feed = np.asarray([0.5, 0.3, 0.1, 0.1], dtype=float)
    mix = epcsaft.ePCSAFTMixture.from_params(
        {
            "MW": np.asarray([18.0e-3, 74.0e-3, 23.0e-3, 35.5e-3]),
            "m": np.asarray([1.1, 1.4, 1.0, 1.0]),
            "s": np.asarray([3.0, 3.4, 3.0, 3.0]),
            "e": np.asarray([180.0, 220.0, 150.0, 150.0]),
            "k_ij": np.zeros((4, 4)),
            "z": np.asarray([0.0, 0.0, 1.0, -1.0]),
            "dielc": np.asarray([80.0, 12.0, 1.0, 1.0]),
        },
        species=species,
    )
    phase_amounts = [
        np.asarray([0.2, 0.1, 0.08, 0.02], dtype=float),
        feed - np.asarray([0.2, 0.1, 0.08, 0.02], dtype=float),
    ]
    volumes = [float(phase_amounts[0].sum() / 80.0), float(phase_amounts[1].sum() / 120.0)]

    payload = _core._native_reactive_electrolyte_two_phase_eos_postsolve(
        mix._native,
        298.15,
        1.013e5,
        [phase.tolist() for phase in phase_amounts],
        volumes,
        3,
        [1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0],
        [float(feed[0] + feed[1]), float(feed[2]), float(feed[3])],
        1,
        [-1.0, 1.0, 0.0, 0.0],
        [float(np.log(0.2))],
        1.0e-12,
        1.0e12,
        1.0e12,
        1.0e-3,
    )

    assert payload["derivative_backend"] == "analytic_cppad"
    assert payload["accepted"] is False
    assert payload["rejection_reason"] == "charge_balance"
    assert payload["conserved_balance_norm"] == pytest.approx(0.0, abs=1.0e-12)
    assert payload["charge_balance_norm"] > 1.0e-12


def test_reactive_two_phase_eos_route_result_uses_native_ipopt_gate() -> None:
    mix = _neutral_binary_mixture()
    temperature = 300.0
    phase_amounts = [
        np.asarray([0.1, 0.4], dtype=float),
        np.asarray([0.2, 0.3], dtype=float),
    ]
    volumes = [float(phase_amounts[0].sum() / 80.0), float(phase_amounts[1].sum() / 120.0)]
    species_totals = phase_amounts[0] + phase_amounts[1]
    target_pressure = mix.state(
        T=temperature,
        rho=phase_amounts[0].sum() / volumes[0],
        x=phase_amounts[0] / phase_amounts[0].sum(),
        phase="liquid",
    ).pressure()

    payload = _core._native_reactive_two_phase_eos_route_result(
        mix._native,
        temperature,
        target_pressure,
        [phase.tolist() for phase in phase_amounts],
        volumes,
        1,
        [1.0, 1.0],
        [float(species_totals.sum())],
        1,
        [-1.0, 1.0],
        [float(np.log(3.0))],
        10,
        1.0e-8,
        0.0,
        "limited-memory",
        20,
        1.0e-8,
        1.0e-6,
        1.0e-6,
        1.0e-3,
        None,
        linear_solver="mumps",
        acceptable_tolerance=9.0e-7,
        constraint_violation_tolerance=8.0e-8,
        dual_infeasibility_tolerance=7.0e-8,
        complementarity_tolerance=6.0e-8,
    )

    assert payload["backend"] == "ipopt"
    assert payload["problem_name"] == "reactive_two_phase_eos"
    assert payload["derivative_backend"] == "analytic_cppad"
    assert payload["exact_gradient_required"] is True
    assert payload["exact_jacobian_required"] is True
    assert payload["phase_count"] == 2
    assert payload["species_count"] == 2
    assert payload["balance_row_count"] == 1
    assert payload["reaction_count"] == 1
    standard_mu = np.asarray(payload["standard_mu_rt"], dtype=float)
    assert standard_mu[1] - standard_mu[0] == pytest.approx(-np.log(3.0))
    if not payload["compiled"]:
        assert payload["ran"] is False
        assert payload["accepted"] is False
        assert payload["status"] == "ipopt_dependency_required"
        assert payload["postsolve"]["accepted"] is False
        return

    assert payload["ran"] is True
    assert payload["status"] in {"accepted", "solver_rejected"}


def test_reactive_two_phase_eos_route_uses_exact_hessian_when_requested() -> None:
    mix = _neutral_binary_mixture()
    payload = _core._native_reactive_two_phase_eos_route_result(
        mix._native,
        300.0,
        1.0e5,
        [[0.1, 0.4], [0.2, 0.3]],
        [0.005, 0.004],
        1,
        [1.0, 1.0],
        [1.0],
        1,
        [-1.0, 1.0],
        [float(np.log(3.0))],
        10,
        1.0e-8,
        0.0,
        "exact",
        20,
        1.0e-8,
        1.0e-6,
        1.0e-6,
        1.0e-3,
        None,
    )

    if not payload["compiled"]:
        assert payload["status"] == "ipopt_dependency_required"
        return

    assert payload["ran"] is True
    assert payload["hessian_approximation"] == "exact"
    assert payload["exact_hessian_available"] is True
    assert payload["hessian_backend"] != "limited-memory"
    assert payload["eval_h_calls"] > 0
    assert payload["solver_status"] != "invalid_number_detected"
    assert payload["last_callback_exception"] == ""


def test_reactive_lle_eos_route_builder_owns_canonical_initial_point() -> None:
    mix = _neutral_binary_mixture()
    temperature = 300.0
    target_pressure = 1.0e5
    feed = np.asarray([0.3, 0.7], dtype=float)

    contract = _core._native_reactive_lle_eos_nlp_contract(
        mix._native,
        temperature,
        target_pressure,
        feed.tolist(),
        1,
        [1.0, 1.0],
        [1.0],
        1,
        [-1.0, 1.0],
        [float(np.log(3.0))],
    )
    initial = np.asarray(contract["initial_point"], dtype=float)
    first_amounts = np.exp(initial[:2])
    second_amounts = np.exp(initial[2:4])
    first = first_amounts / np.sum(first_amounts)
    second = second_amounts / np.sum(second_amounts)

    assert contract["problem_name"] == "reactive_liquid_root_eos"
    assert contract["derivative_backend"] == "cppad_implicit"
    assert contract["density_backend"] == "liquid_pressure_root"
    assert contract["variable_model"] == "log_phase_species_amounts"
    assert contract["variable_count"] == 2 * contract["species_count"]
    assert contract["constraint_count"] == 2
    assert contract["jacobian_nonzero_count"] == 8
    assert contract["balance_row_count"] == 1
    assert contract["reaction_count"] == 1
    assert np.max(np.abs(first - second)) > 1.0e-3
    assert contract["constraint_lower_bounds"][0] == pytest.approx(0.0)
    assert contract["constraint_upper_bounds"][0] == pytest.approx(0.0)
    assert contract["constraint_lower_bounds"][-1] == pytest.approx(1.0e-8)
    assert contract["constraint_upper_bounds"][-1] > 1.0e6
    assert contract["constraints_at_initial"][0] == pytest.approx(0.0)
    assert contract["constraints_at_initial"][-1] >= contract["constraint_lower_bounds"][-1]
    assert np.all(np.asarray(contract["variable_upper_bounds"], dtype=float) < 50.0)
    assert np.asarray(contract["jacobian_values_at_initial"], dtype=float).shape == (
        contract["jacobian_nonzero_count"],
    )
    assert np.count_nonzero(np.asarray(contract["jacobian_values_at_initial"], dtype=float)) > 0

    payload = _core._native_reactive_lle_eos_route_result(
        mix._native,
        temperature,
        target_pressure,
        feed.tolist(),
        1,
        [1.0, 1.0],
        [1.0],
        1,
        [-1.0, 1.0],
        [float(np.log(3.0))],
        10,
        1.0e-8,
        0.0,
        "limited-memory",
        20,
        1.0e-8,
        1.0e-3,
        1.0e-8,
        1.0e-3,
        1.0e-12,
        [0],
        [],
        None,
    )

    assert payload["backend"] == "ipopt"
    assert payload["problem_name"] == "reactive_liquid_root_eos"
    assert payload["derivative_backend"] == "cppad_implicit"
    assert payload["exact_gradient_required"] is True
    assert payload["exact_jacobian_required"] is True
    assert payload["phase_count"] == 2
    assert payload["species_count"] == 2
    assert payload["balance_row_count"] == 1
    assert payload["reaction_count"] == 1
    if not payload["compiled"]:
        assert payload["ran"] is False
        assert payload["accepted"] is False
        assert payload["status"] == "ipopt_dependency_required"
        return

    assert payload["ran"] is True
    assert payload["status"] in {"accepted", "solver_rejected", "postsolve_rejected"}
    if payload["status"] != "solver_rejected":
        assert payload["postsolve"]["density_backend"] == "liquid_pressure_root"


def test_reactive_lle_eos_route_uses_exact_hessian_when_requested() -> None:
    mix = _neutral_binary_mixture()
    payload = _core._native_reactive_lle_eos_route_result(
        mix._native,
        300.0,
        1.0e5,
        [0.3, 0.7],
        1,
        [1.0, 1.0],
        [1.0],
        1,
        [-1.0, 1.0],
        [float(np.log(3.0))],
        10,
        1.0e-8,
        0.0,
        "exact",
        20,
        1.0e-8,
        1.0e-3,
        1.0e-8,
        1.0e-3,
        1.0e-12,
        [0],
        [],
        None,
    )

    if not payload["compiled"]:
        assert payload["status"] == "ipopt_dependency_required"
        return

    assert payload["ran"] is True
    assert payload["hessian_approximation"] == "exact"
    assert payload["exact_hessian_available"] is True
    assert payload["hessian_backend"] != "limited-memory"
    assert payload["eval_h_calls"] > 0
    assert payload["solver_status"] != "invalid_number_detected"
    assert payload["last_callback_exception"] == ""


def test_reactive_electrolyte_lle_eos_route_builder_uses_liquid_root_residual_route() -> None:
    species = ["A", "B", "C+", "D-"]
    feed = np.asarray([0.535, 0.25, 0.1075, 0.1075], dtype=float)
    mix = epcsaft.ePCSAFTMixture.from_params(
        {
            "MW": np.asarray([18.0e-3, 74.0e-3, 23.0e-3, 35.5e-3]),
            "m": np.asarray([1.1, 1.4, 1.0, 1.0]),
            "s": np.asarray([3.0, 3.4, 3.0, 3.0]),
            "e": np.asarray([180.0, 220.0, 150.0, 150.0]),
            "k_ij": np.zeros((4, 4)),
            "z": np.asarray([0.0, 0.0, 1.0, -1.0]),
            "dielc": np.asarray([80.0, 12.0, 1.0, 1.0]),
        },
        species=species,
    )
    charges = np.asarray(mix.parameters["z"], dtype=float)
    balance_matrix = [
        1.0,
        1.0,
        0.0,
        0.0,
        0.0,
        0.0,
        1.0,
        0.0,
        0.0,
        0.0,
        0.0,
        1.0,
    ]
    totals = [float(feed[0] + feed[1]), float(feed[2]), float(feed[3])]
    reaction = [-1.0, 1.0, 0.0, 0.0]

    contract = _core._native_reactive_electrolyte_lle_eos_nlp_contract(
        mix._native,
        298.15,
        1.013e5,
        feed.tolist(),
        3,
        balance_matrix,
        totals,
        1,
        reaction,
        [float(np.log(0.2))],
    )
    initial = np.asarray(contract["initial_point"], dtype=float)
    first = np.exp(initial[:4])
    second = np.exp(initial[4:8])

    assert contract["problem_name"] == "reactive_liquid_root_eos"
    assert contract["derivative_backend"] == "cppad_implicit"
    assert contract["density_backend"] == "liquid_pressure_root"
    assert contract["variable_model"] == "log_phase_species_amounts"
    assert contract["variable_count"] == 2 * contract["species_count"]
    assert contract["constraint_count"] == 6
    assert contract["jacobian_nonzero_count"] == 20
    assert contract["balance_row_count"] == 3
    assert contract["reaction_count"] == 1
    assert contract["constraint_lower_bounds"][:3] == pytest.approx([0.0, 0.0, 0.0])
    assert contract["constraint_upper_bounds"][:3] == pytest.approx([0.0, 0.0, 0.0])
    assert contract["constraint_lower_bounds"][3:5] == pytest.approx([0.0, 0.0])
    assert contract["constraint_upper_bounds"][3:5] == pytest.approx([0.0, 0.0])
    assert contract["constraint_lower_bounds"][-1] == pytest.approx(1.0e-8)
    assert contract["constraint_upper_bounds"][-1] > 1.0e6
    assert contract["constraints_at_initial"][:5] == pytest.approx([0.0, 0.0, 0.0, 0.0, 0.0])
    assert contract["constraints_at_initial"][-1] >= contract["constraint_lower_bounds"][-1]
    assert np.all(np.asarray(contract["variable_upper_bounds"], dtype=float) < 50.0)
    assert np.asarray(contract["jacobian_values_at_initial"], dtype=float).shape == (
        contract["jacobian_nonzero_count"],
    )
    assert np.count_nonzero(np.asarray(contract["jacobian_values_at_initial"], dtype=float)) > 0
    assert np.dot(first, charges) == pytest.approx(0.0, abs=1.0e-14)
    assert np.dot(second, charges) == pytest.approx(0.0, abs=1.0e-14)

    payload = _core._native_reactive_electrolyte_lle_eos_route_result(
        mix._native,
        298.15,
        1.013e5,
        feed.tolist(),
        3,
        balance_matrix,
        totals,
        1,
        reaction,
        [float(np.log(0.2))],
        10,
        1.0e-8,
        0.0,
        "limited-memory",
        20,
        1.0e-8,
        1.0e-3,
        1.0e-8,
        1.0e-3,
        1.0e-12,
        [0],
        [],
        None,
        linear_solver="mumps",
        acceptable_tolerance=9.0e-7,
        constraint_violation_tolerance=8.0e-8,
        dual_infeasibility_tolerance=7.0e-8,
        complementarity_tolerance=6.0e-8,
    )

    assert payload["backend"] == "ipopt"
    assert payload["problem_name"] == "reactive_liquid_root_eos"
    assert payload["derivative_backend"] == "cppad_implicit"
    assert payload["balance_row_count"] == 3
    assert payload["reaction_count"] == 1
    if not payload["compiled"]:
        assert payload["ran"] is False
        assert payload["accepted"] is False
        assert payload["status"] == "ipopt_dependency_required"
        return

    assert payload["ran"] is True
    assert payload["status"] in {"accepted", "solver_rejected", "postsolve_rejected"}
    if payload["status"] != "solver_rejected":
        assert payload["postsolve"]["density_backend"] == "liquid_pressure_root"


def test_reactive_electrolyte_lle_eos_route_uses_exact_hessian_when_requested() -> None:
    species = ["A", "B", "C+", "D-"]
    feed = np.asarray([0.535, 0.25, 0.1075, 0.1075], dtype=float)
    mix = epcsaft.ePCSAFTMixture.from_params(
        {
            "MW": np.asarray([18.0e-3, 74.0e-3, 23.0e-3, 35.5e-3]),
            "m": np.asarray([1.1, 1.4, 1.0, 1.0]),
            "s": np.asarray([3.0, 3.4, 3.0, 3.0]),
            "e": np.asarray([180.0, 220.0, 150.0, 150.0]),
            "k_ij": np.zeros((4, 4)),
            "z": np.asarray([0.0, 0.0, 1.0, -1.0]),
            "dielc": np.asarray([80.0, 12.0, 1.0, 1.0]),
        },
        species=species,
    )
    balance_matrix = [
        1.0,
        1.0,
        0.0,
        0.0,
        0.0,
        0.0,
        1.0,
        0.0,
        0.0,
        0.0,
        0.0,
        1.0,
    ]
    payload = _core._native_reactive_electrolyte_lle_eos_route_result(
        mix._native,
        298.15,
        1.013e5,
        feed.tolist(),
        3,
        balance_matrix,
        [float(feed[0] + feed[1]), float(feed[2]), float(feed[3])],
        1,
        [-1.0, 1.0, 0.0, 0.0],
        [float(np.log(0.2))],
        10,
        1.0e-8,
        0.0,
        "exact",
        20,
        1.0e-8,
        1.0e-3,
        1.0e-8,
        1.0e-3,
        1.0e-12,
        [0],
        [],
        None,
    )

    if not payload["compiled"]:
        assert payload["status"] == "ipopt_dependency_required"
        return

    assert payload["ran"] is True
    assert payload["hessian_approximation"] == "exact"
    assert payload["exact_hessian_available"] is True
    assert payload["hessian_backend"] != "limited-memory"
    assert payload["eval_h_calls"] > 0
    assert payload["solver_status"] != "invalid_number_detected"
    assert payload["last_callback_exception"] == ""


def test_neutral_lle_route_result_records_multistart_seed_attempts_on_failure() -> None:
    mix = _neutral_binary_mixture()
    payload = _core._native_neutral_lle_eos_route_result(
        mix._native,
        298.15,
        1.013e5,
        [0.45, 0.55],
        0,
        1.0e-8,
        0.0,
        "limited-memory",
        2,
        1.0e-8,
        1.0e-8,
        1.0e-8,
        1.0e-3,
        None,
    )

    if not payload["compiled"]:
        assert payload["status"] == "ipopt_dependency_required"
        return

    assert payload["initial_point_strategy"] == "deterministic_multistart"
    assert payload["seed_name"] in {"canonical_shifted_feed", "mirrored_shifted_feed"}
    attempts = payload["seed_attempts"]
    assert len(attempts) >= 2
    assert attempts[0]["seed_name"] == "canonical_shifted_feed"
    assert {attempt["seed_name"] for attempt in attempts} >= {
        "canonical_shifted_feed",
        "mirrored_shifted_feed",
    }
    assert all("status" in attempt for attempt in attempts)
    assert all("iteration_count" in attempt for attempt in attempts)


def test_neutral_bubble_pressure_route_records_multistart_seed_attempts_on_failure() -> None:
    mix = _neutral_binary_mixture()
    payload = _core._native_neutral_bubble_p_eos_route_result(
        mix._native,
        300.0,
        [0.35, 0.65],
        0,
        1.0e-8,
        0.0,
        "limited-memory",
        2,
        1.0e-8,
        1.0e-8,
        1.0e-8,
        1.0e-3,
        None,
    )

    if not payload["compiled"]:
        assert payload["status"] == "ipopt_dependency_required"
        return

    assert payload["initial_point_strategy"] == "deterministic_multistart"
    assert payload["seed_name"] in {
        "canonical_shifted_partner_phase",
        "mirrored_shifted_partner_phase",
    }
    attempts = payload["seed_attempts"]
    assert len(attempts) >= 2
    assert attempts[0]["seed_name"] == "canonical_shifted_partner_phase"
    assert {attempt["seed_name"] for attempt in attempts} >= {
        "canonical_shifted_partner_phase",
        "mirrored_shifted_partner_phase",
    }


def test_neutral_stability_route_records_multistart_seed_attempts_on_failure() -> None:
    mix = _neutral_binary_mixture()
    payload = _core._native_neutral_stability_tpd_route_result(
        mix._native,
        300.0,
        1.0e5,
        [0.3, 0.7],
        "vap",
        "vap",
        0,
        1.0e-8,
        0.0,
        "limited-memory",
        2,
        1.0e-8,
        [],
        None,
    )

    if not payload["compiled"]:
        assert payload["status"] == "ipopt_dependency_required"
        return

    assert payload["initial_point_strategy"] == "deterministic_multistart"
    assert payload["seed_name"] in {"canonical_shifted_feed", "mirrored_shifted_feed"}
    attempts = payload["seed_attempts"]
    assert len(attempts) >= 2
    assert attempts[0]["seed_name"] == "canonical_shifted_feed"
    assert {attempt["seed_name"] for attempt in attempts} >= {
        "canonical_shifted_feed",
        "mirrored_shifted_feed",
    }


def test_electrolyte_lle_route_records_formula_seed_attempts_on_failure() -> None:
    mix, feed = _ascani_electrolyte_mixture()
    payload = _core._native_electrolyte_lle_eos_route_result(
        mix._native,
        298.15,
        1.0e5,
        feed,
        0,
        1.0e-8,
        0.0,
        "limited-memory",
        2,
        1.0e-8,
        1.0e-8,
        1.0e-8,
        1.0e-6,
        0.1,
        None,
    )

    if not payload["compiled"]:
        assert payload["status"] == "ipopt_dependency_required"
        return

    assert payload["initial_point_strategy"] == "deterministic_multistart"
    assert payload["seed_name"] in {"canonical_formula_shift", "mirrored_formula_shift"}
    attempts = payload["seed_attempts"]
    assert len(attempts) >= 2
    assert attempts[0]["seed_name"] == "canonical_formula_shift"
    assert {attempt["seed_name"] for attempt in attempts} >= {
        "canonical_formula_shift",
        "mirrored_formula_shift",
    }


def test_reactive_lle_route_records_multistart_seed_attempts_on_failure() -> None:
    mix = _neutral_binary_mixture()
    payload = _core._native_reactive_lle_eos_route_result(
        mix._native,
        300.0,
        1.0e5,
        [0.3, 0.7],
        1,
        [1.0, 1.0],
        [1.0],
        1,
        [-1.0, 1.0],
        [float(np.log(3.0))],
        0,
        1.0e-8,
        0.0,
        "limited-memory",
        2,
        1.0e-8,
        1.0e-3,
        1.0e-8,
        1.0e-3,
        1.0e-12,
        [0],
        [],
        None,
    )

    if not payload["compiled"]:
        assert payload["status"] == "ipopt_dependency_required"
        return

    assert payload["initial_point_strategy"] == "deterministic_multistart"
    assert payload["seed_name"] in {"canonical_shifted_feed", "mirrored_shifted_feed"}
    attempts = payload["seed_attempts"]
    assert len(attempts) >= 2
    assert attempts[0]["seed_name"] == "canonical_shifted_feed"
    assert {attempt["seed_name"] for attempt in attempts} >= {
        "canonical_shifted_feed",
        "mirrored_shifted_feed",
    }
