from __future__ import annotations

import math

import numpy as np
import pytest

from epcsaft import _core, ePCSAFTMixture


def _neutral_reactive_lle_mixture() -> ePCSAFTMixture:
    params = {
        "MW": np.asarray([32.042e-3, 84.147e-3]),
        "m": np.asarray([1.5255, 2.5303]),
        "s": np.asarray([3.2300, 3.8499]),
        "e": np.asarray([188.90, 278.11]),
        "e_assoc": np.asarray([2899.5, 0.0]),
        "vol_a": np.asarray([0.035176, 0.0]),
        "assoc_scheme": ["2B", None],
        "k_ij": np.asarray([[0.0, 0.051], [0.051, 0.0]]),
        "z": np.asarray([0.0, 0.0]),
        "dielc": np.asarray([33.05, 2.02]),
    }
    return ePCSAFTMixture.from_params(params, species=["Methanol", "Cyclohexane"])


def _electrolyte_reactive_lle_mixture() -> ePCSAFTMixture:
    feed = np.asarray([0.55, 0.40, 0.025, 0.025], dtype=float)
    return ePCSAFTMixture.from_dataset("2022_Ascani", ["H2O", "Butanol", "Na+", "Cl-"], feed, 298.15)


def test_reactive_phase_residual_surface_exposes_single_coupled_state() -> None:
    mix = _neutral_reactive_lle_mixture()
    liq1 = np.asarray([0.05, 0.95], dtype=float)
    liq2 = np.asarray([0.85, 0.15], dtype=float)
    feed = (0.5 * liq1 + 0.5 * liq2).tolist()
    request = {
        "T": 298.15,
        "P": 1.013e5,
        "z": feed,
        "initial_phases": {"liq1": liq1.tolist(), "liq2": liq2.tolist(), "phase_fraction": 0.5},
        "balance_matrix": [1.0, 1.0],
        "balance_rows": 1,
        "total_vector": [1.0],
        "reaction_stoichiometry": [-1.0, 1.0],
        "reaction_rows": 1,
        "log_equilibrium_constants": [math.log(feed[1] / feed[0])],
        "reaction_standard_states": [1],
        "options": {"min_composition": 1.0e-12, "tolerance": 1.0e-10},
    }

    payload = _core._evaluate_reactive_phase_equilibrium_residual_native(mix._native, request)
    diagnostics = payload["diagnostics"]

    assert payload["variable_model"] == "log_phase_species_amounts"
    assert len(payload["variables"]) == 2 * mix.ncomp
    assert payload["jacobian_shape"] == (len(payload["residual"]), len(payload["variables"]))
    assert payload["jacobian_row_major"] == []
    assert payload["jacobian_backend"] == "not_available"
    assert diagnostics["jacobian_available"] is False
    assert diagnostics["ceres_accepted_solve"] is False
    assert diagnostics["residual_surface"] == "native_reactive_phase_equilibrium_coupled_state"
    assert diagnostics["coupling_level"] == "single_native_residual_state"
    assert diagnostics["reaction_and_phase_residuals_share_state"] is True
    assert diagnostics["nonnegative_amounts_enforced_by_transform"] is True
    assert diagnostics["element_balance_residual_size"] == 1
    assert diagnostics["reaction_residual_size_per_phase"] == 1
    assert diagnostics["neutral_phase_equilibrium_residual_size"] == 2
    assert diagnostics["ionic_equilibrium_residual_size"] == 0
    assert diagnostics["phase_charge_residual_size"] == 2
    assert diagnostics["element_balance_norm"] == pytest.approx(0.0, abs=1.0e-10)
    assert diagnostics["phase_distance"] > 0.1
    assert math.isfinite(payload["objective"])
    assert len(payload["reaction_residuals_phase1"]) == 1
    assert len(payload["reaction_residuals_phase2"]) == 1


def test_reactive_phase_reaction_standard_states_define_residual_basis() -> None:
    mix = _neutral_reactive_lle_mixture()
    liq1 = np.asarray([0.05, 0.95], dtype=float)
    liq2 = np.asarray([0.85, 0.15], dtype=float)
    feed = (0.5 * liq1 + 0.5 * liq2).tolist()
    request = {
        "T": 298.15,
        "P": 1.013e5,
        "z": feed,
        "initial_phases": {"liq1": liq1.tolist(), "liq2": liq2.tolist(), "phase_fraction": 0.5},
        "balance_matrix": [1.0, 1.0],
        "balance_rows": 1,
        "total_vector": [1.0],
        "reaction_stoichiometry": [-1.0, 2.0],
        "reaction_rows": 1,
        "log_equilibrium_constants": [0.0],
        "reaction_standard_states": [1],
        "options": {"min_composition": 1.0e-12, "tolerance": 1.0e-10},
    }

    ideal = _core._evaluate_reactive_phase_equilibrium_residual_native(mix._native, request)
    activity_request = dict(request)
    activity_request["reaction_standard_states"] = [0]
    activity = _core._evaluate_reactive_phase_equilibrium_residual_native(mix._native, activity_request)
    concentration_request = dict(request)
    concentration_request["reaction_standard_states"] = [2]
    concentration = _core._evaluate_reactive_phase_equilibrium_residual_native(mix._native, concentration_request)

    x1 = np.asarray(ideal["phase1_composition"], dtype=float)
    ln_phi1 = np.asarray(ideal["phase1_ln_fugacity_coefficient"], dtype=float)
    rho1 = float(ideal["phase1_density"])
    expected_ideal = -math.log(x1[0]) + 2.0 * math.log(x1[1])
    expected_activity = expected_ideal - float(ln_phi1[0]) + 2.0 * float(ln_phi1[1])
    expected_concentration = expected_ideal + math.log(rho1)

    assert ideal["reaction_residuals_phase1"][0] == pytest.approx(expected_ideal, abs=1.0e-12)
    assert activity["reaction_residuals_phase1"][0] == pytest.approx(expected_activity, abs=1.0e-12)
    assert concentration["reaction_residuals_phase1"][0] == pytest.approx(expected_concentration, abs=1.0e-10)
    assert ideal["reaction_residuals_phase1"][0] != pytest.approx(activity["reaction_residuals_phase1"][0])
    assert concentration["reaction_residuals_phase1"][0] != pytest.approx(ideal["reaction_residuals_phase1"][0])


def test_reactive_electrolyte_residual_surface_includes_ion_and_charge_blocks() -> None:
    mix = _electrolyte_reactive_lle_mixture()
    aq = np.asarray([0.798324680201737, 0.016320352824141723, 0.09267748348706063, 0.09267748348706063])
    org = np.asarray([0.37006036048879404, 0.6214918588210971, 0.004223890345054407, 0.004223890345054407])
    beta_org = 0.613766575013417
    feed = ((1.0 - beta_org) * aq + beta_org * org).tolist()
    request = {
        "T": 298.15,
        "P": 1.013e5,
        "z": feed,
        "initial_phases": {"aq": aq.tolist(), "org": org.tolist(), "phase_fraction": beta_org},
        "balance_matrix": np.eye(4, dtype=float).ravel().tolist(),
        "balance_rows": 4,
        "total_vector": feed,
        "reaction_stoichiometry": [-1.0, 1.0, 0.0, 0.0],
        "reaction_rows": 1,
        "log_equilibrium_constants": [0.0],
        "reaction_standard_states": [1],
        "options": {"min_composition": 1.0e-12, "tolerance": 1.0e-10},
    }

    payload = _core._evaluate_reactive_phase_equilibrium_residual_native(mix._native, request)
    diagnostics = payload["diagnostics"]

    assert diagnostics["component_count"] == 4
    assert diagnostics["reaction_count"] == 1
    assert diagnostics["element_balance_residual_size"] == 4
    assert diagnostics["neutral_phase_equilibrium_residual_size"] == 2
    assert diagnostics["ionic_equilibrium_residual_size"] == 1
    assert diagnostics["phase_charge_residual_size"] == 2
    assert diagnostics["element_balance_norm"] == pytest.approx(0.0, abs=1.0e-10)
    assert diagnostics["phase_charge_balance_norm"] == pytest.approx(0.0, abs=1.0e-8)
    assert diagnostics["phase_equilibrium_residual_norm"] >= 0.0
    assert diagnostics["reaction_and_phase_residuals_share_state"] is True
    assert payload["phase_distance"] > 0.1
    assert len(payload["ionic_equilibrium_residuals"]) == 1
