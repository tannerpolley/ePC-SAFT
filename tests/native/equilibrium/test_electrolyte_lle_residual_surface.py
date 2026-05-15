from __future__ import annotations

import math

import numpy as np
import pytest

from epcsaft import _core, ePCSAFTMixture


def _electrolyte_mixture() -> ePCSAFTMixture:
    feed = np.asarray([0.55, 0.40, 0.025, 0.025], dtype=float)
    return ePCSAFTMixture.from_dataset("2022_Ascani", ["H2O", "Butanol", "Na+", "Cl-"], feed, 298.15)


def _initial_request(mix: ePCSAFTMixture) -> dict[str, object]:
    aq = np.asarray([0.798324680201737, 0.016320352824141723, 0.09267748348706063, 0.09267748348706063])
    org = np.asarray([0.37006036048879404, 0.6214918588210971, 0.004223890345054407, 0.004223890345054407])
    beta_org = 0.613766575013417
    feed = ((1.0 - beta_org) * aq + beta_org * org).tolist()
    return {
        "T": 298.15,
        "P": 1.013e5,
        "z": feed,
        "species": mix.species,
        "initial_phases": {"aq": aq.tolist(), "org": org.tolist(), "phase_fraction": beta_org},
        "options": {"max_iterations": 80, "tolerance": 1.0e-8, "min_composition": 1.0e-12},
    }


def test_electrolyte_lle_residual_surface_returns_transformed_payload() -> None:
    mix = _electrolyte_mixture()
    payload = _core._evaluate_electrolyte_lle_residual_native(mix._native, _initial_request(mix))
    diagnostics = payload["diagnostics"]

    assert payload["variable_model"] == "ascani_transformed_salt_pairs"
    assert len(payload["variables"]) == 3
    assert len(payload["lower_bounds"]) == len(payload["variables"])
    assert len(payload["upper_bounds"]) == len(payload["variables"])
    assert payload["jacobian_shape"] == (len(payload["residual"]), len(payload["variables"]))
    assert len(payload["jacobian_row_major"]) == len(payload["residual"]) * len(payload["variables"])
    assert len(payload["gradient"]) == len(payload["variables"])
    assert payload["jacobian_backend"] == "cppad_implicit"
    assert diagnostics["jacobian_available"] is True
    assert diagnostics["residual_surface"] == "native_electrolyte_lle_transformed_variables"
    assert diagnostics["residual_blocks"] == "phase_equilibrium,material_balance"
    assert diagnostics["phase_equilibrium_residual_size"] == 3
    assert diagnostics["material_balance_residual_size"] == 4
    assert diagnostics["residual_size"] == 7
    assert math.isfinite(payload["objective"])
    assert payload["objective"] >= 0.0
    assert payload["material_balance_error"] <= 1.0e-10
    assert payload["charge_balance_error"] <= 1.0e-8
    assert payload["phase_distance"] > 0.1

    charges = np.asarray(mix.parameters["z"], dtype=float)
    assert float(np.dot(payload["aq_composition"], charges)) == pytest.approx(0.0, abs=1.0e-8)
    assert float(np.dot(payload["org_composition"], charges)) == pytest.approx(0.0, abs=1.0e-8)


def test_electrolyte_lle_residual_surface_reuses_explicit_variables() -> None:
    mix = _electrolyte_mixture()
    request = _initial_request(mix)
    initial = _core._evaluate_electrolyte_lle_residual_native(mix._native, request)
    request = dict(request)
    request["variables"] = initial["variables"]
    repeated = _core._evaluate_electrolyte_lle_residual_native(mix._native, request)

    np.testing.assert_allclose(repeated["residual"], initial["residual"], atol=1.0e-10)
    np.testing.assert_allclose(repeated["jacobian_row_major"], initial["jacobian_row_major"], atol=1.0e-10)
    np.testing.assert_allclose(repeated["aq_composition"], initial["aq_composition"], atol=1.0e-12)
    np.testing.assert_allclose(repeated["org_composition"], initial["org_composition"], atol=1.0e-12)
