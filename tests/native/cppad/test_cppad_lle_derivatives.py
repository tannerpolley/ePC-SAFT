from __future__ import annotations

import json

import numpy as np

import epcsaft
from epcsaft import _core


def _electrolyte_mixture() -> epcsaft.ePCSAFTMixture:
    feed = np.asarray([0.55, 0.40, 0.025, 0.025], dtype=float)
    return epcsaft.ePCSAFTMixture.from_dataset("2022_Ascani", ["H2O", "Butanol", "Na+", "Cl-"], feed, 298.15)


def test_electrolyte_lle_residual_surface_reports_cppad_implicit_jacobian() -> None:
    mix = _electrolyte_mixture()
    aq = np.asarray([0.798324680201737, 0.016320352824141723, 0.09267748348706063, 0.09267748348706063])
    org = np.asarray([0.37006036048879404, 0.6214918588210971, 0.004223890345054407, 0.004223890345054407])
    beta_org = 0.613766575013417
    feed = ((1.0 - beta_org) * aq + beta_org * org).tolist()
    request = {
        "T": 298.15,
        "P": 1.013e5,
        "z": feed,
        "species": mix.species,
        "initial_phases": {"aq": aq.tolist(), "org": org.tolist(), "phase_fraction": beta_org},
        "options": {"max_iterations": 80, "tolerance": 1.0e-8, "min_composition": 1.0e-12},
    }

    payload = _core._evaluate_electrolyte_lle_residual_native(mix._native, request)

    assert payload["jacobian_backend"] == "cppad_implicit"
    assert payload["diagnostics"]["jacobian_available"] is True
    assert payload["diagnostics"]["derivative_available"] is True
    assert len(payload["jacobian_row_major"]) == len(payload["residual"]) * len(payload["variables"])
    assert payload["diagnostics"]["residual_surface"] == "native_electrolyte_lle_transformed_variables"
    assert "numerical" + "_derivative" not in json.dumps(payload, default=str).lower()
