from __future__ import annotations

import json

import numpy as np
import pytest

import epcsaft
from epcsaft import _core


def _associating_lle_mixture() -> epcsaft.ePCSAFTMixture:
    return epcsaft.ePCSAFTMixture.from_params(
        {
            "m": np.asarray([1.5255, 2.5303]),
            "s": np.asarray([3.2300, 3.8499]),
            "e": np.asarray([188.90, 278.11]),
            "e_assoc": np.asarray([2899.5, 0.0]),
            "vol_a": np.asarray([0.035176, 0.0]),
            "assoc_scheme": ["2B", None],
            "k_ij": np.asarray([[0.0, 0.051], [0.051, 0.0]]),
        },
        species=["Methanol", "Cyclohexane"],
    )


def _electrolyte_mixture() -> epcsaft.ePCSAFTMixture:
    feed = np.asarray([0.55, 0.40, 0.025, 0.025], dtype=float)
    return epcsaft.ePCSAFTMixture.from_dataset("2022_Ascani", ["H2O", "Butanol", "Na+", "Cl-"], feed, 298.15)


def test_associating_neutral_lle_solves_without_numerical_derivative_derivatives() -> None:
    mix = _associating_lle_mixture()

    result = mix.lle_tp(
        T=298.15,
        P=1.013e5,
        z=[0.5, 0.5],
        options=epcsaft.EquilibriumOptions(max_iterations=240, tolerance=1.0e-8),
    )

    assert result.split_detected is True
    assert result.phase_labels == ["liq1", "liq2"]
    assert result.diagnostics["fugacity_residual_norm"] < 1.0e-8
    assert result.diagnostics["material_balance_error"] < 1.0e-8
    assert result.diagnostics["phase_distance"] > 0.1
    assert result.diagnostics["nonlinear_solver"] == "native_derivative_free_nelder_mead"
    assert result.diagnostics["derivative_backend"] == "not_applicable"
    assert result.diagnostics["derivative_status"] == "not_required"
    payload = json.dumps(result.to_dict(), default=str).lower()
    assert "not_available" not in payload
    assert "numerical_derivative" not in payload


def test_electrolyte_lle_residual_surface_reports_missing_jacobian_without_numeric_derivatives() -> None:
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

    assert payload["jacobian_backend"] == "not_available"
    assert payload["diagnostics"]["jacobian_available"] is False
    assert payload["diagnostics"]["residual_surface"] == "native_electrolyte_lle_transformed_variables"
    assert "numerical_derivative" not in json.dumps(payload, default=str).lower()
