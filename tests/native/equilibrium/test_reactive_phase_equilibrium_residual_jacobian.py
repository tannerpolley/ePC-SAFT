from __future__ import annotations

import math

import numpy as np

from epcsaft import _core
from tests.native.equilibrium.test_reactive_phase_equilibrium_residual_surface import (
    _neutral_reactive_lle_mixture,
)


def _request() -> dict[str, object]:
    liq1 = np.asarray([0.05, 0.95], dtype=float)
    liq2 = np.asarray([0.85, 0.15], dtype=float)
    feed = (0.5 * liq1 + 0.5 * liq2).tolist()
    return {
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
        "options": {"min_composition": 1.0e-12, "tolerance": 1.0e-10, "jacobian_backend": "cppad"},
    }


def test_reactive_phase_residual_jacobian_matches_source_perturbation() -> None:
    mix = _neutral_reactive_lle_mixture()
    request = _request()
    payload = _core._evaluate_reactive_phase_equilibrium_residual_native(mix._native, request)
    diagnostics = payload["diagnostics"]

    assert payload["jacobian_backend"] == "cppad_implicit"
    assert diagnostics["jacobian_available"] is True
    assert diagnostics["solved_state_sensitivity_available"] is True
    assert diagnostics["derivative_backend"] == "cppad_implicit"

    variables = np.asarray(payload["variables"], dtype=float)
    residual = np.asarray(payload["residual"], dtype=float)
    jacobian = np.asarray(payload["jacobian_row_major"], dtype=float).reshape(payload["jacobian_shape"])
    assert jacobian.shape == (residual.size, variables.size)
    assert np.all(np.isfinite(jacobian))

    direction = np.asarray([0.17, -0.11, -0.07, 0.13], dtype=float)
    direction /= np.linalg.norm(direction)
    step = 1.0e-6
    plus_request = dict(request)
    plus_request["variables"] = (variables + step * direction).tolist()
    minus_request = dict(request)
    minus_request["variables"] = (variables - step * direction).tolist()
    plus = _core._evaluate_reactive_phase_equilibrium_residual_native(mix._native, plus_request)
    minus = _core._evaluate_reactive_phase_equilibrium_residual_native(mix._native, minus_request)
    source_delta = (np.asarray(plus["residual"], dtype=float) - np.asarray(minus["residual"], dtype=float)) / (
        2.0 * step
    )

    np.testing.assert_allclose(jacobian @ direction, source_delta, rtol=2.0e-4, atol=2.0e-5)


def test_ideal_reaction_standard_state_jacobian_uses_log_mole_fraction_basis() -> None:
    mix = _neutral_reactive_lle_mixture()
    request = _request()
    request["reaction_stoichiometry"] = [-1.0, 2.0]
    request["log_equilibrium_constants"] = [0.0]
    request["reaction_standard_states"] = [1]

    payload = _core._evaluate_reactive_phase_equilibrium_residual_native(mix._native, request)
    variables = np.asarray(payload["variables"], dtype=float)
    jacobian = np.asarray(payload["jacobian_row_major"], dtype=float).reshape(payload["jacobian_shape"])
    x1 = np.asarray(payload["phase1_composition"], dtype=float)
    x2 = np.asarray(payload["phase2_composition"], dtype=float)

    assert jacobian.shape[1] == variables.size
    expected_phase1 = np.asarray([-1.0, 2.0], dtype=float) - x1
    expected_phase2 = np.asarray([-1.0, 2.0], dtype=float) - x2

    np.testing.assert_allclose(jacobian[1, :2], expected_phase1, rtol=1.0e-12, atol=1.0e-12)
    np.testing.assert_allclose(jacobian[1, 2:], np.zeros(2), atol=1.0e-12)
    np.testing.assert_allclose(jacobian[2, :2], np.zeros(2), atol=1.0e-12)
    np.testing.assert_allclose(jacobian[2, 2:], expected_phase2, rtol=1.0e-12, atol=1.0e-12)
