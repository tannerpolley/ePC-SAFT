from __future__ import annotations

import math

import numpy as np
import pytest

from epcsaft import _core
from tests.native.equilibrium.test_electrolyte_lle_residual_surface import _electrolyte_mixture, _initial_request


def test_electrolyte_lle_native_accepted_solve_uses_ceres_jacobian_route() -> None:
    mix = _electrolyte_mixture()
    request = _initial_request(mix)
    request["kind"] = "electrolyte_lle"

    payload = _core._solve_equilibrium_native(mix._native, request)
    diagnostics = payload["diagnostics"]

    assert payload["backend"] == "electrolyte_lle"
    assert payload["split_detected"] is True
    assert len(payload["phases"]) == 2
    assert diagnostics["solver_backend"] == "ceres"
    assert diagnostics["selected_solver_backend"] == "ceres"
    assert diagnostics["solver_method"] == "ceres_trust_region_residual_solve"
    assert diagnostics["ceres_trust_region_strategy"] == "levenberg_marquardt"
    assert diagnostics["ceres_linear_solver"] == "dense_qr"
    assert diagnostics["ceres_termination_type"] in {"convergence", "no_convergence"}
    assert diagnostics["jacobian_backend"] == "cppad_implicit"
    assert diagnostics["derivative_backend"] == "cppad_implicit"
    assert diagnostics["residual_surface_jacobian_backend"] == "cppad_implicit"
    assert diagnostics["residual_surface_derivative_backend"] == "cppad_implicit"
    assert "local_residual_slope" not in str(diagnostics)
    assert diagnostics["jacobian_available"] is True
    assert diagnostics["derivative_available"] is True
    assert diagnostics["ceres_final_cost"] <= diagnostics["ceres_initial_cost"]
    assert diagnostics["solver_residual_norm"] <= request["options"]["tolerance"]
    assert diagnostics["material_balance_error"] <= 1.0e-10
    assert diagnostics["phase_charge_balance_max_abs"] <= 1.0e-8
    assert diagnostics["gibbs_delta"] < 0.0
    assert diagnostics["phase_distance"] > 0.1

    feed = np.asarray(request["z"], dtype=float)
    reconstructed = np.zeros_like(feed)
    for phase in payload["phases"]:
        fraction = float(phase["phase_fraction"])
        composition = np.asarray(phase["composition"], dtype=float)
        assert math.isfinite(float(phase["density"]))
        assert float(phase["density"]) > 0.0
        assert composition.sum() == pytest.approx(1.0, abs=1.0e-12)
        reconstructed += fraction * composition
    np.testing.assert_allclose(reconstructed, feed, atol=1.0e-10)
