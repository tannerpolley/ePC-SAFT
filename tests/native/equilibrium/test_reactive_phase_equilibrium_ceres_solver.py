from __future__ import annotations

import math

import numpy as np
import pytest

import epcsaft
from epcsaft import _core


def _toy_mixture() -> epcsaft.ePCSAFTMixture:
    return epcsaft.ePCSAFTMixture.from_params(
        {
            "m": np.asarray([1.0, 1.0]),
            "s": np.asarray([3.0, 3.0]),
            "e": np.asarray([200.0, 200.0]),
            "z": np.asarray([0.0, 0.0]),
            "dielc": np.asarray([2.0, 2.0]),
        },
        species=["A", "B"],
    )


def _request(mix: epcsaft.ePCSAFTMixture) -> dict[str, object]:
    feed = np.asarray([0.4, 0.6], dtype=float)
    state = mix.state(T=298.15, P=1.0e5, x=feed, phase="liq")
    ln_phi = state.fugacity_coefficient(natural_log=True)
    log_k = -math.log(feed[0]) - float(ln_phi[0]) + math.log(feed[1]) + float(ln_phi[1])
    return {
        "T": 298.15,
        "P": 1.0e5,
        "z": feed.tolist(),
        "initial_phases": {"liq1": feed.tolist(), "liq2": feed.tolist(), "phase_fraction": 0.5},
        "balance_matrix": [1.0, 1.0],
        "balance_rows": 1,
        "total_vector": [1.0],
        "reaction_stoichiometry": [-1.0, 1.0],
        "reaction_rows": 1,
        "log_equilibrium_constants": [log_k],
        "reaction_standard_states": [1],
        "options": {"max_iterations": 20, "min_composition": 1.0e-12, "tolerance": 1.0e-10},
    }


def test_reactive_phase_native_ceres_solver_reports_coupled_jacobian_route() -> None:
    mix = _toy_mixture()
    payload = _core._solve_reactive_phase_equilibrium_native(mix._native, _request(mix))
    diagnostics = payload["diagnostics"]

    assert payload["backend"] == "reactive_phase_equilibrium"
    assert payload["problem_kind"] == "reactive_phase_equilibrium"
    assert len(payload["phases"]) == 2
    assert diagnostics["solver_backend"] == "ceres"
    assert diagnostics["selected_solver_backend"] == "ceres"
    assert diagnostics["solver_method"] == "ceres_trust_region_coupled_reactive_phase_equilibrium"
    assert diagnostics["ceres_trust_region_strategy"] == "levenberg_marquardt"
    assert diagnostics["ceres_linear_solver"] == "dense_qr"
    assert diagnostics["jacobian_backend"] == "cppad_implicit"
    assert diagnostics["derivative_backend"] == "cppad_implicit"
    assert diagnostics["solved_state_sensitivity_backend"] == "cppad_implicit"
    assert diagnostics["jacobian_available"] is True
    assert diagnostics["solved_state_sensitivity_available"] is True
    assert diagnostics["ceres_accepted_solve"] is True
    assert diagnostics["reaction_and_phase_residuals_share_state"] is True
    assert diagnostics["ceres_final_cost"] <= diagnostics["ceres_initial_cost"] + 1.0e-18
    assert diagnostics["residual_norm"] == pytest.approx(0.0, abs=1.0e-7)
    assert diagnostics["reaction_residual_norm"] == pytest.approx(0.0, abs=1.0e-7)
    assert diagnostics["phase_equilibrium_residual_norm"] == pytest.approx(0.0, abs=1.0e-7)
