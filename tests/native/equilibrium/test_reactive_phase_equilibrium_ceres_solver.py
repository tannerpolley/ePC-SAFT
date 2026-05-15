from __future__ import annotations

import numpy as np
import pytest

import epcsaft
from epcsaft import _core


def _toy_mixture() -> epcsaft.ePCSAFTMixture:
    return epcsaft.ePCSAFTMixture.from_params(
        {
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
        },
        species=["Methanol", "Cyclohexane"],
    )


def _request(mix: epcsaft.ePCSAFTMixture) -> dict[str, object]:
    liq1 = np.asarray([0.11757838279937723, 0.8824216172006228])
    liq2 = np.asarray([0.7985874308392054, 0.20141256916079467])
    beta2 = 0.48813098468607985
    feed = (1.0 - beta2) * liq1 + beta2 * liq2
    return {
        "T": 298.15,
        "P": 1.013e5,
        "z": feed.tolist(),
        "initial_phases": {"liq1": liq1.tolist(), "liq2": liq2.tolist(), "phase_fraction": beta2},
        "balance_matrix": [1.0, 1.0],
        "balance_rows": 1,
        "total_vector": [1.0],
        "reaction_stoichiometry": [-1.0, 1.0],
        "reaction_rows": 1,
        "log_equilibrium_constants": [-0.079259405371],
        "reaction_standard_states": [0],
        "options": {"max_iterations": 80, "min_composition": 1.0e-12, "tolerance": 1.0e-8},
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


def test_reactive_phase_native_ceres_solver_rejects_unusable_solution() -> None:
    mix = _toy_mixture()
    request = _request(mix)
    request["log_equilibrium_constants"] = [10.0]
    request["options"] = dict(request["options"])
    request["options"]["max_iterations"] = 0

    with pytest.raises(_core.NativeSolutionError, match="reactive phase Ceres solve was not accepted"):
        _core._solve_reactive_phase_equilibrium_native(mix._native, request)
