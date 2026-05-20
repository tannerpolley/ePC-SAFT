from __future__ import annotations

import json

import numpy as np
import pytest

import epcsaft
from epcsaft.equilibrium_core.native_results import (
    RouteDiagnosticsView,
    native_route_diagnostics,
    native_route_solved_pressure,
    native_route_solved_temperature,
    native_route_summed_phase_amounts,
    neutral_two_phase_payload_to_result,
    raise_native_route_rejected,
)


def _accepted_native_payload() -> dict[str, object]:
    diagnostics = {
        "derivative_backend": "analytic_cppad",
        "rejection_reason": "accepted",
        "material_balance_norm": 1.0e-12,
        "pressure_consistency_norm": 2.0e-8,
        "chemical_potential_consistency_norm": 3.0e-9,
        "ln_fugacity_consistency_norm": 4.0e-9,
        "phase_distance": 0.4,
    }
    return {
        "accepted": True,
        "backend": "native_equilibrium_nlp",
        "problem_kind": "neutral_two_phase_eos",
        "phase_labels": ["phase_0", "phase_1"],
        "stable": False,
        "split_detected": True,
        "diagnostics": diagnostics,
        "phases": [
            {
                "label": "phase_0",
                "composition": [0.7, 0.3],
                "density": 120.0,
                "temperature": 300.0,
                "pressure": 2.0e5,
                "phase_fraction": 0.5,
                "ln_fugacity_coefficient": [-0.01, 0.02],
                "fugacity_coefficient": np.exp([-0.01, 0.02]).tolist(),
                "diagnostics": {"volume": 1.0e-2},
            },
            {
                "label": "phase_1",
                "composition": [0.1, 0.9],
                "density": 120.0,
                "temperature": 300.0,
                "pressure": 2.0e5,
                "phase_fraction": 0.5,
                "ln_fugacity_coefficient": [0.03, -0.04],
                "fugacity_coefficient": np.exp([0.03, -0.04]).tolist(),
                "diagnostics": {"volume": 1.0e-2},
            },
        ],
    }


def test_neutral_native_payload_converts_to_public_equilibrium_result() -> None:
    result = neutral_two_phase_payload_to_result(_accepted_native_payload())

    assert isinstance(result, epcsaft.EquilibriumResult)
    assert result.backend == "native_equilibrium_nlp"
    assert result.problem_kind == "neutral_two_phase_eos"
    assert result.phase_labels == ["phase_0", "phase_1"]
    assert result.stable is False
    assert result.split_detected is True
    assert result.diagnostics["derivative_backend"] == "analytic_cppad"
    assert result.diagnostics["ln_fugacity_consistency_norm"] == pytest.approx(4.0e-9)
    assert len(result.phases) == 2
    assert result.phases[0].label == "phase_0"
    assert np.allclose(result.phases[0].composition, [0.7, 0.3])
    assert np.allclose(result.phases[0].ln_fugacity_coefficient, [-0.01, 0.02])
    assert np.allclose(result.phases[0].fugacity_coefficient, np.exp(result.phases[0].ln_fugacity_coefficient))
    json.dumps(result.to_dict(), allow_nan=False)


def test_neutral_native_payload_rejection_raises_solution_error() -> None:
    payload = _accepted_native_payload()
    payload["accepted"] = False
    payload["rejection_reason"] = "phase_distance"
    payload["diagnostics"] = {"rejection_reason": "phase_distance", "phase_distance": 0.0}
    payload["phases"] = []

    with pytest.raises(epcsaft.SolutionError, match="rejected") as exc_info:
        neutral_two_phase_payload_to_result(payload)

    assert exc_info.value.diagnostics["rejection_reason"] == "phase_distance"


def test_native_route_diagnostics_merges_postsolve_and_solver_metadata() -> None:
    route = {
        "accepted": False,
        "status": "postsolve_rejected",
        "solver_status": "Solve_Succeeded",
        "application_status": "solve_succeeded",
        "backend": "ipopt",
        "problem_name": "electrolyte_lle_eos",
        "adapter_kind": "native_tnlp_adapter",
        "gradient_approximation": "exact",
        "jacobian_approximation": "exact",
        "hessian_approximation": "exact",
        "hessian_backend": "analytic",
        "scaling_method": "user-scaling",
        "linear_solver_requested": "mumps",
        "linear_solver_selected": "mumps",
        "iteration_count": 7,
        "iteration_history_limit": 3,
        "iteration_history_size": 3,
        "variable_scaling_count": 2,
        "constraint_scaling_count": 1,
        "objective_scaling": 0.5,
        "acceptable_tolerance": 1.0e-7,
        "constraint_violation_tolerance": 2.0e-8,
        "dual_infeasibility_tolerance": 3.0e-8,
        "complementarity_tolerance": 4.0e-8,
        "variable_scaling_min": 0.25,
        "variable_scaling_max": 1.0,
        "constraint_scaling_min": 0.5,
        "constraint_scaling_max": 0.5,
        "exact_hessian_available": True,
        "warm_start_requested": True,
        "warm_start_used": True,
        "initial_point_strategy": "deterministic_seed_sweep",
        "seed_name": "mirrored_formula_shift",
        "seed_attempts": [
            {
                "seed_name": "canonical_formula_shift",
                "status": "solver_rejected",
                "solver_accepted": False,
                "accepted": False,
                "iteration_count": 0,
                "objective": 0.0,
            },
            {
                "seed_name": "mirrored_formula_shift",
                "status": "postsolve_rejected",
                "solver_accepted": True,
                "accepted": False,
                "iteration_count": 7,
                "objective": 0.5,
            },
        ],
        "continuation_state": {
            "variables": [0.2, 0.8],
            "bound_lower_multipliers": [0.0, 0.0],
            "bound_upper_multipliers": [0.0, 0.0],
            "constraint_multipliers": [1.0],
        },
        "iteration_history": [
            {
                "iteration": 0,
                "objective": 1.0,
                "primal_infeasibility": 1.0e-2,
            },
            {
                "iteration": 1,
                "objective": 0.5,
                "primal_infeasibility": 1.0e-6,
            },
        ],
        "exact_gradient_required": True,
        "exact_jacobian_required": True,
        "postsolve": {
            "accepted": False,
            "rejection_reason": "phase_distance",
            "phase_distance": 0.0,
        },
    }

    diagnostics = native_route_diagnostics(route)

    assert diagnostics["accepted"] is False
    assert diagnostics["rejection_reason"] == "phase_distance"
    assert diagnostics["route_status"] == "postsolve_rejected"
    assert diagnostics["solver_status"] == "Solve_Succeeded"
    assert diagnostics["application_status"] == "solve_succeeded"
    assert diagnostics["solver_backend"] == "ipopt"
    assert diagnostics["problem_name"] == "electrolyte_lle_eos"
    assert diagnostics["adapter_kind"] == "native_tnlp_adapter"
    assert diagnostics["hessian_approximation"] == "exact"
    assert diagnostics["hessian_backend"] == "analytic"
    assert diagnostics["scaling_method"] == "user-scaling"
    assert diagnostics["linear_solver_requested"] == "mumps"
    assert diagnostics["linear_solver_selected"] == "mumps"
    assert diagnostics["iteration_count"] == 7
    assert diagnostics["iteration_history_limit"] == 3
    assert diagnostics["iteration_history_size"] == 3
    assert diagnostics["variable_scaling_count"] == 2
    assert diagnostics["constraint_scaling_count"] == 1
    assert diagnostics["objective_scaling"] == pytest.approx(0.5)
    assert diagnostics["acceptable_tolerance"] == pytest.approx(1.0e-7)
    assert diagnostics["constraint_violation_tolerance"] == pytest.approx(2.0e-8)
    assert diagnostics["dual_infeasibility_tolerance"] == pytest.approx(3.0e-8)
    assert diagnostics["complementarity_tolerance"] == pytest.approx(4.0e-8)
    assert diagnostics["variable_scaling_min"] == pytest.approx(0.25)
    assert diagnostics["variable_scaling_max"] == pytest.approx(1.0)
    assert diagnostics["constraint_scaling_min"] == pytest.approx(0.5)
    assert diagnostics["constraint_scaling_max"] == pytest.approx(0.5)
    assert diagnostics["exact_hessian_available"] is True
    assert diagnostics["warm_start_requested"] is True
    assert diagnostics["warm_start_used"] is True
    assert diagnostics["initial_point_strategy"] == "deterministic_seed_sweep"
    assert diagnostics["seed_name"] == "mirrored_formula_shift"
    assert diagnostics["seed_attempts"][0]["seed_name"] == "canonical_formula_shift"
    assert diagnostics["seed_attempts"][1]["status"] == "postsolve_rejected"
    assert diagnostics["continuation_state"]["variables"] == pytest.approx([0.2, 0.8])
    assert diagnostics["iteration_history"][1]["objective"] == pytest.approx(0.5)
    assert diagnostics["exact_gradient_required"] is True
    assert diagnostics["exact_jacobian_required"] is True
    assert diagnostics["postsolve_certification"] == {
        "accepted": False,
        "status": "postsolve_rejected",
        "route_accepted": False,
        "postsolve_accepted": False,
        "solver_accepted": False,
        "stability_checked": False,
        "stability_accepted": False,
        "stability_source": "",
        "failure_reason": "phase_distance",
        "active_residuals_reported": False,
        "hard_constraints_reported": True,
        "physical_admissibility_reported": True,
    }


def test_route_diagnostics_view_exposes_stable_interface() -> None:
    route = {
        "accepted": True,
        "status": "accepted",
        "backend": "ipopt",
        "gradient_approximation": "exact",
        "jacobian_approximation": "exact",
        "hessian_approximation": "limited-memory",
        "exact_gradient_required": True,
        "exact_jacobian_required": True,
        "residual_families": ["phase_equilibrium", "material_balance"],
        "constraint_families": ["phase_pressure_consistency", "phase_distance"],
        "seed_name": "canonical_shifted_feed",
        "seed_attempts": [
            {"seed_name": "canonical_shifted_feed", "status": "accepted", "accepted": True},
        ],
        "postsolve": {"accepted": True},
    }

    diagnostics = native_route_diagnostics(route)
    view = RouteDiagnosticsView(diagnostics)
    result = epcsaft.EquilibriumResult(
        backend="native",
        problem_kind="lle",
        phases=(),
        stable=False,
        split_detected=True,
        diagnostics=diagnostics,
    )

    assert epcsaft.RouteDiagnosticsView is RouteDiagnosticsView
    assert view.route_status == "accepted"
    assert view.solver_backend == "ipopt"
    assert view.route_accepted is True
    assert view.postsolve_accepted is True
    assert view.certification_status == "stability_not_checked"
    assert view.postsolve_certification_accepted is False
    assert view.stability_checked is False
    assert view.gradient_is_exact is True
    assert view.jacobian_is_exact is True
    assert view.hessian_is_exact is False
    assert view.exact_derivatives_required is True
    assert view.residual_families == ("phase_equilibrium", "material_balance")
    assert view.constraint_families == ("phase_pressure_consistency", "phase_distance")
    assert view.selected_seed_name == "canonical_shifted_feed"
    assert view.seed_attempt_count == 1
    assert result.route_diagnostics.residual_families == view.residual_families


def test_native_route_diagnostics_marks_tpd_certified_postsolve() -> None:
    route = {
        "accepted": True,
        "status": "accepted",
        "solver_accepted": True,
        "residual_families": ["phase_equilibrium"],
        "constraint_families": ["material_balance", "phase_distance"],
        "postsolve": {
            "accepted": True,
            "ln_fugacity_consistency_norm": 2.0e-10,
            "material_balance_norm": 3.0e-12,
            "phase_distance": 0.4,
        },
        "tpdf_stability": {
            "accepted": True,
            "status": "accepted",
            "feed_unstable": True,
            "final_phases_stable": True,
        },
    }

    diagnostics = native_route_diagnostics(route)
    view = RouteDiagnosticsView(diagnostics)

    assert diagnostics["tpdf_stability"]["accepted"] is True
    assert diagnostics["postsolve_certification"] == {
        "accepted": True,
        "status": "accepted",
        "route_accepted": True,
        "postsolve_accepted": True,
        "solver_accepted": True,
        "stability_checked": True,
        "stability_accepted": True,
        "stability_source": "tpdf_stability",
        "failure_reason": "",
        "active_residuals_reported": True,
        "hard_constraints_reported": True,
        "physical_admissibility_reported": True,
    }
    assert view.certification_status == "accepted"
    assert view.postsolve_certification_accepted is True
    assert view.stability_checked is True


def test_raise_native_route_rejected_uses_shared_diagnostics() -> None:
    route = {
        "accepted": False,
        "status": "solver_rejected",
        "solver_status": "Maximum_Iterations_Exceeded",
        "postsolve": {"accepted": False, "rejection_reason": "solver_rejected"},
    }

    with pytest.raises(epcsaft.SolutionError, match="Native route rejected") as exc_info:
        raise_native_route_rejected(route, "Native route rejected.")

    assert exc_info.value.diagnostics["route_status"] == "solver_rejected"
    assert exc_info.value.diagnostics["solver_status"] == "Maximum_Iterations_Exceeded"
    assert exc_info.value.diagnostics["rejection_reason"] == "solver_rejected"
    assert exc_info.value.route_diagnostics.route_status == "solver_rejected"


def test_native_route_result_helpers_validate_fixed_temperature_payloads() -> None:
    route = {
        "phase_amounts": [[0.7, 0.3], [0.2, 0.8]],
        "variables": [0.7, 0.3, 0.2, 0.8, 1.0e5],
    }

    assert np.allclose(native_route_summed_phase_amounts(route, 2, "bubble_p"), [0.9, 1.1])
    assert native_route_solved_pressure(route, "bubble_p") == pytest.approx(1.0e5)
    assert native_route_solved_temperature({**route, "variables": [0.7, 0.3, 0.2, 0.8, 300.0]}, "bubble_t") == (
        pytest.approx(300.0)
    )

    with pytest.raises(epcsaft.SolutionError, match="positive P"):
        native_route_solved_pressure({"variables": [0.0]}, "bubble_p")
    with pytest.raises(epcsaft.SolutionError, match="positive T"):
        native_route_solved_temperature({"variables": [0.0]}, "bubble_t")
    with pytest.raises(epcsaft.SolutionError, match="positive feed"):
        native_route_summed_phase_amounts({"phase_amounts": [[1.0, -1.0], [1.0, 0.5]]}, 2, "bubble_p")
