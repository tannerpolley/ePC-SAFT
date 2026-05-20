from __future__ import annotations

from epcsaft.equilibrium_core.native_results import native_route_diagnostics


def test_native_route_diagnostics_normalizes_solver_route_and_seed_contract() -> None:
    route = {
        "compiled": True,
        "adapter_available": True,
        "ran": True,
        "accepted": False,
        "status": "postsolve_rejected",
        "solver_accepted": True,
        "solver_status": "success",
        "application_status": "solve_succeeded",
        "gradient_approximation": "exact",
        "jacobian_approximation": "exact",
        "hessian_approximation": "exact",
        "exact_gradient_required": True,
        "exact_jacobian_required": True,
        "exact_hessian_available": True,
        "seed_name": "mirrored_formula_shift",
        "seed_attempts": [
            {
                "seed_name": "canonical_formula_shift",
                "status": "solver_rejected",
                "solver_accepted": False,
                "accepted": False,
            },
            {
                "seed_name": "mirrored_formula_shift",
                "status": "postsolve_rejected",
                "solver_accepted": True,
                "accepted": False,
            },
        ],
        "postsolve": {
            "accepted": False,
            "rejection_reason": "phase_distance",
        },
    }

    diagnostics = native_route_diagnostics(route)

    assert diagnostics["compiled"] is True
    assert diagnostics["adapter_available"] is True
    assert diagnostics["solver_ran"] is True
    assert diagnostics["route_accepted"] is False
    assert diagnostics["postsolve_accepted"] is False
    assert diagnostics["gradient_is_exact"] is True
    assert diagnostics["jacobian_is_exact"] is True
    assert diagnostics["hessian_is_exact"] is True
    assert diagnostics["exact_derivatives_required"] is True
    assert diagnostics["seed_attempt_count"] == 2
    assert diagnostics["seed_attempt_solver_accepted_count"] == 1
    assert diagnostics["seed_attempt_route_accepted_count"] == 0
    assert diagnostics["seed_attempts"][0]["route_status"] == "solver_rejected"
    assert diagnostics["seed_attempts"][0]["selected_seed"] is False
    assert diagnostics["seed_attempts"][1]["route_status"] == "postsolve_rejected"
    assert diagnostics["seed_attempts"][1]["selected_seed"] is True
    assert diagnostics["seed_attempts"][1]["gradient_is_exact"] is True


def test_native_route_diagnostics_preserves_route_metadata_contract() -> None:
    route = {
        "accepted": True,
        "status": "accepted",
        "variable_model": "log_phase_species_amounts_plus_log_density",
        "density_backend": "explicit_log_density_pressure_constraint",
        "residual_families": [
            "conserved_balance",
            "reaction_stationarity",
            "phase_equilibrium",
            "phase_charge",
        ],
        "constraint_families": [
            "conserved_balance",
            "phase_charge",
            "phase_pressure_consistency",
            "phase_distance",
        ],
        "postsolve": {
            "accepted": True,
        },
    }

    diagnostics = native_route_diagnostics(route)

    assert diagnostics["variable_model"] == "log_phase_species_amounts_plus_log_density"
    assert diagnostics["density_backend"] == "explicit_log_density_pressure_constraint"
    assert diagnostics["residual_families"] == [
        "conserved_balance",
        "reaction_stationarity",
        "phase_equilibrium",
        "phase_charge",
    ]
    assert diagnostics["constraint_families"] == [
        "conserved_balance",
        "phase_charge",
        "phase_pressure_consistency",
        "phase_distance",
    ]
