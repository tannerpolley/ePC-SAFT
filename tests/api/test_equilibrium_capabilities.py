from __future__ import annotations

import epcsaft


def test_equilibrium_capabilities_expose_derivative_policy() -> None:
    policy = epcsaft.capabilities()["equilibrium"]["derivative_policy"]

    assert policy["finite_difference_backend_available"] is False
    assert policy["unsupported_derivative_status"] == "backend_unavailable"
    assert policy["auto_policy"] == "analytic_or_cppad_or_implicit_where_available_else_backend_unavailable"
    assert "finite_difference" not in {str(item).lower() for item in policy["accepted_derivative_backends"]}
    assert {
        "analytic",
        "cppad",
        "analytic_implicit",
        "cppad_implicit",
        "legacy_eigen_forward",
        "backend_unavailable",
    }.issubset(set(policy["accepted_derivative_backends"]))
    assert {
        "thermodynamic_backend",
        "solver_backend",
        "derivative_backend",
        "derivative_status",
        "backend_unavailable_reason",
        "solved_internal_states",
        "derivative_backend_by_block",
        "implicit_sensitivity_blocks",
        "residual_norm_by_block",
        "association_solver_status",
    }.issubset(set(policy["diagnostic_fields"]))


def test_reactive_speciation_capabilities_include_activity_standard_states() -> None:
    reactive = epcsaft.capabilities()["equilibrium"]["reactive_speciation"]

    assert {
        "ideal_mole_fraction",
        "mole_fraction_activity",
        "thermodynamic_activity",
        "concentration",
        "apparent",
    }.issubset(set(reactive["jacobian_auto_supported_standard_states"]))
    assert reactive["derivative_gap_status"] == "activity_derivatives_not_used_by_fixed_point_outer_iteration"
