from __future__ import annotations

import epcsaft


def test_equilibrium_capabilities_expose_derivative_policy() -> None:
    policy = epcsaft.capabilities()["equilibrium"]["derivative_policy"]

    assert policy["numerical_derivative_backend_available"] is False
    assert policy["unsupported_derivative_status"] == "not_available"
    assert policy["auto_policy"] == "analytic_or_cppad_or_implicit_where_available_else_not_available"
    assert "numerical_derivative" not in {str(item).lower() for item in policy["accepted_derivative_backends"]}
    assert {
        "analytic",
        "cppad",
        "analytic_implicit",
        "cppad_implicit",
        "not_available",
    }.issubset(set(policy["accepted_derivative_backends"]))
    assert {
        "thermodynamic_backend",
        "solver_backend",
        "derivative_backend",
        "derivative_status",
        "not_available_reason",
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
    assert reactive["derivative_gap_status"] == "implicit_sensitivity_available_for_reaction_constant_response"


def test_reactive_phase_equilibrium_capabilities_state_reaction_scope() -> None:
    reactive = epcsaft.capabilities()["equilibrium"]["reactive_phase_equilibrium"]

    assert reactive["available"] is True
    assert {"reactive_lle", "reactive_electrolyte_lle"}.issubset(set(reactive["methods"]))
    assert reactive["problem_class"] == "ReactivePhaseEquilibriumProblem"
    assert reactive["solver_dependency"] == "ceres"
    assert "same_phase_activity_reaction" in reactive["supported_reaction_scopes"]
    assert "phase_tagged_cross_phase_quotient" in reactive["supported_reaction_scopes"]
    assert reactive["unsupported_reaction_scopes"] == []
    assert reactive["cross_phase_reaction_quotients"]["available"] is True
    assert reactive["cross_phase_reaction_quotients"]["api"] == "ReactionDefinition.phase_stoichiometry"
