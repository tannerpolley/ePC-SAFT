from __future__ import annotations

import epcsaft


def test_equilibrium_capabilities_expose_derivative_policy() -> None:
    policy = epcsaft.capabilities()["equilibrium"]["derivative_policy"]

    assert policy["unsupported_derivative_behavior"] == "raise"
    assert policy["auto_policy"] == "analytic_or_cppad_or_implicit_else_raise"
    assert "numerical" + "_derivative" not in {str(item).lower() for item in policy["accepted_derivative_backends"]}
    assert {
        "analytic",
        "cppad",
        "analytic_implicit",
        "cppad_implicit",
    }.issubset(set(policy["accepted_derivative_backends"]))
    assert {
        "thermodynamic_backend",
        "solver_backend",
        "derivative_backend",
        "solved_internal_states",
        "derivative_backend_by_block",
        "implicit_sensitivity_blocks",
        "residual_norm_by_block",
        "association_coupling",
    }.issubset(set(policy["diagnostic_fields"]))
    removed_association_status = "association" + "_solver" + "_status"
    assert removed_association_status not in policy["diagnostic_fields"]


def test_reactive_speciation_capabilities_gate_nonideal_standard_states() -> None:
    capabilities = epcsaft.capabilities()
    reactive = capabilities["equilibrium"]["reactive_speciation"]
    ipopt = capabilities["optimizers"]["ipopt"]

    assert reactive["available"] is ipopt["available"]
    assert reactive["backend"] == "native_ipopt_equilibrium_nlp_required"
    assert reactive["status"] == ("available" if ipopt["available"] else "route_pending")
    assert reactive["sweep_available"] is ipopt["available"]
    assert reactive["continuation_state_available"] is ipopt["available"]
    assert reactive["jacobian_auto_supported_standard_states"] == ["ideal_mole_fraction"]
    assert set(reactive["route_gated_standard_states"]) == {
        "mole_fraction_activity",
        "thermodynamic_activity",
        "concentration",
        "apparent",
    }
    assert "derivative_gap_status" not in reactive


def test_reactive_phase_equilibrium_capabilities_state_reaction_scope() -> None:
    reactive = epcsaft.capabilities()["equilibrium"]["reactive_phase_equilibrium"]

    assert reactive["available"] is False
    assert reactive["status"] == "route_pending"
    assert reactive["backend"] == "native_ipopt_equilibrium_nlp_required"
    assert {"reactive_lle", "reactive_electrolyte_lle"}.issubset(set(reactive["methods"]))
    assert reactive["problem_class"] == "ReactivePhaseEquilibriumProblem"
    assert "same_phase_activity_reaction" in reactive["supported_reaction_scopes"]
    assert "phase_tagged_cross_phase_quotient" in reactive["supported_reaction_scopes"]
    removed_scope_field = "unsupported" + "_reaction" + "_scopes"
    assert removed_scope_field not in reactive
    assert reactive["cross_phase_reaction_quotients"]["available"] is True
    assert reactive["cross_phase_reaction_quotients"]["status"] == "validated_for_pending_ipopt_route"
    assert reactive["cross_phase_reaction_quotients"]["api"] == "ReactionDefinition.phase_stoichiometry"
