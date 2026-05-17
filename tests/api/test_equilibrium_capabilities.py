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
        "association_solver_status",
    }.issubset(set(policy["diagnostic_fields"]))


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
    assert reactive["derivative_gap_status"] == "activity_and_concentration_routes_pending_eos_derivative_nlp_blocks"


def test_reactive_phase_equilibrium_capabilities_state_reaction_scope() -> None:
    reactive = epcsaft.capabilities()["equilibrium"]["reactive_phase_equilibrium"]

    assert reactive["available"] is False
    assert reactive["status"] == "route_pending"
    assert reactive["backend"] == "native_ipopt_equilibrium_nlp_required"
    assert {"reactive_lle", "reactive_electrolyte_lle"}.issubset(set(reactive["methods"]))
    assert reactive["problem_class"] == "ReactivePhaseEquilibriumProblem"
    assert "same_phase_activity_reaction" in reactive["supported_reaction_scopes"]
    assert "phase_tagged_cross_phase_quotient" in reactive["supported_reaction_scopes"]
    assert reactive["unsupported_reaction_scopes"] == []
    assert reactive["cross_phase_reaction_quotients"]["available"] is True
    assert reactive["cross_phase_reaction_quotients"]["status"] == "validated_for_pending_ipopt_route"
    assert reactive["cross_phase_reaction_quotients"]["api"] == "ReactionDefinition.phase_stoichiometry"


def test_neutral_lle_capability_requires_native_ipopt_route() -> None:
    capabilities = epcsaft.capabilities()
    ipopt = capabilities["optimizers"]["ipopt"]
    neutral_lle = capabilities["equilibrium"]["neutral_lle_flash"]

    assert neutral_lle["available"] is ipopt["available"]
    assert neutral_lle["status"] == ("available" if ipopt["available"] else "route_pending")
    assert neutral_lle["backend"] == "native_ipopt_equilibrium_nlp"
    assert neutral_lle["methods"] == ["lle_flash", "lle_tp"]


def test_neutral_tp_flash_capability_requires_native_ipopt_route() -> None:
    capabilities = epcsaft.capabilities()
    ipopt = capabilities["optimizers"]["ipopt"]
    neutral_tp = capabilities["equilibrium"]["neutral_tp_flash"]

    assert neutral_tp["available"] is ipopt["available"]
    assert neutral_tp["status"] == ("available" if ipopt["available"] else "route_pending")
    assert neutral_tp["backend"] == "native_ipopt_equilibrium_nlp"
    assert neutral_tp["methods"] == ["tp_flash", "flash_tp"]


def test_stability_capabilities_require_native_ipopt_route() -> None:
    equilibrium = epcsaft.capabilities()["equilibrium"]

    neutral = equilibrium["neutral_stability"]
    assert neutral["available"] is False
    assert neutral["status"] == "route_pending"
    assert neutral["backend"] == "native_ipopt_equilibrium_nlp_required"
    assert neutral["methods"] == ["stability", "stability_tp"]

    electrolyte = equilibrium["electrolyte_stability"]
    assert electrolyte["available"] is False
    assert electrolyte["status"] == "route_pending"
    assert electrolyte["backend"] == "native_ipopt_equilibrium_nlp_required"
    assert electrolyte["methods"] == ["electrolyte_stability", "electrolyte_stability_tp"]


def test_electrolyte_lle_capability_requires_native_ipopt_route() -> None:
    capabilities = epcsaft.capabilities()
    ipopt = capabilities["optimizers"]["ipopt"]
    electrolyte_lle = capabilities["equilibrium"]["electrolyte_lle"]

    assert electrolyte_lle["available"] is ipopt["available"]
    assert electrolyte_lle["status"] == ("available" if ipopt["available"] else "route_pending")
    assert electrolyte_lle["backend"] == "native_ipopt_equilibrium_nlp"
    assert electrolyte_lle["methods"] == ["electrolyte_lle", "electrolyte_lle_tp"]
