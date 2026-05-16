from __future__ import annotations

import epcsaft


def test_equilibrium_capabilities_expose_derivative_policy() -> None:
    policy = epcsaft.capabilities()["equilibrium"]["derivative_policy"]

    assert policy["unsupported_derivative_behavior"] == "raise"
    assert policy["auto_policy"] == "analytic_or_cppad_or_implicit_else_raise"
    assert "numerical_derivative" not in {str(item).lower() for item in policy["accepted_derivative_backends"]}
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
        "derivative_status",
        "solved_internal_states",
        "derivative_backend_by_block",
        "implicit_sensitivity_blocks",
        "residual_norm_by_block",
        "association_solver_status",
    }.issubset(set(policy["diagnostic_fields"]))


def test_reactive_speciation_capabilities_include_activity_standard_states() -> None:
    capabilities = epcsaft.capabilities()
    reactive = capabilities["equilibrium"]["reactive_speciation"]
    ipopt = capabilities["optimizers"]["ipopt"]

    assert reactive["available"] is ipopt["available"]
    assert reactive["backend"] == "native_ipopt_equilibrium_nlp_required"
    assert reactive["status"] == ("available" if ipopt["available"] else "route_pending")
    assert reactive["previous_solver_disabled"] == "native_chemical_equilibrium_residual_route"
    assert reactive["sweep_available"] is ipopt["available"]
    assert reactive["continuation_state_available"] is ipopt["available"]
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

    assert reactive["available"] is False
    assert reactive["status"] == "route_pending"
    assert reactive["backend"] == "native_ipopt_equilibrium_nlp_required"
    assert {"reactive_lle", "reactive_electrolyte_lle"}.issubset(set(reactive["methods"]))
    assert reactive["problem_class"] == "ReactivePhaseEquilibriumProblem"
    assert reactive["previous_solver_disabled"] == "ceres_coupled_residual_route"
    assert "same_phase_activity_reaction" in reactive["supported_reaction_scopes"]
    assert "phase_tagged_cross_phase_quotient" in reactive["supported_reaction_scopes"]
    assert reactive["unsupported_reaction_scopes"] == []
    assert reactive["cross_phase_reaction_quotients"]["available"] is True
    assert reactive["cross_phase_reaction_quotients"]["status"] == "validated_for_pending_ipopt_route"
    assert reactive["cross_phase_reaction_quotients"]["api"] == "ReactionDefinition.phase_stoichiometry"


def test_neutral_lle_capability_requires_native_ipopt_route() -> None:
    neutral_lle = epcsaft.capabilities()["equilibrium"]["neutral_lle_flash"]

    assert neutral_lle["available"] is False
    assert neutral_lle["status"] == "route_pending"
    assert neutral_lle["backend"] == "native_ipopt_equilibrium_nlp_required"
    assert neutral_lle["methods"] == ["lle_flash", "lle_tp"]
    assert neutral_lle["previous_solver_disabled"] == "ceres_residual_lle_route"


def test_neutral_tp_flash_capability_requires_native_ipopt_route() -> None:
    neutral_tp = epcsaft.capabilities()["equilibrium"]["neutral_tp_flash"]

    assert neutral_tp["available"] is False
    assert neutral_tp["status"] == "route_pending"
    assert neutral_tp["backend"] == "native_ipopt_equilibrium_nlp_required"
    assert neutral_tp["methods"] == ["tp_flash", "flash_tp"]
    assert neutral_tp["previous_solver_disabled"] == "native_tp_flash_route"


def test_stability_capabilities_require_native_ipopt_route() -> None:
    equilibrium = epcsaft.capabilities()["equilibrium"]

    neutral = equilibrium["neutral_stability"]
    assert neutral["available"] is False
    assert neutral["status"] == "route_pending"
    assert neutral["backend"] == "native_ipopt_equilibrium_nlp_required"
    assert neutral["methods"] == ["stability", "stability_tp"]
    assert neutral["previous_solver_disabled"] == "native_tpd_stability_route"

    electrolyte = equilibrium["electrolyte_stability"]
    assert electrolyte["available"] is False
    assert electrolyte["status"] == "route_pending"
    assert electrolyte["backend"] == "native_ipopt_equilibrium_nlp_required"
    assert electrolyte["methods"] == ["electrolyte_stability", "electrolyte_stability_tp"]
    assert electrolyte["previous_solver_disabled"] == "native_tpd_electrolyte_stability_route"


def test_electrolyte_lle_capability_requires_native_ipopt_route() -> None:
    electrolyte_lle = epcsaft.capabilities()["equilibrium"]["electrolyte_lle"]

    assert electrolyte_lle["available"] is False
    assert electrolyte_lle["status"] == "route_pending"
    assert electrolyte_lle["backend"] == "native_ipopt_equilibrium_nlp_required"
    assert electrolyte_lle["methods"] == ["electrolyte_lle", "electrolyte_lle_tp"]
    assert electrolyte_lle["previous_solver_disabled"] == "ceres_residual_electrolyte_lle_route"
