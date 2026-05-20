"""Registered capability evidence and cheap validation lanes.

This module is intentionally side-effect free so command-line wrappers can
import the registry without importing the native extension.
"""

from __future__ import annotations

from typing import Final

REACTIVE_SPECIATION_STANDARD_STATES: Final[tuple[str, ...]] = (
    "ideal_mole_fraction",
    "mole_fraction_activity",
    "concentration",
)

IPOPT_EQUILIBRIUM_ROUTE_EVIDENCE: Final[tuple[dict[str, object], ...]] = (
    {
        "key": "reactive_speciation",
        "public_routes": tuple(
            f"reactive_speciation:{state}" for state in REACTIVE_SPECIATION_STANDARD_STATES
        ),
        "payload": {
            "backend": "native_ipopt_equilibrium_nlp",
            "sweep_available_from_ipopt": True,
            "activity_output_modes": ("auto", "always", "never"),
            "jacobian_auto_policy": "ideal_analytic_nonideal_cppad_explicit_density_else_raise",
            "jacobian_auto_supported_standard_states": REACTIVE_SPECIATION_STANDARD_STATES,
            "implemented_standard_states": REACTIVE_SPECIATION_STANDARD_STATES,
            "auto_request": "implemented_standard_states_route_to_native_ipopt",
            "solver_backends": ("auto", "ipopt"),
            "explicit_ipopt_request": "implemented_standard_states_route_to_native_ipopt_when_compiled",
            "ipopt_formulation": "thermodynamic_constrained_nlp",
            "ideal_speciation_nlp_available_from_ipopt": True,
            "nonideal_speciation_nlp_available_from_ipopt": True,
            "nonideal_derivative_backend": "cppad_explicit_density",
            "mixed_standard_state_policy": "raise_until_single_objective_is_specified",
        },
    },
    {
        "key": "neutral_tp_flash",
        "public_routes": ("neutral_tp_flash",),
        "payload": {
            "backend": "native_ipopt_equilibrium_nlp",
            "methods": ("tp_flash", "flash_tp"),
            "ipopt_formulation": "thermodynamic_constrained_nlp",
        },
    },
    {
        "key": "neutral_stability",
        "public_routes": ("neutral_stability",),
        "payload": {
            "backend": "native_ipopt_equilibrium_nlp",
            "methods": ("stability_tp",),
            "ipopt_formulation": "thermodynamic_constrained_nlp",
            "route": "tangent_plane_distance",
        },
    },
    {
        "key": "electrolyte_stability",
        "public_routes": ("electrolyte_stability",),
        "payload": {
            "backend": "native_ipopt_equilibrium_nlp",
            "methods": ("electrolyte_stability_tp",),
            "ipopt_formulation": "thermodynamic_constrained_nlp",
            "route": "charge_constrained_tangent_plane_distance",
        },
    },
    {
        "key": "neutral_lle_flash",
        "public_routes": ("neutral_lle_flash",),
        "payload": {
            "backend": "native_ipopt_equilibrium_nlp",
            "methods": ("lle_flash", "lle_tp"),
            "ipopt_formulation": "thermodynamic_constrained_nlp",
        },
    },
    {
        "key": "neutral_bubble_dew",
        "public_routes": ("neutral_bubble_p", "neutral_bubble_t", "neutral_dew_p", "neutral_dew_t"),
        "payload": {
            "backend": "native_ipopt_equilibrium_nlp",
            "methods": ("bubble_p", "bubble_t", "dew_p", "dew_t"),
        },
    },
    {
        "key": "electrolyte_lle",
        "public_routes": ("electrolyte_lle",),
        "payload": {
            "backend": "native_ipopt_equilibrium_nlp",
            "methods": ("electrolyte_lle", "electrolyte_lle_tp"),
            "solver_backends": ("auto", "ipopt"),
            "explicit_ipopt_request": "routes_to_native_ipopt_when_compiled",
            "ipopt_formulation": "thermodynamic_constrained_nlp",
        },
    },
    {
        "key": "electrolyte_bubble_pressure",
        "public_routes": ("electrolyte_bubble_pressure",),
        "payload": {
            "backend": "native_ipopt_equilibrium_nlp",
            "scope": "fixed liquid composition with neutral vapor species; ions remain liquid-only",
        },
    },
)

EQUILIBRIUM_PROBLEM_OBJECT_CLASSES: Final[tuple[str, ...]] = (
    "TPFlash",
    "StabilityAnalysis",
    "BubblePoint",
    "DewPoint",
    "LLEProblem",
    "ElectrolyteLLEProblem",
    "ElectrolyteBubblePoint",
    "ReactiveSpeciationProblem",
    "ReactivePhaseEquilibriumProblem",
    "ReactiveElectrolyteBubbleProblem",
)

DERIVATIVE_COVERAGE_ROWS: Final[tuple[dict[str, object], ...]] = (
    {
        "row_family": "regression",
        "subsystem": "regression",
        "quantity": "pure_neutral_parameters",
        "derivative": "objective_jacobian",
        "backend": "cppad_implicit",
        "supported": True,
        "classification": "production_supported",
        "reason": "validated CppAD plus implicit density-sensitivity regression slice",
        "tests": ("tests/native/ceres/test_ceres_pure_regression.py",),
    },
    {
        "row_family": "regression",
        "subsystem": "regression",
        "quantity": "binary_kij",
        "derivative": "objective_jacobian",
        "backend": "cppad_implicit",
        "supported": True,
        "classification": "production_supported",
        "reason": "validated Ceres route with CppAD and implicit density/association sensitivities",
        "tests": (
            "tests/native/ceres/test_ceres_binary_regression.py",
            "tests/regression/literature/test_literature_binary_kij_regression.py",
        ),
    },
    {
        "row_family": "electrolyte_property",
        "subsystem": "electrolyte",
        "quantity": "ssmds_born_liquid",
        "derivative": "parameter_sensitivity",
        "backend": "analytic",
        "supported": True,
        "classification": "production_supported",
        "reason": "liquid electrolyte SSM+DS Born derivatives are analytic; vapor Born derivatives are not claimed",
        "tests": ("tests/api/runtime/test_runtime_ionic_methods.py",),
    },
)

EQUILIBRIUM_ROUTE_DERIVATIVE_EVIDENCE: Final[tuple[dict[str, object], ...]] = (
    {
        "row_family": "equilibrium",
        "subsystem": "native_ipopt",
        "quantity": "neutral_two_phase_routes",
        "derivative": "lagrangian_hessian",
        "backend": "cppad_phase_blocks",
        "supported": True,
        "classification": "production_supported",
        "reason": "neutral TP flash, LLE, stability, bubble, and dew routes require exact gradients/Jacobians and default to exact Hessians",
        "tests": (
            "tests/equilibrium/core/test_api.py",
            "tests/equilibrium/core/test_bubble_dew.py",
            "tests/native/equilibrium/test_route_builders.py",
        ),
    },
    {
        "row_family": "equilibrium",
        "subsystem": "native_ipopt",
        "quantity": "electrolyte_lle_and_stability",
        "derivative": "lagrangian_hessian",
        "backend": "cppad_explicit_density",
        "supported": True,
        "classification": "production_supported",
        "reason": "electrolyte LLE and charge-constrained stability routes report exact Hessian metadata and reject non-certified postsolves",
        "tests": (
            "tests/api/equilibrium/test_electrolyte_lle_problem_native_ipopt.py",
            "tests/equilibrium/electrolyte/test_electrolyte_lle_solver_contracts.py",
            "tests/native/equilibrium/test_route_builders.py",
        ),
    },
    {
        "row_family": "equilibrium",
        "subsystem": "native_ipopt",
        "quantity": "reactive_lle_and_reactive_electrolyte_lle",
        "derivative": "lagrangian_hessian",
        "backend": "cppad_explicit_density",
        "supported": True,
        "classification": "route_builder_supported_capability_pending",
        "reason": "native reactive LLE route builders expose exact Hessian paths, while production capability registration still waits for benchmark and acceptance evidence",
        "tests": (
            "tests/native/equilibrium/test_reactive_phase_equilibrium_residual_jacobian.py",
            "tests/native/equilibrium/test_route_builders.py",
            "tests/api/reactive/test_reactive_phase_equilibrium_problem_routes_native.py",
        ),
    },
)

REGRESSION_CAPABILITY_KEYS: Final[tuple[str, ...]] = (
    "pure_neutral",
    "pure_ion",
    "binary_pair",
    "mea_co2_h2o_electrolyte_benchmark",
    "reactive_electrolyte_residuals",
    "reactive_electrolyte_batch_context",
)

NATIVE_CONTRACT_TEST_TARGETS: Final[tuple[str, ...]] = (
    "tests/native/equilibrium/test_route_metadata_contracts.py",
    "tests/native/equilibrium/test_native_route_diagnostics_contract.py",
)

GENERIC_TEST_TARGETS: Final[tuple[str, ...]] = (
    "tests/api/package/test_package_main.py::test_python_m_epcsaft_reports_package_and_core_status",
    "tests/api/runtime/test_runtime_exports_and_metadata.py::test_package_exports_are_available",
    "tests/api/runtime/test_runtime_neutral_methods.py::test_neutral_scalar_methods_return_expected_values",
    "tests/api/runtime/test_runtime_ionic_methods.py::test_ionic_activity_and_solution_methods_return_expected_values",
    "tests/api/parameters/test_parameter_templates.py::test_runtime_options_accept_cppad_modes_and_preserve_explicit_overrides",
    "tests/api/regression/test_regression_api_native_backends.py::test_public_pure_neutral_regression_is_robust_to_distinct_initial_guesses",
    "tests/regression/core/test_hydrocarbon.py::test_methane_reference_parameters_keep_native_objective_pinned",
    "tests/equilibrium/core/test_vle.py::test_tp_flash_builds_one_native_route_request_before_ipopt_gate",
    "tests/equilibrium/core/test_lle.py::test_lle_flash_builds_one_native_route_request_before_ipopt_gate",
    "tests/equilibrium/core/test_stability.py::test_stability_uses_native_ipopt_route_after_validation",
    (
        "tests/equilibrium/electrolyte/test_electrolyte_lle_smokes.py::"
        "test_electrolyte_lle_builds_native_route_before_ipopt_gate"
    ),
    (
        "tests/workflows/validation/equilibrium_core/test_electrolyte_thermo_diagnostics.py::"
        "test_khudaida_fixture_loads_charge_neutral_explicit_ions"
    ),
    "tests/native/contracts/test_equilibrium_native_contracts.py::test_native_equilibrium_entrypoint_is_exposed",
    "tests/native/runtime/test_runtime_density_closure.py::test_pressure_based_and_density_based_states_match_for_neutral_system",
    "tests/native/contracts/test_equation_registry.py::test_equation_registry_outputs_are_synced",
    *NATIVE_CONTRACT_TEST_TARGETS,
    "tests/workflows/repo/test_project_structure.py",
    "tests/workflows/repo/test_run_pytest.py::test_list_slices_exits_without_running_pytest",
    "tests/workflows/repo/test_workflow_entrypoints.py::test_docs_make_confidence_suite_the_default_runtime_check",
)
CONFIDENCE_TEST_TARGETS: Final[tuple[str, ...]] = (
    *GENERIC_TEST_TARGETS,
    "tests/native/runtime/test_runtime_density_closure.py::test_pressure_based_and_density_based_states_match_for_ionic_system",
    "tests/native/runtime/test_runtime_contribution_contracts.py::test_native_residual_helmholtz_and_compressibility_contributions_match_neutral_contract",
    "tests/native/contracts/test_equilibrium_native_contracts.py::test_public_tp_flash_requires_native_ipopt_route",
)
EQUILIBRIUM_CONFIDENCE_TEST_TARGETS: Final[tuple[str, ...]] = (
    (
        "tests/workflows/validation/equilibrium_core/test_electrolyte_lle_confidence.py::"
        "test_khudaida_benchmark_fixture_loads_charge_neutral_cases"
    ),
    (
        "tests/workflows/validation/equilibrium_core/test_electrolyte_thermo_diagnostics.py::"
        "test_khudaida_package_tieline_fixed_phase_residual_is_internally_consistent"
    ),
)
EQUILIBRIUM_API_TEST_TARGETS: Final[tuple[str, ...]] = (
    "tests/equilibrium/core/test_vle.py::test_tp_flash_builds_one_native_route_request_before_ipopt_gate",
    "tests/equilibrium/core/test_lle.py::test_lle_flash_builds_one_native_route_request_before_ipopt_gate",
    "tests/equilibrium/core/test_stability.py::test_stability_uses_native_ipopt_route_after_validation",
    (
        "tests/equilibrium/electrolyte/test_electrolyte_lle_smokes.py::"
        "test_electrolyte_lle_builds_native_route_before_ipopt_gate"
    ),
    "tests/api/runtime/test_runtime_exports_and_metadata.py::test_runtime_build_info_and_capabilities_are_json_like",
    (
        "tests/api/reactive/test_reactive_speciation_results.py::"
        "test_solve_reactive_speciation_activity_coupled_state_requires_native_ipopt_route"
    ),
    "tests/api/reactive/test_reactive_speciation_options.py::test_reactive_speciation_options_expose_jacobian_backend_selector",
    "tests/api/reactive/test_reactive_speciation_options.py::test_reactive_speciation_requested_ipopt_routes_ideal_speciation_when_compiled",
    "tests/api/reactive/test_reactive_electrolyte_bubble_setup.py",
    "tests/api/reactive/test_reactive_electrolyte_bubble_results.py",
    (
        "tests/native/equilibrium/test_chemical_equilibrium_native_api.py::"
        "test_native_chemical_equilibrium_residual_evaluator_uses_analytic_jacobian_by_default"
    ),
    (
        "tests/native/equilibrium/test_chemical_equilibrium_native_errors.py::"
        "test_native_chemical_equilibrium_residual_evaluator_rejects_removed_backend"
    ),
)
RUNTIME_TEST_TARGETS: Final[tuple[str, ...]] = (
    "tests/api/runtime/test_runtime_exports_and_metadata.py",
    "tests/api/runtime/test_runtime_neutral_methods.py",
    "tests/api/runtime/test_runtime_ionic_methods.py",
    "tests/native/runtime/test_runtime_density_closure.py",
    "tests/native/runtime/test_runtime_contribution_contracts.py",
    "tests/native/runtime/test_runtime_cache_contracts.py",
)
API_TEST_TARGETS: Final[tuple[str, ...]] = (
    "tests/api/runtime/test_runtime_exports_and_metadata.py",
    "tests/api/runtime/test_runtime_neutral_methods.py",
    "tests/api/runtime/test_runtime_ionic_methods.py",
    "tests/api/parameters/test_parameter_templates.py",
    "tests/api/regression/test_regression_api_public_contracts.py",
    "tests/api/regression/test_regression_api_native_backends.py",
    "tests/api/regression/test_regression_api_results_and_errors.py",
)
NATIVE_TEST_TARGETS: Final[tuple[str, ...]] = (
    "tests/native/runtime/test_runtime_density_closure.py",
    "tests/native/runtime/test_runtime_contribution_contracts.py",
    "tests/native/runtime/test_runtime_cache_contracts.py",
)
PROFILE_TEST_TARGETS: Final[tuple[str, ...]] = ("tests/profile/test_runtime_profile.py",)
FULL_PROFILE_TEST_TARGETS: Final[tuple[str, ...]] = (
    "tests/profile/test_runtime_profile.py",
    "tests/profile/test_miac_profile.py",
    "tests/profile/test_regression_profile.py",
)

TEST_SLICES: Final[dict[str, dict[str, object]]] = {
    "generic": {"targets": GENERIC_TEST_TARGETS, "cheap_by_default": True},
    "all": {"targets": ("tests",), "cheap_by_default": False},
    "confidence": {"targets": CONFIDENCE_TEST_TARGETS, "cheap_by_default": False},
    "equilibrium-confidence": {"targets": EQUILIBRIUM_CONFIDENCE_TEST_TARGETS, "cheap_by_default": False},
    "equilibrium-api": {"targets": EQUILIBRIUM_API_TEST_TARGETS, "cheap_by_default": True},
    "runtime": {"targets": RUNTIME_TEST_TARGETS, "cheap_by_default": True},
    "api": {"targets": API_TEST_TARGETS, "cheap_by_default": True},
    "native": {"targets": NATIVE_TEST_TARGETS, "cheap_by_default": True},
    "native-contracts": {"targets": NATIVE_CONTRACT_TEST_TARGETS, "cheap_by_default": True},
    "profile": {"targets": PROFILE_TEST_TARGETS, "cheap_by_default": False},
    "profile-full": {"targets": FULL_PROFILE_TEST_TARGETS, "cheap_by_default": False},
}

VALIDATION_LANES: Final[dict[str, dict[str, object]]] = {
    "quick": {
        "commands": (("scripts/dev/doctor.py",), ("run_pytest.py", "-q")),
        "cheap_by_default": True,
        "evidence": "doctor plus generic pytest slice",
    },
    "confidence": {
        "commands": (("scripts/dev/doctor.py",), ("run_pytest.py", "--confidence", "-q")),
        "cheap_by_default": False,
        "evidence": "doctor plus native runtime confidence pytest slice",
    },
    "docs": {
        "commands": (("-m", "sphinx", "-b", "html", "docs", "build/docs-html"),),
        "cheap_by_default": False,
        "evidence": "Sphinx documentation build",
    },
    "full": {
        "commands": (("scripts/dev/doctor.py",), ("run_pytest.py", "--all", "-q")),
        "cheap_by_default": False,
        "evidence": "doctor plus exhaustive historical pytest suite",
    },
    "ceres-cppad": {
        "commands": (
            ("scripts/dev/build_epcsaft.py", "--profile", "full"),
            (
                "run_pytest.py",
                "tests/native/ceres/test_ceres_pure_regression.py",
                "tests/native/ceres/test_ceres_binary_regression.py",
                "-q",
            ),
        ),
        "cheap_by_default": False,
        "evidence": "Ceres and CppAD build profile plus focused native regression tests",
    },
}


def registry_targets(name: str) -> tuple[str, ...]:
    """Return the registered pytest targets for a named slice."""

    return tuple(str(target) for target in TEST_SLICES[name]["targets"])


def validation_lane_commands(name: str) -> tuple[tuple[str, ...], ...]:
    """Return executable commands for a named validation lane."""

    return tuple(tuple(str(part) for part in command) for command in VALIDATION_LANES[name]["commands"])


def registered_ipopt_public_routes() -> list[str]:
    """Return public Ipopt route labels registered by executable route evidence."""

    return [
        str(route)
        for evidence in IPOPT_EQUILIBRIUM_ROUTE_EVIDENCE
        for route in evidence["public_routes"]  # type: ignore[index]
    ]
