"""Runtime metadata and capability discovery for downstream applications."""

from __future__ import annotations

import json
import os
import platform
import re
import subprocess
import sys
from datetime import datetime, timezone
from functools import lru_cache
from importlib import metadata
from pathlib import Path
from urllib.parse import unquote, urlparse
from urllib.request import url2pathname


def _package_version() -> str:
    try:
        return metadata.version("epcsaft")
    except metadata.PackageNotFoundError:
        pyproject = Path(__file__).resolve().parents[2] / "pyproject.toml"
        if pyproject.exists():
            match = re.search(
                r'(?m)^version\s*=\s*"([^"]+)"',
                pyproject.read_text(encoding="utf-8", errors="replace"),
            )
            if match:
                return match.group(1)
        return "0+unknown"


__version__ = _package_version()


def _direct_url_payload() -> dict:
    try:
        text = metadata.distribution("epcsaft").read_text("direct_url.json")
    except metadata.PackageNotFoundError:
        return {}
    if not text:
        return {}
    try:
        payload = json.loads(text)
    except json.JSONDecodeError:
        return {}
    return payload if isinstance(payload, dict) else {}


def _path_from_file_url(url: str | None) -> Path | None:
    if not url:
        return None
    parsed = urlparse(str(url))
    if parsed.scheme != "file":
        return None
    raw_path = unquote(parsed.path)
    if parsed.netloc:
        raw_path = f"//{parsed.netloc}{raw_path}"
    return Path(url2pathname(raw_path))


def _source_checkout_from_package() -> Path | None:
    package_path = Path(__file__).resolve()
    for candidate in (package_path.parents[2], package_path.parents[1]):
        if (candidate / ".git").exists():
            return candidate
    return None


def _source_checkout_from_direct_url(payload: dict) -> Path | None:
    source = _path_from_file_url(payload.get("url"))
    if source is None:
        return None
    if source.is_file():
        source = source.parent
    for candidate in (source, *source.parents):
        if (candidate / ".git").exists():
            return candidate
    return source if source.exists() else None


def _git_commit(source_root: Path | None) -> str:
    if source_root is None or not source_root.exists():
        return "unknown"
    try:
        completed = subprocess.run(
            ["git", "-C", str(source_root), "rev-parse", "--short=12", "HEAD"],
            capture_output=True,
            text=True,
            check=False,
            timeout=2.0,
        )
    except (OSError, subprocess.SubprocessError):
        return "unknown"
    commit = completed.stdout.strip()
    return commit if completed.returncode == 0 and commit else "unknown"


def _native_extension_path() -> Path | None:
    try:
        from . import _core
    except (ImportError, OSError):
        return None
    return Path(_core.__file__).resolve()


def _platform_label() -> str:
    if os.name == "nt":
        return sys.platform
    return platform.platform()


def _machine_label() -> str:
    if os.name == "nt":
        return os.environ.get("PROCESSOR_ARCHITECTURE", "unknown")
    return platform.machine()


def _native_cppad_backend_info() -> dict[str, object]:
    try:
        from . import _core
    except (ImportError, OSError):
        return {
            "backend": "cppad",
            "status": "not_configured",
            "compiled": False,
            "available": False,
        }
    try:
        smoke = _core._native_cppad_smoke()
    except AttributeError:
        return {
            "backend": "cppad",
            "status": "not_configured",
            "compiled": False,
            "available": False,
        }
    status = str(smoke.get("status", "not_configured"))
    compiled = bool(smoke.get("cppad_compiled", False))
    return {
        "backend": "cppad",
        "status": status,
        "compiled": compiled,
        "available": status == "enabled_available" and compiled,
    }


def _native_ceres_backend_info() -> dict[str, object]:
    try:
        from . import _core
    except (ImportError, OSError):
        return {
            "backend": "ceres",
            "status": "not_configured",
            "compiled": False,
            "available": False,
        }
    try:
        smoke = _core._native_ceres_smoke()
    except AttributeError:
        return {
            "backend": "ceres",
            "status": "not_configured",
            "compiled": False,
            "available": False,
        }
    status = str(smoke.get("status", "not_configured"))
    compiled = bool(smoke.get("compiled", False))
    return {
        "backend": "ceres",
        "status": status,
        "compiled": compiled,
        "available": status == "enabled_available" and compiled,
    }


def _native_ipopt_backend_info() -> dict[str, object]:
    try:
        from . import _core
    except (ImportError, OSError):
        return {
            "backend": "ipopt",
            "status": "not_configured",
            "compiled": False,
            "available": False,
            "adapter_available": False,
            "adapter_kind": "native_tnlp_adapter",
            "adapter_source_available": False,
            "hessian_strategy": "limited_memory",
            "requires_exact_gradient": True,
            "requires_exact_jacobian": True,
            "requires_exact_hessian": False,
        }
    try:
        smoke = _core._native_ipopt_smoke()
    except AttributeError:
        return {
            "backend": "ipopt",
            "status": "not_configured",
            "compiled": False,
            "available": False,
            "adapter_available": False,
            "adapter_kind": "native_tnlp_adapter",
            "adapter_source_available": False,
            "hessian_strategy": "limited_memory",
            "requires_exact_gradient": True,
            "requires_exact_jacobian": True,
            "requires_exact_hessian": False,
        }
    status = str(smoke.get("status", "not_configured"))
    compiled = bool(smoke.get("compiled", False))
    return {
        "backend": "ipopt",
        "status": status,
        "compiled": compiled,
        "available": status == "enabled_available" and compiled,
        "adapter_available": bool(smoke.get("adapter_available", False)),
        "adapter_kind": str(smoke.get("adapter_kind", "native_tnlp_adapter")),
        "adapter_source_available": bool(smoke.get("adapter_source_available", False)),
        "hessian_strategy": str(smoke.get("hessian_strategy", "limited_memory")),
        "requires_exact_gradient": bool(smoke.get("requires_exact_gradient", True)),
        "requires_exact_jacobian": bool(smoke.get("requires_exact_jacobian", True)),
        "requires_exact_hessian": bool(smoke.get("requires_exact_hessian", False)),
    }


def _mtime_utc(path: Path) -> str:
    return datetime.fromtimestamp(path.stat().st_mtime, timezone.utc).isoformat()


@lru_cache(maxsize=1)
def runtime_build_info() -> dict[str, object]:
    """Return JSON-like package, source, and native-extension metadata."""

    direct_url = _direct_url_payload()
    source_root = _source_checkout_from_package() or _source_checkout_from_direct_url(direct_url)
    native_path = _native_extension_path()
    ceres = _native_ceres_backend_info()
    direct_url_info = direct_url.get("dir_info") if isinstance(direct_url.get("dir_info"), dict) else {}
    return {
        "package_version": __version__,
        "source_git_commit": _git_commit(source_root),
        "source_root": None if source_root is None else str(source_root),
        "direct_url": direct_url.get("url"),
        "editable": bool(direct_url_info.get("editable", False)),
        "package_file": str(Path(__file__).resolve()),
        "native_extension": None if native_path is None else str(native_path),
        "native_extension_available": native_path is not None,
        "native_extension_mtime_utc": None if native_path is None else _mtime_utc(native_path),
        "python": sys.version.split()[0],
        "platform": _platform_label(),
        "machine": _machine_label(),
        "optional_dependencies": {
            "ceres": ceres,
            "cppad": _native_cppad_backend_info(),
            "ipopt": _native_ipopt_backend_info(),
        },
    }


def _dependency_capability(
    dependency: dict[str, object],
    *,
    production: bool = False,
    reason: str | None = None,
    **extra: object,
) -> dict[str, object]:
    available = bool(dependency.get("available", False))
    payload: dict[str, object] = {
        **dependency,
        "production": bool(production and available),
    }
    if reason is not None:
        payload["reason"] = reason
    elif not available:
        payload["reason"] = "dependency_not_compiled"
    elif not payload["production"]:
        payload["reason"] = "not_validated_for_production"
    payload.update(extra)
    return payload


def _derivative_coverage_capabilities(cppad: dict[str, object], ceres: dict[str, object]) -> dict[str, object]:
    cppad_available = bool(cppad.get("available", False))
    ceres_available = bool(ceres.get("available", False))
    coverage_rows = [
        {
            "row_family": "regression",
            "subsystem": "regression",
            "quantity": "pure_neutral_parameters",
            "derivative": "objective_jacobian",
            "backend": "cppad_implicit",
            "supported": True,
            "not_applicable": False,
            "classification": "production_supported",
            "reason": "validated CppAD plus implicit density-sensitivity regression slice",
            "tests": ["tests/native/ceres/test_ceres_pure_regression.py"],
        },
        {
            "row_family": "regression",
            "subsystem": "regression",
            "quantity": "binary_kij",
            "derivative": "objective_jacobian",
            "backend": "cppad" if cppad_available and ceres_available else "not_available",
            "supported": False,
            "not_applicable": False,
            "classification": "blocker",
            "reason": (
                "requires Ceres and CppAD compiled"
                if not (cppad_available and ceres_available)
                else "explicit Ceres/CppAD residual Jacobians lack production validation"
            ),
            "tests": [
                "tests/native/contracts/test_ceres_cppad_build_contract.py",
                "tests/regression/literature/test_literature_binary_kij_regression.py",
            ],
        },
        {
            "row_family": "regression",
            "subsystem": "regression",
            "quantity": "binary_khb_ij",
            "derivative": "parameter_sensitivity",
            "backend": "not_available",
            "supported": False,
            "not_applicable": False,
            "classification": "blocker",
            "blocker_detail": "blocker_requires_implicit_association_sensitivity",
            "reason": "k_hb_ij is a real binary association parameter family, but full property/regression derivatives require implicit association site-fraction sensitivities",
            "future_owner": "Task C",
            "parameters": ["k_hb_ij"],
            "tests": [
                "tests/native/contracts/test_ceres_cppad_build_contract.py",
                "tests/native/contracts/test_derivative_coverage_matrix.py",
            ],
        },
        {
            "row_family": "solved_state",
            "subsystem": "equilibrium",
            "quantity": "association_site_fractions",
            "derivative": "direct_cppad_recording",
            "backend": "not_available",
            "supported": False,
            "not_applicable": True,
            "classification": "out_of_scope",
            "reason": "active association uses solved site fractions; production derivative is implicit",
            "tests": ["tests/native/contracts/test_association_implicit_derivative_contract.py"],
        },
        {
            "row_family": "solved_state",
            "subsystem": "equilibrium",
            "quantity": "bubble_pressure",
            "derivative": "root_implicit_sensitivity",
            "backend": "not_available",
            "supported": False,
            "not_applicable": False,
            "classification": "blocker",
            "reason": "implicit bubble-pressure sensitivity lacks production validation",
            "tests": ["tests/native/cppad/test_cppad_bubble_derivatives.py"],
        },
        {
            "row_family": "electrolyte_property",
            "subsystem": "electrolyte",
            "quantity": "ssmds_born_liquid",
            "derivative": "parameter_sensitivity",
            "backend": "analytic",
            "supported": True,
            "not_applicable": False,
            "classification": "production_supported",
            "reason": "liquid electrolyte SSM+DS Born derivatives are analytic; vapor Born derivatives are not claimed",
            "tests": ["tests/api/runtime/test_runtime_ionic_methods.py"],
        },
    ]
    return {
        "derivative_coverage_matrix_available": True,
        "minimum_columns": [
            "row_family",
            "subsystem",
            "quantity",
            "derivative",
            "backend",
            "supported",
            "not_applicable",
            "classification",
            "reason",
            "tests",
        ],
        "rows": coverage_rows,
        "association_direct_cppad_recording": {
            "available": False,
            "production": False,
            "reason": "active association uses solved site fractions; production derivative is implicit",
        },
        "association_implicit_sensitivities": {
            "available": True,
            "production": True,
            "scope": "validated association solved-state reporting and implicit-sensitivity diagnostics",
        },
        "density_root_implicit_sensitivities": {
            "available": False,
            "production": False,
            "reason": "not_available",
        },
        "speciation_implicit_sensitivities": {
            "available": True,
            "production": True,
            "scope": "standard-state slices with analytic residual Jacobians",
        },
        "bubble_pressure_implicit_sensitivities": {
            "available": False,
            "production": False,
            "reason": "not_available",
        },
        "born_ssmds_liquid_derivatives": {
            "available": True,
            "production": True,
            "phase_scope": "liquid_electrolyte_only",
            "vapor_support": False,
        },
        "regression_ceres_explicit_cppad_jacobians": {
            "available": bool(cppad_available and ceres_available),
            "production": False,
            "reason": (
                "not_validated_for_production" if cppad_available and ceres_available else "dependency_not_compiled"
            ),
            "requires": ["ceres", "cppad"],
        },
        "regression_ceres_implicit_jacobians": {
            "available": False,
            "production": False,
            "reason": "not_available",
        },
    }


def capabilities() -> dict[str, object]:
    """Return structured availability flags for high-level package workflows."""

    ceres = dict(runtime_build_info()["optional_dependencies"]["ceres"])  # type: ignore[index]
    cppad = dict(runtime_build_info()["optional_dependencies"]["cppad"])  # type: ignore[index]
    ipopt = dict(runtime_build_info()["optional_dependencies"]["ipopt"])  # type: ignore[index]
    ipopt_public_routes = ["reactive_speciation:ideal_mole_fraction"]
    ipopt_route_available = bool(ipopt.get("available", False))
    cppad_capability = _dependency_capability(
        cppad,
        production=False,
        scope="package-wide AD substrate",
        production_eos_coverage=False,
    )
    ceres_capability = _dependency_capability(
        ceres,
        production=bool(ceres.get("available", False)),
        scope="native optimizer backend for explicitly supported regression paths",
        native_hot_loop=bool(ceres.get("available", False)),
        production_routes=["regression:pure_neutral", "regression:pure_ion", "regression:binary_pair"],
    )
    derivative_coverage = _derivative_coverage_capabilities(cppad, ceres)
    return {
        "native_extension": bool(runtime_build_info()["native_extension_available"]),
        "derivatives": {
            "numerical_derivative": {
                "available": False,
                "production": False,
                "reason": "numerical_derivative_derivatives_forbidden",
            },
            "cppad": cppad_capability,
            "coverage_matrix": derivative_coverage,
            "ssmds_born_derivatives": {
                "available": True,
                "production": True,
                "backend": "analytic",
                "phase_scope": "liquid_electrolyte_only",
                "parameters": ["d_born", "f_solv"],
                "vapor_support": False,
            },
            "property_derivative_result_apis": {
                "available": True,
                "result_shape": [
                    "supported",
                    "backend",
                    "derivative_backend",
                    "message",
                    "value",
                    "jacobian",
                    "outputs",
                    "variables",
                    "shape",
                ],
                "backend_labels": [
                    "analytic",
                    "cppad",
                    "analytic_implicit",
                    "cppad_implicit",
                    "not_available",
                ],
                "numerical_derivative_backend_available": False,
                "parameter_families": {
                    "production_supported": [
                        "m",
                        "sigma",
                        "epsilon",
                        "k_ij",
                        "l_ij",
                        "d_born",
                        "f_solv",
                        "relative_permittivity",
                    ],
                    "blocker_requires_implicit_association_sensitivity": {
                        "parameters": ["k_hb_ij"],
                        "future_owner": "Task C",
                    },
                },
                "state_methods": [
                    "pressure_density_derivative_result",
                    "pressure_composition_derivative_result",
                    "pressure_parameter_derivative_result",
                    "density_pressure_derivative_result",
                    "ares_composition_derivative_result",
                    "chemical_potential_composition_derivative_result",
                    "chemical_potential_parameter_derivative_result",
                    "ln_fugacity_composition_derivative_result",
                    "ln_fugacity_parameter_derivative_result",
                    "activity_composition_derivative_result",
                    "activity_parameter_derivative_result",
                    "relative_permittivity_composition_derivative_result",
                    "relative_permittivity_parameter_derivative_result",
                    "derivative_coverage_matrix",
                ],
            },
        },
        "optimizers": {
            "ceres": ceres_capability,
            "ipopt": {
                **ipopt,
                "solver_backend": "ipopt",
                "production": ipopt_route_available,
                "scope": "native Ipopt dependency for production equilibrium NLP routes",
                "hessian_strategies": ["limited_memory"],
                "formulations": ["thermodynamic_constrained_nlp"],
                "adapter_available": bool(ipopt.get("adapter_available", False)),
                "adapter_source_available": bool(ipopt.get("adapter_source_available", False)),
                "adapter_kind": ipopt.get("adapter_kind", "native_tnlp_adapter"),
                "public_routes": ipopt_public_routes,
                "full_constrained_nlp_available": ipopt_route_available,
                "default_auto_uses_ipopt": False,
                "exact_hessian_available": False,
                "hessian_includes_second_residual_derivatives": False,
            },
        },
        "equilibrium": {
            "derivative_policy": {
                "numerical_derivative_backend_available": False,
                "accepted_derivative_backends": [
                    "auto",
                    "analytic",
                    "cppad",
                    "analytic_implicit",
                    "cppad_implicit",
                ],
                "auto_policy": "analytic_or_cppad_or_implicit_else_raise",
                "unsupported_derivative_behavior": "raise",
                "diagnostic_fields": [
                    "thermodynamic_backend",
                    "solver_backend",
                    "derivative_backend",
                    "derivative_status",
                    "residual_norm",
                    "best_state_available",
                    "best_state",
                    "row_failure_count",
                    "solved_internal_states",
                    "derivative_backend_by_block",
                    "implicit_sensitivity_blocks",
                    "residual_norm_by_block",
                    "association_solver_status",
                ],
            },
            "neutral_tp_flash": {"available": True, "backend": "native"},
            "neutral_lle_flash": {
                "available": True,
                "backend": "native",
                "solver": "ceres_trust_region_residual_solve",
                "solver_dependency": "ceres",
                "ceres_available": bool(ceres["available"]),
                "derivative_backend": "not_applicable",
                "stability_precheck": "neutral_tpd",
                "anti_trivial_solution_strategy": "phase_fraction_and_phase_distance_gate",
                "numerical_derivative_backend_available": False,
            },
            "neutral_stability": {"available": True, "backend": "native"},
            "neutral_bubble_dew": {
                "available": False,
                "backend": "native_ipopt_equilibrium_nlp_required",
                "methods": ["bubble_p", "bubble_t", "dew_p", "dew_t"],
                "status": "route_pending",
            },
            "electrolyte_lle": {
                "available": True,
                "backend": "native",
                "solver_backends": ["auto"],
                "ipopt_available": bool(ipopt["available"]),
                "explicit_ipopt_request": "raises_until_native_adapter_is_implemented",
                "default_auto_uses_ipopt": False,
                "ipopt_formulation": "thermodynamic_constrained_nlp",
                "full_constrained_nlp_available": False,
            },
            "electrolyte_bubble_pressure": {
                "available": True,
                "backend": "native",
                "scope": "fixed liquid composition with neutral vapor species; ions remain liquid-only",
            },
            "reactive_electrolyte_bubble": {
                "available": True,
                "backend": "native",
                "scope": "native chemical speciation with fixed-liquid native bubble-pressure handoff and explicit partial-pressure diagnostics",
            },
            "reactive_phase_equilibrium": {
                "available": True,
                "backend": "native_coupled_reactive_phase_equilibrium",
                "methods": ["reactive_lle", "reactive_electrolyte_lle"],
                "problem_class": "ReactivePhaseEquilibriumProblem",
                "solver_dependency": "ceres",
                "ceres_available": bool(ceres["available"]),
                "reaction_semantics": "reaction residuals support both per-phase same-stoichiometry reactions and phase-tagged cross-phase quotients",
                "supported_reaction_scopes": ["same_phase_activity_reaction", "phase_tagged_cross_phase_quotient"],
                "unsupported_reaction_scopes": [],
                "cross_phase_reaction_quotients": {
                    "available": True,
                    "status": "supported",
                    "api": "ReactionDefinition.phase_stoichiometry",
                    "phase_terms": ["phase1/liq1/aq", "phase2/liq2/org"],
                },
            },
            "reactive_speciation": {
                "available": True,
                "backend": "native chemical-equilibrium residual with ePC-SAFT activity/fugacity terms",
                "sweep_available": True,
                "continuation_state_available": True,
                "activity_output_modes": ["auto", "always", "never"],
                "jacobian_auto_policy": "native_analytic_log_amount_residual_jacobian_with_implicit_sensitivity",
                "jacobian_auto_supported_standard_states": [
                    "ideal_mole_fraction",
                    "mole_fraction_activity",
                    "thermodynamic_activity",
                    "concentration",
                    "apparent",
                ],
                "derivative_gap_status": "implicit_sensitivity_available_for_reaction_constant_response",
                "explicit_cppad_request_raises_until_implemented": True,
                "solver_backends": ["auto", "ipopt"],
                "ipopt_available": bool(ipopt["available"]),
                "explicit_ipopt_request": "ideal_mole_fraction_routes_to_native_ipopt_when_compiled",
                "ipopt_routes": ipopt_public_routes,
                "default_auto_uses_ipopt": False,
                "ipopt_formulation": "thermodynamic_constrained_nlp",
                "full_constrained_nlp_available": ipopt_route_available,
            },
            "repeated_state_properties": {
                "available": True,
                "helpers": ["evaluate_fugacity_coefficients", "evaluate_fugacity_coefficients_batch"],
                "density_seed_aliases": ["rho_guess", "rho_seed"],
            },
            "problem_objects": {
                "available": True,
                "backend": "typed_python_dispatch_to_existing_native_methods",
                "classes": [
                    "TPFlash",
                    "StabilityAnalysis",
                    "BubblePoint",
                    "DewPoint",
                    "LLEProblem",
                    "ElectrolyteLLEProblem",
                    "ElectrolyteBubblePoint",
                    "ReactiveSpeciationProblem",
                    "ReactiveElectrolyteBubbleProblem",
                ],
                "entrypoint": "mixture.solve_equilibrium(problem)",
            },
            "contribution_maps": {
                "available": True,
                "backend": "native_term_payloads",
                "families": ["hard_chain", "dispersion", "association", "ionic", "born"],
                "inactive_terms_explicit": True,
                "activity_coefficient_term_decomposition_available": False,
            },
            "dataset_validation": {
                "available": True,
                "helper": "validate_dataset_bundle",
                "scope": "external parameter bundle structure and reaction/species consistency checks",
            },
        },
        "regression": {
            "pure_neutral": {"available": True, "backend": "native_ceres"},
            "pure_ion": {"available": True, "backend": "native_ceres"},
            "binary_pair": {"available": True, "backend": "native_ceres"},
            "mea_co2_h2o_electrolyte_benchmark": {
                "available": True,
                "backend": "native",
                "scope": "fixed-composition benchmark, not reactive bubble-pressure fitting",
            },
            "reactive_electrolyte_residuals": {
                "available": True,
                "backend": "python_orchestrated_native_solvers",
                "scope": "fixed-shape residual evaluator for downstream-owned coupled reactive electrolyte regression loops",
            },
            "reactive_electrolyte_batch_context": {
                "available": True,
                "backend": "python_batched_native_solvers",
                "fit_status_contract": {
                    "available": True,
                    "statuses": ["converged", "max_iterations", "line_search_failed", "failed_rows"],
                    "fields": [
                        "status",
                        "termination_reason",
                        "objective_initial",
                        "objective_final",
                        "gradient_norm",
                        "step_norm",
                    ],
                    "public_placeholder_statuses": [],
                },
                "bounded_mixed_pressure_speciation_regression": {
                    "available": True,
                    "status": "production",
                    "optimizer": "bounded_gauss_newton_least_squares",
                    "supports_pressure_targets": True,
                    "supports_speciation_targets": True,
                    "supports_activity_targets": True,
                    "supports_fugacity_targets": True,
                    "supports_density_targets": True,
                    "supports_relative_permittivity_targets": True,
                    "supports_bounds": True,
                    "native_hot_loop": False,
                    "ceres": {
                        "available": bool(ceres["available"]),
                        "production": False,
                        "reason": ("reactive batch regression is python-orchestrated bounded Gauss-Newton, not Ceres"),
                    },
                    "thermodynamic_backend": "native",
                    "python_role": "row orchestration, bounded step control, diagnostics",
                },
                "classes": [
                    "ReactiveElectrolyteBatch",
                    "ReactiveElectrolyteRegressionContext",
                    "ReactiveRegressionObjective",
                ],
                "methods": ["evaluate_objective"],
                "benchmark_commands": [
                    "uv run python scripts/benchmarks/benchmark_reactive_regression.py --warmup 3 --repeat 10 --json build/benchmarks/reactive_regression_main.json",
                    "uv run python scripts/benchmarks/benchmark_reactive_regression.py --case reactive_regression_objective_tiny --warmup 3 --repeat 20 --json build/benchmarks/reactive_regression_objective_main.json",
                    "uv run python scripts/benchmarks/benchmark_reactive_regression.py --case reactive_regression_parameter_perturbation --warmup 3 --repeat 20 --json build/benchmarks/reactive_regression_perturbation_main.json",
                ],
            },
        },
    }


__git_commit__ = str(runtime_build_info()["source_git_commit"])
