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
        pass
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


def _parse_cmake_cache(cache_path: Path | None) -> dict[str, str]:
    if cache_path is None or not cache_path.exists():
        return {}
    pattern = re.compile(r"^(?P<key>[^:]+):[^=]*=(?P<value>.*)$")
    entries: dict[str, str] = {}
    for line in cache_path.read_text(encoding="utf-8", errors="replace").splitlines():
        match = pattern.match(line)
        if match:
            entries[match.group("key")] = match.group("value")
    return entries


def _cmake_bool(entries: dict[str, str], key: str) -> bool | None:
    value = entries.get(key)
    if value is None:
        return None
    return value.upper() in {"ON", "YES", "TRUE", "1"}


def _native_dependency_profile(
    cache_path: Path | None, enabled_key: str, available_key: str, use_system_key: str
) -> dict[str, object]:
    entries = _parse_cmake_cache(cache_path)
    enabled = _cmake_bool(entries, enabled_key)
    found = _cmake_bool(entries, available_key)
    use_system = _cmake_bool(entries, use_system_key)
    if enabled is None:
        return {
            "enabled": False,
            "use_system": False,
            "found": False,
            "available": False,
            "status": "not_configured",
        }
    if not enabled:
        return {
            "enabled": False,
            "use_system": bool(use_system),
            "found": bool(found),
            "available": False,
            "status": "disabled",
        }
    if found is None:
        return {
            "enabled": True,
            "use_system": bool(use_system),
            "found": False,
            "available": False,
            "status": "not_detected",
        }
    if not found:
        return {
            "enabled": True,
            "use_system": bool(use_system),
            "found": False,
            "available": False,
            "status": "missing",
        }
    return {
        "enabled": True,
        "use_system": bool(use_system),
        "found": True,
        "available": True,
        "status": "available",
    }


def _native_dependency_status(source_root: Path | None) -> dict[str, dict[str, object]]:
    cache_path = None if source_root is None else source_root / "build" / "dev" / "CMakeCache.txt"
    return {
        "ceres": _native_dependency_profile(
            cache_path, "EPCSAFT_ENABLE_CERES", "EPCSAFT_CERES_AVAILABLE", "EPCSAFT_USE_SYSTEM_CERES"
        ),
        "cppad": _native_dependency_profile(
            cache_path, "EPCSAFT_ENABLE_CPPAD", "EPCSAFT_CPPAD_AVAILABLE", "EPCSAFT_USE_SYSTEM_CPPAD"
        ),
    }


def _finite_difference_debug_enabled() -> bool:
    return os.environ.get("EPCSAFT_ALLOW_FINITE_DIFFERENCE_DEBUG") == "1"


def _native_extension_path() -> Path | None:
    try:
        from . import _core
    except Exception:
        return None
    return Path(_core.__file__).resolve()


def _mtime_utc(path: Path) -> str:
    return datetime.fromtimestamp(path.stat().st_mtime, timezone.utc).isoformat()


@lru_cache(maxsize=1)
def runtime_build_info() -> dict[str, object]:
    """Return JSON-like package, source, and native-extension metadata."""

    from .ipopt_backend import cyipopt_backend_info

    direct_url = _direct_url_payload()
    source_root = _source_checkout_from_package() or _source_checkout_from_direct_url(direct_url)
    native_path = _native_extension_path()
    cyipopt = cyipopt_backend_info()
    native_dependencies = _native_dependency_status(source_root)
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
        "platform": platform.platform(),
        "machine": platform.machine(),
        "optional_dependencies": {
            "cyipopt": cyipopt,
        },
        "native_dependencies": native_dependencies,
    }


def capabilities() -> dict[str, object]:
    """Return structured availability flags for high-level package workflows."""

    build_info = runtime_build_info()
    cyipopt = dict(build_info["optional_dependencies"]["cyipopt"])  # type: ignore[index]
    native_dependencies = build_info["native_dependencies"]  # type: ignore[index]
    return {
        "native_extension": bool(build_info["native_extension_available"]),
        "native_dependencies": native_dependencies,
        "optimizers": {
            "ipopt": {
                **cyipopt,
                "solver_backend": "ipopt",
                "scope": "explicit opt-in cyipopt backend; native Newton/default backends remain defaults",
                "hessian_strategies": ["gauss_newton", "lbfgs"],
                "formulations": ["bound_constrained_residual_minimization"],
                "full_constrained_nlp_available": False,
                "default_auto_uses_ipopt": False,
                "exact_hessian_available": False,
                "hessian_includes_second_residual_derivatives": False,
            }
        },
        "equilibrium": {
            "neutral_tp_flash": {"available": True, "backend": "native"},
            "neutral_lle_flash": {"available": True, "backend": "native"},
            "neutral_stability": {"available": True, "backend": "native"},
            "neutral_bubble_dew": {
                "available": True,
                "backend": "native_state_fugacity_with_python_scalar_root",
                "methods": ["bubble_p", "bubble_t", "dew_p", "dew_t"],
                "status": "production",
            },
            "electrolyte_lle": {
                "available": True,
                "backend": "native",
                "solver_backends": ["auto", "newton", "ipopt"],
                "ipopt_available": bool(cyipopt["available"]),
                "default_auto_uses_ipopt": False,
                "ipopt_formulation": "bound_constrained_residual_minimization",
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
                "scope": "native chemical speciation followed by native fixed-liquid electrolyte bubble pressure",
            },
            "reactive_speciation": {
                "available": True,
                "backend": "native activity evaluations",
                "sweep_available": True,
                "continuation_state_available": True,
                "activity_output_modes": ["auto", "always", "never"],
                "jacobian_auto_policy": "analytic_ideal_else_backend_unavailable",
                "jacobian_auto_supported_standard_states": [
                    "ideal_mole_fraction",
                ],
                "finite_difference_requires_explicit_request": True,
                "finite_difference_debug_gate": "EPCSAFT_ALLOW_FINITE_DIFFERENCE_DEBUG=1",
                "finite_difference_debug_enabled": _finite_difference_debug_enabled(),
                "explicit_autodiff_raises_when_unavailable": True,
                "solver_backends": ["auto", "newton", "ipopt"],
                "ipopt_available": bool(cyipopt["available"]),
                "default_auto_uses_ipopt": False,
                "ipopt_formulation": "bound_constrained_residual_minimization",
                "full_constrained_nlp_available": False,
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
                "scope": "external parameter bundle structure and reaction/species compatibility checks",
            },
        },
        "regression": {
            "pure_neutral": {"available": True, "backend": "native"},
            "pure_ion": {"available": True, "backend": "native"},
            "binary_pair": {"available": True, "backend": "native"},
            "native_dependencies": native_dependencies,
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
                    "canonical_statuses_target": [
                        "converged",
                        "max_iterations",
                        "line_search_failed",
                        "singular_jacobian",
                        "all_rows_failed",
                        "nonfinite_objective",
                        "bounds_inconsistent",
                        "invalid_input",
                        "backend_unavailable",
                    ],
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
                    "status": "partial_native_ceres_thermodynamic_slice",
                    "optimizer": "native_ceres_thermodynamic_regression_for_supported_slice",
                    "compatibility_optimizer": "python_bounded_gauss_newton_compatibility",
                    "issue53_native_production_ready": False,
                    "supports_pressure_targets": True,
                    "supports_speciation_targets": True,
                    "supports_activity_targets": True,
                    "supports_fugacity_targets": True,
                    "supports_density_targets": True,
                    "supports_relative_permittivity_targets": True,
                    "supports_bounds": True,
                    "native_hot_loop": True,
                    "native_optimizer_boundary": True,
                    "thermodynamic_backend": "native",
                    "required_native_dependencies": ["ceres", "cppad"],
                    "public_default_backend": "native",
                    "compatibility_backend": "python_compat",
                    "production_blockers": [
                        "Ceres thermodynamic hot loop currently supports only reactive-speciation rows with reaction logK parameters, ideal standard states, and speciation targets",
                        "reactive electrolyte bubble, Born/SSM+DS d_born/f_solv, activity, pressure, and activity-coupled standard-state derivatives still report backend_unavailable",
                    ],
                    "missing_bubble_derivative_residuals": [
                        "log-pressure unknown for bubble pressure",
                        "vapor-composition unknowns and normalization",
                        "liquid fixed-composition fugacity/activity state",
                        "vapor submixture fugacity state",
                        "fugacity equality residuals log(y_i)+log(phi_i_vap)-log(x_i)-log(phi_i_liq)",
                        "pressure target residual through solved P",
                    ],
                    "missing_derivative_call_graph": [
                        "NativeThermoCeresCostFunction::Evaluate",
                        "evaluate_native_thermo_regression_rows",
                        "chemical_equilibrium_native/evaluate_chemical_equilibrium_residual_native",
                        "activity_coefficients",
                        "ePCSAFTStateNative::activity_coefficient_native",
                        "residual_chemical_potential_result_cpp",
                        "composition_contribution_result_cpp",
                        "ares_contributions_cpp",
                        "born_intermediate_state_cpp and dadx_born_cpp for Born/SSM+DS d_born and f_solv",
                        "solve_density_scoped when activity or concentration standards require pressure closure",
                    ],
                    "generic_binary_regression_path": "use fit_binary_pair(...) for k_ij, l_ij, and k_hb_ij against direct binary VLE composition data",
                    "python_role": "request serialization, diagnostics, optional compatibility optimizer",
                },
                "native_ceres_thermodynamic_regression": {
                    "available": True,
                    "backend": "native",
                    "optimizer_backend": "ceres_when_EPCSAFT_ENABLE_CERES",
                    "derivative_backend": "analytic_implicit",
                    "native_hot_loop": True,
                    "python_objective_used": False,
                    "finite_difference_used": False,
                    "supported_slice": {
                        "row_modes": ["reactive_speciation"],
                        "parameter_kinds": ["reaction_equilibrium_constant", "log_equilibrium_constant"],
                        "target_families": ["speciation"],
                        "reaction_standard_states": ["ideal_mole_fraction"],
                    },
                    "blocked_parameter_kinds": {
                        "born_radius": "requires activity/fugacity/EOS parameter sensitivities through Born/SSM+DS and density closure",
                        "born_diameter": "requires activity/fugacity/EOS parameter sensitivities through Born/SSM+DS and density closure",
                        "f_solv": "requires activity/fugacity/EOS parameter sensitivities through Born/SSM+DS solvation-shell terms and density closure",
                        "solvation_factor": "requires activity/fugacity/EOS parameter sensitivities through Born/SSM+DS solvation-shell terms and density closure",
                    },
                    "unsupported_status": "backend_unavailable",
                    "benchmark_command": "uv run python scripts/benchmark_native_ceres_thermo_regression.py --warmup 1 --repeat 3",
                },
                "classes": [
                    "ReactiveElectrolyteBatch",
                    "ReactiveElectrolyteRegressionContext",
                    "ReactiveRegressionObjective",
                ],
                "methods": ["evaluate_objective", "finite_difference_jacobian"],
                "benchmark_commands": [
                    "uv run python scripts/benchmark_reactive_regression.py --warmup 3 --repeat 10 --json build/benchmarks/reactive_regression_main.json",
                    "uv run python scripts/benchmark_reactive_regression.py --case reactive_regression_objective_tiny --warmup 3 --repeat 20 --json build/benchmarks/reactive_regression_objective_main.json",
                    "uv run python scripts/benchmark_reactive_regression.py --case reactive_regression_parameter_perturbation --warmup 3 --repeat 20 --json build/benchmarks/reactive_regression_perturbation_main.json",
                ],
            },
            "native_residual_record_regression": {
                "available": True,
                "backend": "native",
                "optimizer": "analytic_linear_native",
                "status": "production_contract_slice",
                "supports_bounds": True,
                "supports_fixed_shape_residuals": True,
                "supports_row_diagnostics": True,
                "supports_penalty_residuals": True,
                "production_finite_difference_allowed": False,
                "methods": [
                    "native_regression_contract_schema",
                    "evaluate_native_regression_residual_records",
                    "solve_native_regression_residual_records",
                ],
            },
        },
    }


__git_commit__ = str(runtime_build_info()["source_git_commit"])
