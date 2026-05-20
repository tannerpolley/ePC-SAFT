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

from .capability_evidence import (
    DERIVATIVE_COVERAGE_ROWS,
    EQUILIBRIUM_PROBLEM_OBJECT_CLASSES,
    IPOPT_EQUILIBRIUM_ROUTE_EVIDENCE,
    REGRESSION_CAPABILITY_KEYS,
    TEST_SLICES,
    VALIDATION_LANES,
    registered_ipopt_public_routes,
)


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
            "status": "required_native_extension_missing",
            "required": True,
            "compiled": False,
            "available": False,
        }
    try:
        smoke = _core._native_cppad_smoke()
    except AttributeError:
        return {
            "backend": "cppad",
            "status": "required_cppad_smoke_missing",
            "required": True,
            "compiled": False,
            "available": False,
        }
    status = str(smoke.get("status", "required_cppad_status_missing"))
    compiled = bool(smoke.get("cppad_compiled", False))
    return {
        "backend": "cppad",
        "status": status,
        "required": True,
        "compiled": compiled,
        "available": status == "enabled_available" and compiled,
    }


def _native_ceres_backend_info() -> dict[str, object]:
    try:
        from . import _core
    except (ImportError, OSError):
        return {
            "backend": "ceres",
            "status": "required_native_extension_missing",
            "required": True,
            "compiled": False,
            "available": False,
        }
    try:
        smoke = _core._native_ceres_smoke()
    except AttributeError:
        return {
            "backend": "ceres",
            "status": "required_ceres_smoke_missing",
            "required": True,
            "compiled": False,
            "available": False,
        }
    status = str(smoke.get("status", "required_ceres_status_missing"))
    compiled = bool(smoke.get("compiled", False))
    return {
        "backend": "ceres",
        "status": status,
        "required": True,
        "compiled": compiled,
        "available": status == "enabled_available" and compiled,
    }


def _native_ipopt_backend_info() -> dict[str, object]:
    try:
        from . import _core
    except (ImportError, OSError):
        return {
            "backend": "ipopt",
            "status": "native_extension_missing",
            "required": False,
            "compiled": False,
            "available": False,
            "adapter_available": False,
            "adapter_kind": "native_tnlp_adapter",
            "adapter_source_available": False,
            "requires_exact_gradient": True,
            "requires_exact_jacobian": True,
        }
    try:
        smoke = _core._native_ipopt_smoke()
    except AttributeError:
        return {
            "backend": "ipopt",
            "status": "ipopt_probe_missing",
            "required": False,
            "compiled": False,
            "available": False,
            "adapter_available": False,
            "adapter_kind": "native_tnlp_adapter",
            "adapter_source_available": False,
            "requires_exact_gradient": True,
            "requires_exact_jacobian": True,
        }
    status = str(smoke.get("status", "ipopt_probe_missing"))
    compiled = bool(smoke.get("compiled", False))
    return {
        "backend": "ipopt",
        "status": status,
        "required": False,
        "compiled": compiled,
        "available": status == "enabled_available" and compiled,
        "adapter_available": bool(smoke.get("adapter_available", False)),
        "adapter_kind": str(smoke.get("adapter_kind", "native_tnlp_adapter")),
        "adapter_source_available": bool(smoke.get("adapter_source_available", False)),
        "requires_exact_gradient": bool(smoke.get("requires_exact_gradient", True)),
        "requires_exact_jacobian": bool(smoke.get("requires_exact_jacobian", True)),
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
        "native_dependencies": {
            "ceres": ceres,
            "cppad": _native_cppad_backend_info(),
            "ipopt": _native_ipopt_backend_info(),
        },
    }


def _capability_value(value: object) -> object:
    if isinstance(value, tuple):
        return [_capability_value(item) for item in value]
    if isinstance(value, dict):
        return {str(key): _capability_value(item) for key, item in value.items()}
    return value


def _registered_ipopt_public_routes() -> list[str]:
    return registered_ipopt_public_routes()


def _registered_equilibrium_route_capabilities(
    *,
    ipopt_route_available: bool,
    ipopt: dict[str, object],
) -> dict[str, dict[str, object]]:
    routes: dict[str, dict[str, object]] = {}
    for evidence in IPOPT_EQUILIBRIUM_ROUTE_EVIDENCE:
        payload = _capability_value(evidence["payload"])
        assert isinstance(payload, dict)
        route_payload: dict[str, object] = {
            "available": ipopt_route_available,
            **payload,
        }
        if evidence["key"] == "reactive_speciation":
            route_payload["sweep_available"] = ipopt_route_available
            route_payload["ipopt_available"] = bool(ipopt["available"])
            route_payload["ipopt_routes"] = list(evidence["public_routes"])
            route_payload["ideal_speciation_nlp_available"] = ipopt_route_available
            route_payload["nonideal_speciation_nlp_available"] = ipopt_route_available
            for key in (
                "sweep_available_from_ipopt",
                "ideal_speciation_nlp_available_from_ipopt",
                "nonideal_speciation_nlp_available_from_ipopt",
            ):
                route_payload.pop(key, None)
        elif evidence["key"] == "electrolyte_lle":
            route_payload["ipopt_available"] = bool(ipopt["available"])
        routes[str(evidence["key"])] = route_payload
    return routes


def _capability_evidence_summary(derivative_coverage: dict[str, object]) -> dict[str, object]:
    derivative_rows = derivative_coverage.get("rows", [])
    return {
        "source": "registered_capability_evidence",
        "equilibrium_keys": [str(item["key"]) for item in IPOPT_EQUILIBRIUM_ROUTE_EVIDENCE],
        "ipopt_public_routes": _registered_ipopt_public_routes(),
        "problem_object_classes": list(EQUILIBRIUM_PROBLEM_OBJECT_CLASSES),
        "regression_keys": list(REGRESSION_CAPABILITY_KEYS),
        "derivative_row_count": len(derivative_rows) if isinstance(derivative_rows, list) else 0,
        "pytest_slices": list(TEST_SLICES),
        "validation_lanes": list(VALIDATION_LANES),
        "cheap_validation_lanes": [
            name for name, lane in VALIDATION_LANES.items() if bool(lane["cheap_by_default"])
        ],
    }


def _derivative_coverage_capabilities(cppad: dict[str, object], ceres: dict[str, object]) -> dict[str, object]:
    cppad_available = bool(cppad.get("available", False))
    ceres_available = bool(ceres.get("available", False))
    coverage_rows = [_capability_value(row) for row in DERIVATIVE_COVERAGE_ROWS]
    return {
        "derivative_coverage_matrix_available": True,
        "implemented_routes_only": True,
        "minimum_columns": [
            "row_family",
            "subsystem",
            "quantity",
            "derivative",
            "backend",
            "supported",
            "classification",
            "reason",
            "tests",
        ],
        "rows": coverage_rows,
        "association_implicit_sensitivities": {
            "available": True,
            "production": True,
            "scope": "validated association solved-state reporting and implicit-sensitivity diagnostics",
        },
        "density_root_implicit_sensitivities": {
            "available": True,
            "production": True,
            "scope": "validated CppAD-implicit density sensitivities used by native regression Jacobians",
        },
        "speciation_implicit_sensitivities": {
            "available": True,
            "production": True,
            "scope": "ideal analytic/CppAD and nonideal CppAD-implicit solved-state sensitivities",
        },
        "born_ssmds_liquid_derivatives": {
            "available": True,
            "production": True,
            "phase_scope": "liquid_electrolyte_only",
            "vapor_support": False,
        },
        "regression_ceres_jacobians": {
            "available": bool(cppad_available and ceres_available),
            "production": bool(cppad_available and ceres_available),
            "routes": ["pure_neutral_parameters", "binary_kij"],
            "backends": ["cppad_implicit"],
            "requires": ["ceres", "cppad"],
        },
    }


def capabilities() -> dict[str, object]:
    """Return structured availability flags for high-level package workflows."""

    build_info = runtime_build_info()
    native_dependencies = build_info["native_dependencies"]  # type: ignore[index]
    ceres = dict(native_dependencies["ceres"])  # type: ignore[index]
    cppad = dict(native_dependencies["cppad"])  # type: ignore[index]
    ipopt = dict(native_dependencies["ipopt"])  # type: ignore[index]
    ipopt_public_routes = _registered_ipopt_public_routes()
    ipopt_route_available = bool(ipopt.get("available", False))
    ceres_available = bool(ceres.get("available", False))
    cppad_capability = {
        **cppad,
        "scope": "package-wide AD substrate; production derivative routes are listed in coverage_matrix",
    }
    ceres_capability: dict[str, object] = {
        **ceres,
        "production": ceres_available,
        "scope": "native optimizer backend for explicitly supported regression paths",
        "native_hot_loop": ceres_available,
        "production_routes": ["regression:pure_neutral", "regression:pure_ion", "regression:binary_pair"],
    }
    if not ceres_available:
        ceres_capability["reason"] = "required_native_dependency_missing"
    derivative_coverage = _derivative_coverage_capabilities(cppad, ceres)
    equilibrium_route_capabilities = _registered_equilibrium_route_capabilities(
        ipopt_route_available=ipopt_route_available,
        ipopt=ipopt,
    )
    return {
        "native_extension": bool(build_info["native_extension_available"]),
        "capability_evidence": _capability_evidence_summary(derivative_coverage),
        "derivatives": {
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
                    "cppad_explicit_density",
                ],
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
                "formulations": ["thermodynamic_constrained_nlp"],
                "adapter_available": bool(ipopt.get("adapter_available", False)),
                "adapter_source_available": bool(ipopt.get("adapter_source_available", False)),
                "adapter_kind": ipopt.get("adapter_kind", "native_tnlp_adapter"),
                "public_routes": ipopt_public_routes,
            },
        },
        "equilibrium": {
            "derivative_policy": {
                "accepted_derivative_backends": [
                    "auto",
                    "analytic",
                    "cppad",
                    "analytic_implicit",
                    "cppad_implicit",
                ],
                "auto_policy": "analytic_or_cppad_or_implicit_else_raise",
            },
            **equilibrium_route_capabilities,
            "repeated_state_properties": {
                "available": True,
                "helpers": ["evaluate_fugacity_coefficients", "evaluate_fugacity_coefficients_batch"],
                "density_seed_parameter": "rho_guess",
            },
            "problem_objects": {
                "available": True,
                "backend": "public_python_facade",
                "classes": list(EQUILIBRIUM_PROBLEM_OBJECT_CLASSES),
                "entrypoint": "mixture.solve_equilibrium(problem)",
            },
            "contribution_maps": {
                "available": True,
                "backend": "native_term_payloads",
                "families": ["hard_chain", "dispersion", "association", "ionic", "born"],
                "inactive_terms_explicit": True,
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
                "backend": "structured_residual_evaluation",
                "scope": "fixed-shape residual evaluator for native thermodynamic calls, not a production optimizer",
            },
            # AlgID: reactive_electrolyte_batch_residual_context
            "reactive_electrolyte_batch_context": {
                "available": True,
                "backend": "batch_residual_evaluation_context",
                "mixed_pressure_speciation_residual_context": {
                    "available": True,
                    "production_optimizer": False,
                    "optimizer": None,
                    "fit_role": "diagnostic residual context, not a production optimizer",
                    "supports_pressure_targets": True,
                    "supports_speciation_targets": True,
                    "supports_activity_targets": True,
                    "supports_fugacity_targets": True,
                    "supports_density_targets": True,
                    "supports_relative_permittivity_targets": True,
                    "validates_parameter_bounds": True,
                    "native_hot_loop": False,
                    "ceres": {
                        "available": bool(ceres["available"]),
                        "production": False,
                        "reason": "native Ceres reactive batch optimizer not registered",
                    },
                    "thermodynamic_backend": "native",
                    "python_role": "row orchestration and diagnostics",
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
                ],
            },
        },
    }


__git_commit__ = str(runtime_build_info()["source_git_commit"])
