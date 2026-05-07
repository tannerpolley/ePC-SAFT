"""Optional cyipopt adapter surface for opt-in NLP solvers."""

from __future__ import annotations

import importlib
import os
from pathlib import Path
from typing import Any

import numpy as np

from ._types import InputError

_IPOPT_DLL_HANDLES: list[Any] = []
_IPOPT_DLL_PATHS_PREPARED: set[str] = set()
_IPOPT_PATHS_PREPENDED: set[str] = set()


def _candidate_ipopt_dll_dirs() -> list[Path]:
    candidates: list[Path] = []
    explicit = os.environ.get("EPCSAFT_IPOPT_DLL_DIR")
    if explicit:
        candidates.append(Path(explicit))
    ipopt_win_dir = os.environ.get("IPOPTWINDIR")
    if ipopt_win_dir:
        candidates.append(Path(ipopt_win_dir) / "bin")
    conda_prefix = os.environ.get("CONDA_PREFIX")
    if conda_prefix:
        candidates.append(Path(conda_prefix) / "Library" / "bin")
    return candidates


def _prepare_ipopt_dll_search_path() -> None:
    add_dll_directory = getattr(os, "add_dll_directory", None)
    for candidate in _candidate_ipopt_dll_dirs():
        try:
            resolved = candidate.resolve()
        except OSError:
            continue
        key = str(resolved)
        if not resolved.is_dir():
            continue
        if key not in _IPOPT_PATHS_PREPENDED:
            os.environ["PATH"] = key + os.pathsep + os.environ.get("PATH", "")
            _IPOPT_PATHS_PREPENDED.add(key)
        if add_dll_directory is not None and key not in _IPOPT_DLL_PATHS_PREPARED:
            _IPOPT_DLL_HANDLES.append(add_dll_directory(key))
            _IPOPT_DLL_PATHS_PREPARED.add(key)


def _cyipopt_import_error() -> str:
    _prepare_ipopt_dll_search_path()
    try:
        importlib.import_module("cyipopt")
    except Exception as exc:  # pragma: no cover - depends on local optional install state.
        return f"{type(exc).__name__}: {exc}"
    return ""


def cyipopt_available() -> bool:
    """Return whether the optional cyipopt module can be imported."""

    return _cyipopt_import_error() == ""


def cyipopt_backend_info() -> dict[str, Any]:
    """Return JSON-like metadata for the optional cyipopt backend."""

    error = _cyipopt_import_error()
    return {
        "available": error == "",
        "backend": "cyipopt",
        "module": "cyipopt",
        "optional_extra": "ipopt",
        "dependency_group": "ipopt",
        "import_error": error,
    }


def require_cyipopt(route: str) -> None:
    """Raise a package InputError when an explicit IPOPT request cannot run."""

    if cyipopt_available():
        return
    info = cyipopt_backend_info()
    reason = str(info.get("import_error") or "cyipopt is not installed")
    raise InputError(
        "cyipopt is required when solver_backend='ipopt' is requested for {}. "
        "Install the optional IPOPT dependency with the package 'ipopt' extra or the uv 'ipopt' group. "
        "The requested solve was not run with the Newton/default backend. Import diagnostic: {}".format(
            route,
            reason,
        )
    )


def unsupported_ipopt_route(route: str) -> None:
    """Reject routes that have not yet been wired to native residual callbacks."""

    require_cyipopt(route)
    raise InputError(
        "solver_backend='ipopt' is opt-in but is not wired for {} yet. "
        "The cyipopt adapter is available, but this route still needs native residual/Jacobian callbacks before it "
        "can run without reimplementing thermodynamics in Python.".format(route)
    )


class _NativeLeastSquaresProblem:
    def __init__(self, evaluator: Any, hessian_strategy: str) -> None:
        self._evaluator = evaluator
        self._hessian_strategy = str(hessian_strategy)
        self._cache_x: np.ndarray | None = None
        self._cache_payload: dict[str, Any] | None = None
        self.evaluation_count = 0

    def _evaluate(self, x: Any) -> dict[str, Any]:
        x_array = np.asarray(x, dtype=float)
        if self._cache_x is not None and np.array_equal(x_array, self._cache_x):
            assert self._cache_payload is not None
            return self._cache_payload
        payload = self._evaluator(x_array)
        self._cache_x = np.array(x_array, dtype=float, copy=True)
        self._cache_payload = payload
        self.evaluation_count += 1
        return payload

    def objective(self, x: Any) -> float:
        return float(self._evaluate(x)["objective"])

    def gradient(self, x: Any) -> np.ndarray:
        return np.asarray(self._evaluate(x)["gradient"], dtype=float)

    def constraints(self, x: Any) -> np.ndarray:
        _ = x
        return np.asarray([], dtype=float)

    def jacobian(self, x: Any) -> np.ndarray:
        _ = x
        return np.asarray([], dtype=float)

    def hessianstructure(self) -> tuple[np.ndarray, np.ndarray]:
        n = int(np.asarray(self._evaluate_initial_variables(), dtype=float).size)
        rows: list[int] = []
        cols: list[int] = []
        for row in range(n):
            for col in range(row + 1):
                rows.append(row)
                cols.append(col)
        return np.asarray(rows, dtype=int), np.asarray(cols, dtype=int)

    def hessian(self, x: Any, lagrange: Any, objective_factor: float) -> np.ndarray:
        _ = lagrange
        payload = self._evaluate(x)
        shape = tuple(int(v) for v in payload["jacobian_shape"])
        jac = np.asarray(payload["jacobian_row_major"], dtype=float).reshape(shape)
        hess = float(objective_factor) * (jac.T @ jac)
        rows, cols = self.hessianstructure()
        return np.asarray([hess[int(r), int(c)] for r, c in zip(rows, cols)], dtype=float)

    def _evaluate_initial_variables(self) -> np.ndarray:
        payload = self._cache_payload
        if payload is None:
            raise RuntimeError("initial variables are not cached yet")
        return np.asarray(payload["variables"], dtype=float)


def _solve_bound_constrained_least_squares(
    *,
    route: str,
    initial_payload: dict[str, Any],
    evaluator: Any,
    max_iterations: int,
    tolerance: float,
    hessian_strategy: str,
) -> dict[str, Any]:
    require_cyipopt(route)
    cyipopt = importlib.import_module("cyipopt")
    x0 = np.asarray(initial_payload["variables"], dtype=float)
    lb = np.asarray(initial_payload["lower_bounds"], dtype=float)
    ub = np.asarray(initial_payload["upper_bounds"], dtype=float)
    problem_obj = _NativeLeastSquaresProblem(evaluator, hessian_strategy)
    problem_obj._cache_x = np.array(x0, dtype=float, copy=True)
    problem_obj._cache_payload = initial_payload
    problem = cyipopt.Problem(
        n=int(x0.size),
        m=0,
        problem_obj=problem_obj,
        lb=lb,
        ub=ub,
        cl=np.asarray([], dtype=float),
        cu=np.asarray([], dtype=float),
    )
    problem.add_option("print_level", 0)
    problem.add_option("max_iter", int(max_iterations))
    problem.add_option("tol", float(tolerance))
    if hessian_strategy == "lbfgs":
        problem.add_option("hessian_approximation", "limited-memory")
    try:
        x_opt, info = problem.solve(x0)
    finally:
        if hasattr(problem, "close"):
            problem.close()
    info_dict = dict(info)
    status_msg = info_dict.get("status_msg", "")
    if isinstance(status_msg, bytes):
        info_dict["status_msg"] = status_msg.decode(errors="replace")
    final_payload = evaluator(np.asarray(x_opt, dtype=float))
    return {
        "x": np.asarray(x_opt, dtype=float),
        "info": info_dict,
        "evaluation": final_payload,
        "initial_evaluation": initial_payload,
        "evaluation_count": int(problem_obj.evaluation_count),
        "hessian_backend": "lbfgs" if hessian_strategy == "lbfgs" else "gauss_newton",
    }


def _ipopt_status_success(info: dict[str, Any]) -> bool:
    try:
        return int(info.get("status", 999)) in {0, 1}
    except Exception:
        return False


def _hessian_metadata(strategy: str) -> dict[str, Any]:
    if strategy == "lbfgs":
        return {
            "exact_hessian_available": False,
            "hessian_strategy": "lbfgs",
            "hessian_kind": "ipopt_limited_memory",
            "hessian_includes_second_residual_derivatives": False,
            "hessian_structure": "ipopt_internal_limited_memory",
            "hessian_callback_available": False,
            "sparse_hessian_available": False,
            "hessian_available": False,
        }
    return {
        "exact_hessian_available": False,
        "hessian_strategy": "gauss_newton",
        "hessian_kind": "approximate_least_squares_gauss_newton",
        "hessian_includes_second_residual_derivatives": False,
        "hessian_structure": "dense_lower_triangular",
        "hessian_callback_available": True,
        "sparse_hessian_available": False,
        "hessian_available": True,
    }


def _solver_selection_metadata(route: str, selected_backend: str, reason: str) -> dict[str, Any]:
    return {
        "solver_route": route,
        "selected_solver_backend": selected_backend,
        "solver_selection_reason": reason,
        "formulation": "bound_constrained_residual_minimization",
        "full_constrained_nlp": False,
    }


def _electrolyte_lle_native_seed(
    *,
    mixture: Any,
    T: float,
    P: float,
    feed: Any,
    species: list[str],
    options: Any,
) -> tuple[dict[str, Any] | None, str]:
    from . import _core

    request = {
        "kind": "electrolyte_lle",
        "T": float(T),
        "P": float(P),
        "z": np.asarray(feed, dtype=float).tolist(),
        "species": species,
        "initial_phases": None,
        "options": {
            "max_iterations": int(options.max_iterations),
            "tolerance": float(options.tolerance),
            "damping": float(options.damping),
            "min_composition": float(options.min_composition),
            "include_phase_diagnostics": bool(options.include_phase_diagnostics),
            "stability_precheck": bool(options.stability_precheck),
            "density_diagnostics": str(options.density_diagnostics),
            "experimental_coupled_density_lle": bool(options.experimental_coupled_density_lle),
            "jacobian_backend": str(options.jacobian_backend),
            "solver_backend": "newton",
            "hessian_strategy": str(options.hessian_strategy),
        },
    }
    try:
        payload = _core._solve_equilibrium_native(mixture._native, request)
    except Exception as exc:
        return None, f"{type(exc).__name__}: {exc}"
    phases = list(payload.get("phases") or [])
    if len(phases) != 2:
        return None, "native seed did not return two phases"
    by_label = {str(item.get("label", "")): item for item in phases}
    aq = by_label.get("aq", phases[0])
    org = by_label.get("org", phases[1])
    return {
        "aq": list(aq["composition"]),
        "org": list(org["composition"]),
        "phase_fraction": float(org["phase_fraction"]),
    }, "native_transformed_newton_seed"


def solve_electrolyte_lle_ipopt(**kwargs: Any) -> Any:
    """Solve electrolyte LLE as a bound-constrained residual minimization with cyipopt."""

    require_cyipopt("electrolyte_lle")
    from . import _core
    from ._types import SolutionError
    from .equilibrium import EquilibriumPhase, EquilibriumResult

    mixture = kwargs["mixture"]
    options = kwargs["options"]
    species = list(getattr(mixture, "species", []))
    initial_phases = kwargs.get("initial_phases")
    seed_source = "user_initial_phases" if initial_phases is not None else "native_transformed_newton_seed"
    seed_failure = ""
    if initial_phases is None:
        initial_phases, seed_failure = _electrolyte_lle_native_seed(
            mixture=mixture,
            T=float(kwargs["T"]),
            P=float(kwargs["P"]),
            feed=kwargs["feed"],
            species=species,
            options=options,
        )
        if initial_phases is None:
            seed_source = "native_residual_default_seed"
    request = {
        "T": float(kwargs["T"]),
        "P": float(kwargs["P"]),
        "z": np.asarray(kwargs["feed"], dtype=float).tolist(),
        "species": species,
        "initial_phases": initial_phases,
        "options": {
            "max_iterations": int(options.max_iterations),
            "tolerance": float(options.tolerance),
            "damping": float(options.damping),
            "min_composition": float(options.min_composition),
            "include_phase_diagnostics": bool(options.include_phase_diagnostics),
            "stability_precheck": bool(options.stability_precheck),
            "density_diagnostics": str(options.density_diagnostics),
            "experimental_coupled_density_lle": bool(options.experimental_coupled_density_lle),
            "jacobian_backend": str(options.jacobian_backend),
        },
    }

    def evaluate(x: np.ndarray) -> dict[str, Any]:
        payload_request = dict(request)
        payload_request["variables"] = np.asarray(x, dtype=float).tolist()
        return _core._evaluate_electrolyte_lle_residual_native(mixture._native, payload_request)

    initial_payload = _core._evaluate_electrolyte_lle_residual_native(mixture._native, request)
    solve = _solve_bound_constrained_least_squares(
        route="electrolyte_lle",
        initial_payload=initial_payload,
        evaluator=evaluate,
        max_iterations=int(options.max_iterations),
        tolerance=float(options.tolerance),
        hessian_strategy=str(options.hessian_strategy),
    )
    final = solve["evaluation"]
    residual_norm = float(max((abs(v) for v in final["residual"]), default=0.0))
    material_error = float(final["material_balance_error"])
    charge_error = float(final["charge_balance_error"])
    phase_distance = float(final["phase_distance"])
    residual_gate_success = residual_norm <= float(options.tolerance)
    physical_gate_success = (
        material_error <= max(float(options.tolerance), 1.0e-10)
        and charge_error <= 1.0e-8
        and phase_distance > max(1.0e-4, 100.0 * float(options.min_composition))
    )
    accepted = bool(residual_gate_success and physical_gate_success)
    ipopt_success = _ipopt_status_success(solve["info"])
    diagnostics = dict(final.get("diagnostics") or {})
    diagnostics.update(
        {
            "solver_method": "cyipopt_bound_min_residual",
            "solver_language": "python_cyipopt_native_callbacks",
            "native_entrypoint": "_evaluate_electrolyte_lle_residual_native",
            "requested_solver_backend": "ipopt",
            "hessian_backend": solve["hessian_backend"],
            "ipopt_status": int(solve["info"].get("status", 999)),
            "ipopt_status_msg": str(solve["info"].get("status_msg", "")),
            "ipopt_success": ipopt_success,
            "ipopt_objective": float(final["objective"]),
            "ipopt_evaluation_count": int(solve["evaluation_count"]),
            "solver_residual_norm": residual_norm,
            "fugacity_residual_norm": residual_norm,
            "material_balance_error": material_error,
            "charge_balance_error": charge_error,
            "phase_distance": phase_distance,
            "residual_gate_success": bool(residual_gate_success),
            "physical_gate_success": bool(physical_gate_success),
            "accepted": accepted,
            "accepted_despite_ipopt_status": bool(accepted and not ipopt_success),
            "acceptance_gate": "ipopt_min_residual" if accepted else "ipopt_min_residual_failed",
            "variable_model": str(final["variable_model"]),
            "ipopt_seed_source": seed_source,
            "ipopt_seed_failure": seed_failure,
        }
    )
    diagnostics.update(_hessian_metadata(str(solve["hessian_backend"])))
    diagnostics.update(_solver_selection_metadata("electrolyte_lle", "ipopt", "explicit_request"))
    diagnostics.update(dict(kwargs.get("feed_diagnostics") or {}))
    if not accepted:
        raise SolutionError("electrolyte LLE IPOPT solve did not satisfy acceptance gates", diagnostics)
    aq = EquilibriumPhase(
        label="aq",
        composition=final["aq_composition"],
        density=final["aq_density"],
        temperature=float(kwargs["T"]),
        pressure=float(kwargs["P"]),
        phase_fraction=1.0 - float(final["phase_fraction_org"]),
        ln_fugacity_coefficient=final["aq_ln_fugacity_coefficient"],
    )
    org = EquilibriumPhase(
        label="org",
        composition=final["org_composition"],
        density=final["org_density"],
        temperature=float(kwargs["T"]),
        pressure=float(kwargs["P"]),
        phase_fraction=float(final["phase_fraction_org"]),
        ln_fugacity_coefficient=final["org_ln_fugacity_coefficient"],
    )
    return EquilibriumResult(
        backend="electrolyte_lle_ipopt",
        problem_kind="electrolyte_lle_flash",
        phases=(aq, org),
        stable=False,
        split_detected=True,
        diagnostics=diagnostics,
    )


def solve_reactive_speciation_ipopt(**kwargs: Any) -> Any:
    """Solve reactive speciation as a bound-constrained residual minimization with cyipopt."""

    require_cyipopt("reactive_speciation")
    from . import _core
    from .reactive_speciation import (
        _REACTION_STANDARD_STATES,
        _json_like,
        _named_reaction_residuals,
        _reaction_standard_state_summary,
        ReactiveSpeciationResult,
    )

    species = list(kwargs["species"])
    mixture_factory = kwargs["mixture_factory"]
    initial_x = np.asarray(kwargs["initial_x"], dtype=float)
    options = kwargs["options"]
    mixture = mixture_factory(initial_x, float(kwargs["T"]), float(kwargs["P"]))
    native = getattr(mixture, "_native", None)
    if native is None:
        raise InputError("IPOPT reactive speciation requires mixture_factory to return an ePCSAFTMixture.")
    reaction_matrix = np.asarray(
        [[float(reaction.stoichiometry.get(label, 0.0)) for label in species] for reaction in kwargs["reactions"]],
        dtype=float,
    )
    request = {
        "T": float(kwargs["T"]),
        "P": float(kwargs["P"]),
        "initial_x": initial_x.tolist(),
        "balance_matrix": np.asarray(kwargs["balance_matrix"], dtype=float).reshape(-1).tolist(),
        "balance_rows": int(kwargs["balance_matrix"].shape[0]),
        "total_vector": np.asarray(kwargs["total_vector"], dtype=float).tolist(),
        "reaction_stoichiometry": reaction_matrix.reshape(-1).tolist(),
        "reaction_rows": int(reaction_matrix.shape[0]),
        "log_equilibrium_constants": [float(reaction.log_equilibrium_constant) for reaction in kwargs["reactions"]],
        "reaction_standard_states": [
            int(_REACTION_STANDARD_STATES[reaction.standard_state]) for reaction in kwargs["reactions"]
        ],
        "options": {
            "max_iterations": int(options.max_iterations),
            "tolerance": float(options.tolerance),
            "damping": float(options.damping),
            "min_mole_fraction": float(options.min_mole_fraction),
            "finite_difference_step": float(options.finite_difference_step),
            "jacobian_backend": str(options.jacobian_backend),
            "phase": str(options.phase),
            "activity_output": str(options.activity_output),
        },
    }

    def evaluate(x: np.ndarray) -> dict[str, Any]:
        payload_request = dict(request)
        payload_request["variables"] = np.asarray(x, dtype=float).tolist()
        return _core._evaluate_chemical_equilibrium_residual_native(native, payload_request)

    initial_payload = _core._evaluate_chemical_equilibrium_residual_native(native, request)
    solve = _solve_bound_constrained_least_squares(
        route="reactive_speciation",
        initial_payload=initial_payload,
        evaluator=evaluate,
        max_iterations=int(options.max_iterations),
        tolerance=float(options.tolerance),
        hessian_strategy=str(options.hessian_strategy),
    )
    final = solve["evaluation"]
    x = {label: float(value) for label, value in zip(species, final["composition"])}
    activity_coefficients = {label: float(value) for label, value in zip(species, final["activity_coefficients"])}
    mass_balance_residuals = {
        name: float(value) for name, value in zip(kwargs["balance_names"], final["mass_balance_residuals"])
    }
    reaction_residuals = [float(value) for value in final["reaction_residuals"]]
    named_reaction_residuals = _named_reaction_residuals(kwargs["reactions"], reaction_residuals)
    mass_tolerance = options.mass_tolerance if options.mass_tolerance is not None else options.tolerance
    charge_tolerance = options.charge_tolerance if options.charge_tolerance is not None else options.tolerance
    reaction_tolerance = options.reaction_tolerance if options.reaction_tolerance is not None else options.tolerance
    mass_residual_norm = float(max((abs(value) for value in mass_balance_residuals.values()), default=0.0))
    charge_residual = float(final["charge_residual"])
    reaction_residual_norm = float(max((abs(value) for value in reaction_residuals), default=0.0))
    residual_family_success = (
        mass_residual_norm <= mass_tolerance
        and abs(charge_residual) <= charge_tolerance
        and reaction_residual_norm <= reaction_tolerance
    )
    mass_balance_gate_success = mass_residual_norm <= mass_tolerance
    charge_balance_gate_success = abs(charge_residual) <= charge_tolerance
    reaction_residual_gate_success = reaction_residual_norm <= reaction_tolerance
    ipopt_success = _ipopt_status_success(solve["info"])
    diagnostics = dict(final.get("diagnostics") or {})
    activity_basis = _reaction_standard_state_summary(kwargs["reactions"])
    handoff = dict(diagnostics.get("phase_equilibrium_handoff", {}))
    handoff["composition"] = [float(value) for value in final["composition"]]
    handoff["activity_coefficients"] = [float(value) for value in final["activity_coefficients"]]
    handoff["composition_map"] = dict(x)
    handoff["activity_coefficients_map"] = dict(activity_coefficients)
    handoff["activity_basis"] = activity_basis
    diagnostics["phase_equilibrium_handoff"] = handoff
    diagnostics.update(
        {
            "activity_basis": activity_basis,
            "success": bool(residual_family_success),
            "native_success": bool(residual_family_success),
            "residual_family_success": bool(residual_family_success),
            "message": "converged" if residual_family_success else "reactive speciation IPOPT residual tolerances were not met",
            "backend": "ipopt",
            "solver_method": "cyipopt_bound_min_residual",
            "solver_language": "python_cyipopt_native_callbacks",
            "native_entrypoint": "_evaluate_chemical_equilibrium_residual_native",
            "requested_solver_backend": "ipopt",
            "hessian_backend": solve["hessian_backend"],
            "ipopt_status": int(solve["info"].get("status", 999)),
            "ipopt_status_msg": str(solve["info"].get("status_msg", "")),
            "ipopt_success": ipopt_success,
            "ipopt_objective": float(final["objective"]),
            "ipopt_evaluation_count": int(solve["evaluation_count"]),
            "mass_tolerance": float(mass_tolerance),
            "charge_tolerance": float(charge_tolerance),
            "reaction_tolerance": float(reaction_tolerance),
            "mass_residual_norm": mass_residual_norm,
            "charge_residual_abs": abs(charge_residual),
            "reaction_residual_norm": reaction_residual_norm,
            "mass_balance_gate_success": bool(mass_balance_gate_success),
            "charge_balance_gate_success": bool(charge_balance_gate_success),
            "reaction_residual_gate_success": bool(reaction_residual_gate_success),
            "residual_gate_success": bool(residual_family_success),
            "physical_gate_success": bool(residual_family_success),
            "accepted": bool(residual_family_success),
            "accepted_despite_ipopt_status": bool(residual_family_success and not ipopt_success),
            "named_reaction_residuals": dict(named_reaction_residuals),
            "best_x": dict(x),
            "best_activity_coefficients": dict(activity_coefficients),
        }
    )
    diagnostics.update(_hessian_metadata(str(solve["hessian_backend"])))
    diagnostics.update(_solver_selection_metadata("reactive_speciation", "ipopt", "explicit_request"))
    result = ReactiveSpeciationResult(
        success=bool(residual_family_success),
        message=str(diagnostics["message"]),
        x=x,
        activity_coefficients=activity_coefficients,
        mass_balance_residuals=mass_balance_residuals,
        charge_residual=charge_residual,
        reaction_residuals=reaction_residuals,
        named_reaction_residuals=named_reaction_residuals,
        state_failure_count=int(diagnostics.get("state_failure_count", 0)),
        diagnostics=diagnostics,
    )
    if not result.success and not options.return_best_effort:
        from ._types import SolutionError

        raise SolutionError(result.message, _json_like(diagnostics))
    return result
