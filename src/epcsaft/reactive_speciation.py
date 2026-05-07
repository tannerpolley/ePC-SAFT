"""Homogeneous reactive speciation helpers using ePC-SAFT activity states."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from dataclasses import field
from typing import Any

import numpy as np

from ._types import InputError, SolutionError

_REACTION_STANDARD_STATES = {
    "mole_fraction_activity": 0,
    "ideal_mole_fraction": 1,
    "concentration": 2,
}


@dataclass(frozen=True, slots=True)
class ReactionDefinition:
    """One reaction residual definition for reactive speciation."""

    stoichiometry: Mapping[str, float]
    log_equilibrium_constant: float
    name: str = ""
    standard_state: str = "mole_fraction_activity"

    def __post_init__(self) -> None:
        object.__setattr__(self, "stoichiometry", {str(k): float(v) for k, v in self.stoichiometry.items()})
        object.__setattr__(self, "log_equilibrium_constant", float(self.log_equilibrium_constant))
        object.__setattr__(self, "name", str(self.name))
        standard_state = str(self.standard_state).strip().lower()
        if standard_state not in _REACTION_STANDARD_STATES:
            supported = "', '".join(_REACTION_STANDARD_STATES)
            raise InputError(f"ReactionDefinition.standard_state must be one of '{supported}'.")
        object.__setattr__(self, "standard_state", standard_state)


@dataclass(frozen=True, slots=True)
class ReactiveSpeciationOptions:
    """Numerical controls for homogeneous reactive speciation."""

    max_iterations: int = 50
    tolerance: float = 1.0e-8
    damping: float = 0.5
    min_mole_fraction: float = 1.0e-14
    finite_difference_step: float = 1.0e-6
    jacobian_backend: str = "auto"
    solver_backend: str = "auto"
    hessian_strategy: str = "gauss_newton"
    phase: str = "liq"
    return_best_effort: bool = False
    error_mode: str = "raise"
    activity_output: str = "auto"
    mass_tolerance: float | None = None
    charge_tolerance: float | None = None
    reaction_tolerance: float | None = None


@dataclass(frozen=True, slots=True)
class ReactiveSpeciationResult:
    """Structured result returned by reactive speciation solves."""

    success: bool
    message: str
    x: dict[str, float]
    activity_coefficients: dict[str, float]
    mass_balance_residuals: dict[str, float]
    charge_residual: float
    reaction_residuals: list[float]
    named_reaction_residuals: dict[str, float]
    state_failure_count: int
    diagnostics: dict[str, Any]
    continuation_state: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        object.__setattr__(self, "success", bool(self.success))
        object.__setattr__(self, "message", str(self.message))
        object.__setattr__(self, "x", {str(k): float(v) for k, v in self.x.items()})
        object.__setattr__(
            self, "activity_coefficients", {str(k): float(v) for k, v in self.activity_coefficients.items()}
        )
        object.__setattr__(
            self, "mass_balance_residuals", {str(k): float(v) for k, v in self.mass_balance_residuals.items()}
        )
        object.__setattr__(self, "charge_residual", float(self.charge_residual))
        object.__setattr__(self, "reaction_residuals", [float(v) for v in self.reaction_residuals])
        object.__setattr__(
            self,
            "named_reaction_residuals",
            {str(k): float(v) for k, v in self.named_reaction_residuals.items()},
        )
        object.__setattr__(self, "state_failure_count", int(self.state_failure_count))
        object.__setattr__(self, "diagnostics", dict(self.diagnostics))
        object.__setattr__(self, "continuation_state", _json_like(dict(self.continuation_state)))

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-like result payload."""
        return {
            "success": self.success,
            "message": self.message,
            "x": dict(self.x),
            "activity_coefficients": dict(self.activity_coefficients),
            "mass_balance_residuals": dict(self.mass_balance_residuals),
            "charge_residual": self.charge_residual,
            "reaction_residuals": list(self.reaction_residuals),
            "named_reaction_residuals": dict(self.named_reaction_residuals),
            "state_failure_count": self.state_failure_count,
            "diagnostics": _json_like(self.diagnostics),
            "continuation_state": _json_like(self.continuation_state),
        }


def solve_reactive_speciation(
    *,
    species: Any,
    mixture_factory: Any,
    T: float,
    P: float,
    balances: Mapping[str, Mapping[str, float]],
    totals: Mapping[str, float],
    reactions: Any,
    initial_x: Any,
    options: ReactiveSpeciationOptions | None = None,
    warm_start: Any = None,
) -> ReactiveSpeciationResult:
    """Solve homogeneous activity-coupled reactive speciation."""
    opts = _normalize_options(options)
    labels = [str(label) for label in species]
    if not labels:
        raise InputError("species must include at least one label.")
    initial, initial_source = _initial_composition_from_inputs(
        initial_x=initial_x,
        warm_start=warm_start,
        ncomp=len(labels),
        min_value=opts.min_mole_fraction,
    )
    balance_matrix, total_vector, balance_names = _normalize_balances(labels, balances, totals)
    reaction_defs = _normalize_reactions(labels, reactions)
    if opts.solver_backend == "ipopt":
        from .ipopt_backend import solve_reactive_speciation_ipopt

        return solve_reactive_speciation_ipopt(
            species=labels,
            mixture_factory=mixture_factory,
            T=T,
            P=P,
            balance_matrix=balance_matrix,
            total_vector=total_vector,
            balance_names=balance_names,
            reactions=reaction_defs,
            initial_x=initial,
            initial_x_source=initial_source,
            options=opts,
        )
    return _solve_reactive_speciation_native(
        species=labels,
        mixture_factory=mixture_factory,
        T=T,
        P=P,
        balance_matrix=balance_matrix,
        total_vector=total_vector,
        balance_names=balance_names,
        reactions=reaction_defs,
        initial_x=initial,
        initial_x_source=initial_source,
        options=opts,
    )


def solve_reactive_speciation_sweep(
    *,
    species: Any,
    mixture_factory: Any,
    points: Any,
    balances: Mapping[str, Mapping[str, float]],
    reactions: Any,
    options: ReactiveSpeciationOptions | None = None,
    continuation: str = "auto",
) -> list[ReactiveSpeciationResult]:
    """Solve an ordered reactive-speciation sweep with fixed-shape results."""

    opts = _normalize_options(options)
    mode = str(continuation).strip().lower()
    if mode not in {"auto", "none"}:
        raise InputError("continuation must be 'auto' or 'none'.")
    labels = [str(label) for label in species]
    if not labels:
        raise InputError("species must include at least one label.")
    reaction_defs = _normalize_reactions(labels, reactions)
    results: list[ReactiveSpeciationResult] = []
    warm_start: dict[str, Any] | None = None
    for index, point in enumerate(points):
        try:
            if "T" not in point or "P" not in point or "totals" not in point:
                raise InputError("Each reactive speciation sweep point requires T, P, and totals.")
            point_warm_start = warm_start if mode == "auto" else None
            result = solve_reactive_speciation(
                species=labels,
                mixture_factory=mixture_factory,
                T=float(point["T"]),
                P=float(point["P"]),
                balances=balances,
                totals=point["totals"],
                reactions=reaction_defs,
                initial_x=point.get("initial_x"),
                options=opts,
                warm_start=point_warm_start,
            )
        except Exception as exc:
            if opts.error_mode != "result":
                raise
            result = _structured_failure_result(
                species=labels,
                point=point,
                index=index,
                exc=exc,
                options=opts,
            )
        results.append(result)
        if result.success:
            warm_start = dict(result.continuation_state)
    return results


def _normalize_options(options: ReactiveSpeciationOptions | None) -> ReactiveSpeciationOptions:
    if options is None:
        return ReactiveSpeciationOptions()
    if not isinstance(options, ReactiveSpeciationOptions):
        raise InputError("options must be a ReactiveSpeciationOptions instance.")
    if options.max_iterations < 0:
        raise InputError("ReactiveSpeciationOptions.max_iterations must be non-negative.")
    if options.tolerance <= 0.0 or options.min_mole_fraction <= 0.0 or options.finite_difference_step <= 0.0:
        raise InputError("ReactiveSpeciationOptions tolerances and steps must be positive.")
    if not (0.0 < options.damping <= 1.0):
        raise InputError("ReactiveSpeciationOptions.damping must be in (0, 1].")
    if not isinstance(options.return_best_effort, bool):
        raise InputError("ReactiveSpeciationOptions.return_best_effort must be a bool.")
    error_mode = str(options.error_mode).strip().lower()
    if error_mode not in {"raise", "result"}:
        raise InputError("ReactiveSpeciationOptions.error_mode must be 'raise' or 'result'.")
    activity_output = str(options.activity_output).strip().lower()
    if activity_output not in {"auto", "always", "never"}:
        raise InputError("ReactiveSpeciationOptions.activity_output must be 'auto', 'always', or 'never'.")
    return_best_effort = bool(options.return_best_effort or error_mode == "result")
    jacobian_backend = str(options.jacobian_backend).strip().lower()
    if jacobian_backend in {"numerical", "fd"}:
        jacobian_backend = "finite_difference"
    if jacobian_backend not in {"auto", "autodiff", "finite_difference"}:
        raise InputError(
            "ReactiveSpeciationOptions.jacobian_backend must be 'auto', 'autodiff', or 'finite_difference'."
        )
    solver_backend = str(options.solver_backend).strip().lower()
    if solver_backend not in {"auto", "newton", "ipopt"}:
        raise InputError("ReactiveSpeciationOptions.solver_backend must be 'auto', 'newton', or 'ipopt'.")
    hessian_strategy = str(options.hessian_strategy).strip().lower()
    hessian_aliases = {"gn": "gauss_newton", "gauss-newton": "gauss_newton", "bfgs": "lbfgs"}
    hessian_strategy = hessian_aliases.get(hessian_strategy, hessian_strategy)
    if hessian_strategy not in {"gauss_newton", "lbfgs"}:
        raise InputError("ReactiveSpeciationOptions.hessian_strategy must be 'gauss_newton' or 'lbfgs'.")
    for name in ("mass_tolerance", "charge_tolerance", "reaction_tolerance"):
        value = getattr(options, name)
        if value is not None and value <= 0.0:
            raise InputError(f"ReactiveSpeciationOptions.{name} must be positive when provided.")
    if (
        jacobian_backend == options.jacobian_backend
        and solver_backend == options.solver_backend
        and hessian_strategy == options.hessian_strategy
        and error_mode == options.error_mode
        and activity_output == options.activity_output
        and return_best_effort == options.return_best_effort
    ):
        return options
    return ReactiveSpeciationOptions(
        max_iterations=options.max_iterations,
        tolerance=options.tolerance,
        damping=options.damping,
        min_mole_fraction=options.min_mole_fraction,
        finite_difference_step=options.finite_difference_step,
        jacobian_backend=jacobian_backend,
        solver_backend=solver_backend,
        hessian_strategy=hessian_strategy,
        phase=options.phase,
        return_best_effort=return_best_effort,
        error_mode=error_mode,
        activity_output=activity_output,
        mass_tolerance=options.mass_tolerance,
        charge_tolerance=options.charge_tolerance,
        reaction_tolerance=options.reaction_tolerance,
    )


def _solve_reactive_speciation_native(
    *,
    species: list[str],
    mixture_factory: Any,
    T: float,
    P: float,
    balance_matrix: np.ndarray,
    total_vector: np.ndarray,
    balance_names: list[str],
    reactions: list[ReactionDefinition],
    initial_x: np.ndarray,
    options: ReactiveSpeciationOptions,
    initial_x_source: str = "initial_x",
) -> ReactiveSpeciationResult:
    from . import _core

    mixture = mixture_factory(initial_x, T, P)
    native = getattr(mixture, "_native", None)
    if native is None:
        raise InputError("native reactive speciation backend requires mixture_factory to return an ePCSAFTMixture.")
    if list(getattr(mixture, "species", species)) != species:
        raise InputError("native reactive speciation backend requires mixture species to match the species argument.")
    reaction_matrix = np.asarray(
        [[float(reaction.stoichiometry.get(label, 0.0)) for label in species] for reaction in reactions],
        dtype=float,
    )
    request = {
        "T": float(T),
        "P": float(P),
        "initial_x": np.asarray(initial_x, dtype=float).tolist(),
        "balance_matrix": np.asarray(balance_matrix, dtype=float).reshape(-1).tolist(),
        "balance_rows": int(balance_matrix.shape[0]),
        "total_vector": np.asarray(total_vector, dtype=float).tolist(),
        "reaction_stoichiometry": reaction_matrix.reshape(-1).tolist(),
        "reaction_rows": int(reaction_matrix.shape[0]),
        "log_equilibrium_constants": [float(reaction.log_equilibrium_constant) for reaction in reactions],
        "reaction_standard_states": [int(_REACTION_STANDARD_STATES[reaction.standard_state]) for reaction in reactions],
        "options": {
            "max_iterations": int(options.max_iterations),
            "tolerance": float(options.tolerance),
            "damping": float(options.damping),
            "min_mole_fraction": float(options.min_mole_fraction),
            "finite_difference_step": float(options.finite_difference_step),
            "jacobian_backend": str(options.jacobian_backend),
            "solver_backend": str(options.solver_backend),
            "hessian_strategy": str(options.hessian_strategy),
            "phase": str(options.phase),
            "activity_output": str(options.activity_output),
        },
    }
    payload = _core._solve_chemical_equilibrium_native(native, request)
    x = {label: float(value) for label, value in zip(species, payload["composition"])}
    activity_coefficients = {label: float(value) for label, value in zip(species, payload["activity_coefficients"])}
    mass_balance_residuals = {
        name: float(value) for name, value in zip(balance_names, payload["mass_balance_residuals"])
    }
    reaction_residuals = [float(value) for value in payload["reaction_residuals"]]
    named_reaction_residuals = _named_reaction_residuals(reactions, reaction_residuals)
    mass_tolerance = options.mass_tolerance if options.mass_tolerance is not None else options.tolerance
    charge_tolerance = options.charge_tolerance if options.charge_tolerance is not None else options.tolerance
    reaction_tolerance = options.reaction_tolerance if options.reaction_tolerance is not None else options.tolerance
    mass_residual_norm = float(max((abs(value) for value in mass_balance_residuals.values()), default=0.0))
    reaction_residual_norm = float(max((abs(value) for value in reaction_residuals), default=0.0))
    charge_residual = float(payload["charge_residual"])
    residual_family_success = (
        mass_residual_norm <= mass_tolerance
        and abs(charge_residual) <= charge_tolerance
        and reaction_residual_norm <= reaction_tolerance
    )
    diagnostics = dict(payload["diagnostics"])
    activity_basis = _reaction_standard_state_summary(reactions)
    handoff = dict(diagnostics.get("phase_equilibrium_handoff", {}))
    handoff.setdefault("composition", [float(value) for value in payload["composition"]])
    handoff.setdefault("activity_coefficients", [float(value) for value in payload["activity_coefficients"]])
    handoff["composition_map"] = dict(x)
    handoff["activity_coefficients_map"] = dict(activity_coefficients)
    handoff["activity_basis"] = activity_basis
    diagnostics["phase_equilibrium_handoff"] = handoff
    diagnostics["reaction_standard_states"] = [reaction.standard_state for reaction in reactions]
    diagnostics.update(
        {
            "activity_basis": activity_basis,
            "success": bool(payload["success"] and residual_family_success),
            "native_success": bool(payload["success"]),
            "residual_family_success": bool(residual_family_success),
            "message": str(payload["message"]),
            "backend": "native",
            "activity_output": str(options.activity_output),
            "initial_x_source": str(initial_x_source),
            "continuation_used": str(initial_x_source) != "initial_x",
            "requested_solver_backend": str(options.solver_backend),
            "requested_hessian_strategy": str(options.hessian_strategy),
            "selected_solver_backend": "native",
            "solver_selection_reason": "default_native" if options.solver_backend == "auto" else "explicit_request",
            "default_auto_uses_ipopt": False,
            "mass_tolerance": float(mass_tolerance),
            "charge_tolerance": float(charge_tolerance),
            "reaction_tolerance": float(reaction_tolerance),
            "mass_residual_norm": mass_residual_norm,
            "charge_residual_abs": abs(charge_residual),
            "reaction_residual_norm": reaction_residual_norm,
            "named_reaction_residuals": dict(named_reaction_residuals),
            "best_x": dict(x),
            "best_activity_coefficients": dict(activity_coefficients),
        }
    )
    success = bool(payload["success"] and residual_family_success)
    message = str(payload["message"])
    if bool(payload["success"]) and not residual_family_success:
        message = "reactive speciation residual family tolerances were not met"
        diagnostics["message"] = message
    result = ReactiveSpeciationResult(
        success=success,
        message=message,
        x=x,
        activity_coefficients=activity_coefficients,
        mass_balance_residuals=mass_balance_residuals,
        charge_residual=charge_residual,
        reaction_residuals=reaction_residuals,
        named_reaction_residuals=named_reaction_residuals,
        state_failure_count=int(diagnostics.get("state_failure_count", 0)),
        diagnostics=diagnostics,
        continuation_state=_continuation_state(
            x=x,
            T=T,
            P=P,
            diagnostics=diagnostics,
        ),
    )
    if not success and not options.return_best_effort:
        raise SolutionError(message, _json_like(diagnostics))
    return result


def _named_reaction_residuals(reactions: list[ReactionDefinition], reaction_residuals: list[float]) -> dict[str, float]:
    names: list[str] = []
    counts: dict[str, int] = {}
    for index, reaction in enumerate(reactions):
        base = reaction.name.strip() if reaction.name.strip() else f"reaction_{index}"
        count = counts.get(base, 0)
        counts[base] = count + 1
        names.append(base if count == 0 else f"{base}_{count}")
    return {name: float(value) for name, value in zip(names, reaction_residuals)}


def _reaction_standard_state_summary(reactions: list[ReactionDefinition]) -> str:
    standard_states = [reaction.standard_state for reaction in reactions]
    if not standard_states:
        return "mole_fraction"
    first = standard_states[0]
    if any(value != first for value in standard_states[1:]):
        return "mixed_standard_state"
    if first == "mole_fraction_activity":
        return "mole_fraction"
    return first


def _initial_composition_from_inputs(
    *,
    initial_x: Any,
    warm_start: Any,
    ncomp: int,
    min_value: float,
) -> tuple[np.ndarray, str]:
    if warm_start is not None:
        if isinstance(warm_start, Mapping):
            for key in ("composition", "x"):
                if key in warm_start:
                    value = warm_start[key]
                    if isinstance(value, Mapping):
                        value = list(value.values())
                    return _normalize_composition(value, ncomp, min_value), "previous_successful_result"
        return _normalize_composition(warm_start, ncomp, min_value), "warm_start"
    if initial_x is None:
        raise InputError("initial_x is required when no warm_start composition is supplied.")
    return _normalize_composition(initial_x, ncomp, min_value), "initial_x"


def _continuation_state(*, x: Mapping[str, float], T: float, P: float, diagnostics: Mapping[str, Any]) -> dict[str, Any]:
    composition = {str(label): float(value) for label, value in x.items()}
    return {
        "composition": composition,
        "T": float(T),
        "P": float(P),
        "density_solve_count": int(diagnostics.get("density_solve_count", 0)),
        "activity_evaluation_count": int(diagnostics.get("activity_evaluation_count", 0)),
    }


def _structured_failure_result(
    *,
    species: list[str],
    point: Mapping[str, Any],
    index: int,
    exc: Exception,
    options: ReactiveSpeciationOptions,
) -> ReactiveSpeciationResult:
    try:
        x_array = _normalize_composition(point.get("initial_x"), len(species), options.min_mole_fraction)
    except Exception:
        x_array = np.full(len(species), 1.0 / max(len(species), 1), dtype=float)
    x = {label: float(value) for label, value in zip(species, x_array)}
    diagnostics = {
        "success": False,
        "structured_failure": True,
        "sweep_index": int(index),
        "message": str(exc),
        "exception_type": type(exc).__name__,
        "backend": "python_sweep",
        "selected_solver_backend": "not_run",
        "initial_x_source": "failure_payload",
        "continuation_used": False,
    }
    return ReactiveSpeciationResult(
        success=False,
        message=str(exc),
        x=x,
        activity_coefficients={},
        mass_balance_residuals={},
        charge_residual=0.0,
        reaction_residuals=[],
        named_reaction_residuals={},
        state_failure_count=0,
        diagnostics=diagnostics,
        continuation_state={},
    )


def _normalize_composition(value: Any, ncomp: int, min_value: float) -> np.ndarray:
    x = np.asarray(value, dtype=float).flatten()
    if x.size != int(ncomp):
        raise InputError("initial_x length must match species length.")
    if np.any(~np.isfinite(x)) or np.any(x < -min_value):
        raise InputError("initial_x values must be finite and non-negative.")
    x = np.clip(x, 0.0, None)
    total = float(np.sum(x))
    if total <= 0.0:
        raise InputError("initial_x must have a positive sum.")
    return x / total


def _normalize_balances(
    species: list[str],
    balances: Mapping[str, Mapping[str, float]],
    totals: Mapping[str, float],
) -> tuple[np.ndarray, np.ndarray, list[str]]:
    rows: list[list[float]] = []
    names: list[str] = []
    total_values: list[float] = []
    for name, coeffs in balances.items():
        row = [0.0 for _ in species]
        for label, coeff in coeffs.items():
            if str(label) not in species:
                raise InputError(f"Unknown species '{label}' in balance '{name}'.")
            row[species.index(str(label))] = float(coeff)
        if name not in totals:
            raise InputError(f"Missing total for balance '{name}'.")
        rows.append(row)
        names.append(str(name))
        total_values.append(float(totals[name]))
    if not rows:
        raise InputError("At least one material balance is required.")
    return np.asarray(rows, dtype=float), np.asarray(total_values, dtype=float), names


def _normalize_reactions(species: list[str], reactions: Any) -> list[ReactionDefinition]:
    out: list[ReactionDefinition] = []
    for reaction in reactions:
        if not isinstance(reaction, ReactionDefinition):
            raise InputError("reactions must contain ReactionDefinition instances.")
        for label in reaction.stoichiometry:
            if label not in species:
                raise InputError(f"Unknown species '{label}' in reaction stoichiometry.")
        out.append(reaction)
    return out


def _json_like(value: Any) -> Any:
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, dict):
        return {str(key): _json_like(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_like(item) for item in value]
    if isinstance(value, np.generic):
        return value.item()
    return value
