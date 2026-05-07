"""Homogeneous reactive speciation helpers using ePC-SAFT activity states."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
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
    phase: str = "liq"
    return_best_effort: bool = False
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
) -> ReactiveSpeciationResult:
    """Solve homogeneous activity-coupled reactive speciation."""
    opts = _normalize_options(options)
    labels = [str(label) for label in species]
    if not labels:
        raise InputError("species must include at least one label.")
    initial = _normalize_composition(initial_x, len(labels), opts.min_mole_fraction)
    balance_matrix, total_vector, balance_names = _normalize_balances(labels, balances, totals)
    reaction_defs = _normalize_reactions(labels, reactions)
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
        options=opts,
    )


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
    jacobian_backend = str(options.jacobian_backend).strip().lower()
    if jacobian_backend in {"numerical", "fd"}:
        jacobian_backend = "finite_difference"
    if jacobian_backend not in {"auto", "autodiff", "finite_difference"}:
        raise InputError("ReactiveSpeciationOptions.jacobian_backend must be 'auto', 'autodiff', or 'finite_difference'.")
    for name in ("mass_tolerance", "charge_tolerance", "reaction_tolerance"):
        value = getattr(options, name)
        if value is not None and value <= 0.0:
            raise InputError(f"ReactiveSpeciationOptions.{name} must be positive when provided.")
    if jacobian_backend == options.jacobian_backend:
        return options
    return ReactiveSpeciationOptions(
        max_iterations=options.max_iterations,
        tolerance=options.tolerance,
        damping=options.damping,
        min_mole_fraction=options.min_mole_fraction,
        finite_difference_step=options.finite_difference_step,
        jacobian_backend=jacobian_backend,
        phase=options.phase,
        return_best_effort=options.return_best_effort,
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
        "reaction_standard_states": [
            int(_REACTION_STANDARD_STATES[reaction.standard_state]) for reaction in reactions
        ],
        "options": {
            "max_iterations": int(options.max_iterations),
            "tolerance": float(options.tolerance),
            "damping": float(options.damping),
            "min_mole_fraction": float(options.min_mole_fraction),
            "finite_difference_step": float(options.finite_difference_step),
            "jacobian_backend": str(options.jacobian_backend),
            "phase": str(options.phase),
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
    diagnostics["reaction_standard_states"] = [reaction.standard_state for reaction in reactions]
    diagnostics.update(
        {
            "success": bool(payload["success"] and residual_family_success),
            "native_success": bool(payload["success"]),
            "residual_family_success": bool(residual_family_success),
            "message": str(payload["message"]),
            "backend": "native",
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
