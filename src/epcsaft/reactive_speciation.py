"""Homogeneous reactive speciation helpers using ePC-SAFT activity states."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from typing import Any

import numpy as np

from ._types import InputError, SolutionError


@dataclass(frozen=True, slots=True)
class ReactionDefinition:
    """One reaction residual definition for reactive speciation."""

    stoichiometry: Mapping[str, float]
    log_equilibrium_constant: float
    name: str = ""

    def __post_init__(self) -> None:
        object.__setattr__(self, "stoichiometry", {str(k): float(v) for k, v in self.stoichiometry.items()})
        object.__setattr__(self, "log_equilibrium_constant", float(self.log_equilibrium_constant))
        object.__setattr__(self, "name", str(self.name))


@dataclass(frozen=True, slots=True)
class ReactiveSpeciationOptions:
    """Numerical controls for homogeneous reactive speciation."""

    max_iterations: int = 50
    tolerance: float = 1.0e-8
    damping: float = 0.5
    min_mole_fraction: float = 1.0e-14
    finite_difference_step: float = 1.0e-6
    phase: str = "liq"


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
    log_n = np.log(np.clip(initial, opts.min_mole_fraction, None))
    state_failure_count = 0
    history: list[dict[str, Any]] = []
    best_payload: dict[str, Any] | None = None

    def evaluate(candidate_log_n: np.ndarray) -> dict[str, Any]:
        nonlocal state_failure_count
        n = np.exp(candidate_log_n)
        x = n / float(np.sum(n))
        try:
            state, charges = _state_for(mixture_factory, x, labels, T, P, opts.phase)
            activity = state.activity_coefficient(species=labels)
        except Exception as exc:
            state_failure_count += 1
            raise exc
        gamma = np.asarray([float(activity[label]) for label in labels], dtype=float)
        mass_residual = balance_matrix @ n - total_vector
        charge_residual = 0.0 if charges is None else float(np.dot(charges, n))
        log_activity = np.log(np.clip(x * gamma, opts.min_mole_fraction, None))
        reaction_residual = np.asarray(
            [
                sum(float(nu) * log_activity[labels.index(label)] for label, nu in reaction.stoichiometry.items())
                - reaction.log_equilibrium_constant
                for reaction in reaction_defs
            ],
            dtype=float,
        )
        residual = np.concatenate([mass_residual, np.asarray([charge_residual]), reaction_residual])
        residual_norm = float(np.max(np.abs(residual))) if residual.size else 0.0
        return {
            "n": n,
            "x": x,
            "activity": activity,
            "mass_residual": mass_residual,
            "charge_residual": charge_residual,
            "reaction_residual": reaction_residual,
            "residual": residual,
            "residual_norm": residual_norm,
        }

    for iteration in range(0, opts.max_iterations + 1):
        payload = evaluate(log_n)
        best_payload = payload
        history.append({"iteration": iteration, "residual_norm": payload["residual_norm"]})
        if payload["residual_norm"] <= opts.tolerance:
            return _result_from_payload(
                labels,
                balance_names,
                payload,
                state_failure_count,
                "converged",
                True,
                history,
                opts,
            )
        if iteration == opts.max_iterations:
            break
        jac = _finite_difference_jacobian(evaluate, log_n, payload["residual"], opts.finite_difference_step)
        step, *_ = np.linalg.lstsq(jac, -payload["residual"], rcond=None)
        log_n = log_n + opts.damping * step

    assert best_payload is not None
    diagnostics = _diagnostics(history, state_failure_count, opts)
    diagnostics["success"] = False
    diagnostics["message"] = "reactive speciation did not converge"
    diagnostics["final_residual_norm"] = best_payload["residual_norm"]
    raise SolutionError("reactive speciation did not converge", diagnostics)


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
    return options


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


def _state_for(
    mixture_factory: Any, x: np.ndarray, species: list[str], T: float, P: float, phase: str
) -> tuple[Any, np.ndarray | None]:
    obj = mixture_factory(x, T, P)
    if hasattr(obj, "state"):
        state = obj.state(T=T, P=P, x=x, phase=phase)
        charges = np.asarray(obj.parameters.get("z", []), dtype=float).flatten()
        if charges.size == 0:
            charges = None
        elif charges.size != len(species):
            raise InputError("mixture_factory returned a mixture with a charge vector that does not match species.")
        return state, charges
    if not hasattr(obj, "activity_coefficient"):
        raise InputError("mixture_factory must return an ePCSAFTMixture or an ePCSAFTState-like object.")
    mix = getattr(obj, "mixture", None)
    charges = None
    if mix is not None:
        raw = np.asarray(mix.parameters.get("z", []), dtype=float).flatten()
        charges = raw if raw.size else None
    return obj, charges


def _finite_difference_jacobian(evaluate: Any, log_n: np.ndarray, residual: np.ndarray, step: float) -> np.ndarray:
    jac = np.zeros((residual.size, log_n.size), dtype=float)
    for idx in range(log_n.size):
        shifted = log_n.copy()
        shifted[idx] += step
        jac[:, idx] = (evaluate(shifted)["residual"] - residual) / step
    return jac


def _result_from_payload(
    species: list[str],
    balance_names: list[str],
    payload: dict[str, Any],
    state_failure_count: int,
    message: str,
    success: bool,
    history: list[dict[str, Any]],
    options: ReactiveSpeciationOptions,
) -> ReactiveSpeciationResult:
    return ReactiveSpeciationResult(
        success=success,
        message=message,
        x={label: float(value) for label, value in zip(species, payload["x"])},
        activity_coefficients={str(k): float(v) for k, v in payload["activity"].items()},
        mass_balance_residuals={name: float(value) for name, value in zip(balance_names, payload["mass_residual"])},
        charge_residual=float(payload["charge_residual"]),
        reaction_residuals=[float(value) for value in payload["reaction_residual"]],
        state_failure_count=state_failure_count,
        diagnostics=_diagnostics(history, state_failure_count, options),
    )


def _diagnostics(
    history: list[dict[str, Any]], state_failure_count: int, options: ReactiveSpeciationOptions
) -> dict[str, Any]:
    return {
        "success": True,
        "iterations": int(history[-1]["iteration"]) if history else 0,
        "history": _json_like(history),
        "state_failure_count": int(state_failure_count),
        "tolerance": float(options.tolerance),
        "phase": str(options.phase),
    }


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
