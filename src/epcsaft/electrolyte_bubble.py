"""Electrolyte bubble-pressure helpers built on public state evaluations."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from typing import Any

import numpy as np

from ._types import InputError, SolutionError


@dataclass(frozen=True, slots=True)
class ElectrolyteBubbleOptions:
    """Numerical controls for electrolyte bubble-pressure solves."""

    initial_pressure: float = 1.0e5
    min_pressure: float = 1.0
    max_pressure: float = 1.0e8
    max_iterations: int = 80
    max_vapor_iterations: int = 30
    max_bracket_expansions: int = 40
    tolerance: float = 1.0e-6
    vapor_tolerance: float = 1.0e-10
    pressure_factor: float = 2.0
    min_composition: float = 1.0e-14
    charge_tolerance: float = 1.0e-8
    return_best_effort: bool = False
    initial_y_vap: Mapping[str, float] | None = None


@dataclass(frozen=True, slots=True)
class ElectrolyteBubbleResult:
    """Structured result returned by electrolyte bubble-pressure calculations."""

    success: bool
    message: str
    P: float
    y_vap: dict[str, float]
    x_liq: np.ndarray
    ln_phi_liq: dict[str, float]
    ln_phi_vap: dict[str, float]
    fugacity_residual: dict[str, float]
    fugacity_residual_norm: float
    charge_residual: float
    partial_pressures: dict[str, float]
    diagnostics: dict[str, Any]

    def __post_init__(self) -> None:
        object.__setattr__(self, "success", bool(self.success))
        object.__setattr__(self, "message", str(self.message))
        object.__setattr__(self, "P", float(self.P))
        object.__setattr__(self, "y_vap", {str(k): float(v) for k, v in self.y_vap.items()})
        object.__setattr__(self, "x_liq", np.asarray(self.x_liq, dtype=float))
        object.__setattr__(self, "ln_phi_liq", {str(k): float(v) for k, v in self.ln_phi_liq.items()})
        object.__setattr__(self, "ln_phi_vap", {str(k): float(v) for k, v in self.ln_phi_vap.items()})
        object.__setattr__(self, "fugacity_residual", {str(k): float(v) for k, v in self.fugacity_residual.items()})
        object.__setattr__(self, "fugacity_residual_norm", float(self.fugacity_residual_norm))
        object.__setattr__(self, "charge_residual", float(self.charge_residual))
        object.__setattr__(self, "partial_pressures", {str(k): float(v) for k, v in self.partial_pressures.items()})
        object.__setattr__(self, "diagnostics", dict(self.diagnostics))

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-like result payload."""
        return {
            "success": self.success,
            "message": self.message,
            "P": self.P,
            "y_vap": dict(self.y_vap),
            "x_liq": self.x_liq.tolist(),
            "ln_phi_liq": dict(self.ln_phi_liq),
            "ln_phi_vap": dict(self.ln_phi_vap),
            "fugacity_residual": dict(self.fugacity_residual),
            "fugacity_residual_norm": self.fugacity_residual_norm,
            "charge_residual": self.charge_residual,
            "partial_pressures": dict(self.partial_pressures),
            "diagnostics": _json_like(self.diagnostics),
        }


def electrolyte_bubble_pressure(
    mixture: Any,
    *,
    T: float,
    x_liq: Any = None,
    z: Any = None,
    volatile_species: Any = None,
    vapor_species: Any = None,
    nonvolatile_species: Any = None,
    options: ElectrolyteBubbleOptions | None = None,
) -> ElectrolyteBubbleResult:
    """Solve electrolyte bubble pressure with neutral vapor species and liquid-only ions."""
    opts = _normalize_options(options)
    temperature = _positive_scalar(T, "T")
    composition = _normalize_liquid_composition(x_liq if x_liq is not None else z, mixture.ncomp, opts)
    species = list(mixture.species)
    charges = _mixture_charges(mixture)
    if not np.any(np.abs(charges) > 1.0e-12):
        raise InputError("electrolyte_bubble_pressure requires an ion-containing mixture.")
    charge_residual = float(np.dot(composition, charges))
    if abs(charge_residual) > opts.charge_tolerance:
        raise InputError("electrolyte_bubble_pressure x_liq must be charge neutral.")

    vapor_labels = _normalize_species_subset(species, vapor_species if vapor_species is not None else volatile_species)
    if not vapor_labels:
        raise InputError("electrolyte_bubble_pressure requires at least one vapor species.")
    volatile_labels = _normalize_species_subset(
        species, volatile_species if volatile_species is not None else vapor_labels
    )
    vapor_idx = [species.index(label) for label in vapor_labels]
    volatile_idx = [species.index(label) for label in volatile_labels]
    if any(abs(float(charges[idx])) > 1.0e-12 for idx in vapor_idx):
        raise InputError("electrolyte_bubble_pressure vapor_species must be neutral vapor species.")
    if any(abs(float(charges[idx])) > 1.0e-12 for idx in volatile_idx):
        raise InputError("electrolyte_bubble_pressure volatile_species must be neutral species.")

    nonvolatile_labels = (
        _normalize_species_subset(species, nonvolatile_species)
        if nonvolatile_species is not None
        else [label for label, charge in zip(species, charges) if abs(float(charge)) > 1.0e-12]
    )
    overlap = set(vapor_labels).intersection(nonvolatile_labels)
    if overlap:
        raise InputError(
            "electrolyte_bubble_pressure vapor_species cannot also be nonvolatile: {}".format(sorted(overlap))
        )

    vapor_mixture = _build_neutral_submixture(mixture, vapor_idx)
    initial_y_seed = _normalize_vapor_seed(vapor_labels, opts.initial_y_vap, opts)
    history: list[dict[str, Any]] = []
    state_failure_count = 0
    best_point: dict[str, Any] | None = None

    def evaluate(pressure: float, y_seed: np.ndarray | None = None) -> dict[str, Any]:
        nonlocal best_point, state_failure_count
        try:
            out = _bubble_objective(
                mixture,
                vapor_mixture,
                T=temperature,
                P=pressure,
                x_liq=composition,
                species=species,
                vapor_species=vapor_labels,
                vapor_indices=vapor_idx,
                options=opts,
                y_seed=y_seed,
            )
        except Exception as exc:
            state_failure_count += 1
            raise exc
        objective = float(out["objective"])
        residual = _fugacity_residual(out, vapor_labels)
        residual_norm = float(max((abs(value) for value in residual.values()), default=0.0))
        y_vap = {label: float(value) for label, value in zip(vapor_labels, out["y_vap"])}
        partial_pressures = {label: float(pressure * y) for label, y in y_vap.items()}
        history.append(
            {
                "P": float(pressure),
                "logP": float(np.log(pressure)),
                "objective": objective,
                "y_vap": y_vap,
            }
        )
        if np.isfinite(objective) and (best_point is None or abs(objective) < abs(float(best_point["objective"]))):
            best_point = {
                "P": float(pressure),
                "objective": objective,
                "y_vap": y_vap,
                "partial_pressures": partial_pressures,
                "fugacity_residual_norm": residual_norm,
                "payload": out,
            }
        return out

    try:
        bracket = _find_pressure_bracket(evaluate, opts, initial_y_seed)
    except SolutionError as exc:
        diagnostics = dict(exc.diagnostics or {})
        message = exc.message
        diagnostics.update(
            {
                "success": False,
                "message": message,
                "state_failure_count": int(state_failure_count),
                "vapor_history": _bounded_history(history),
                "tolerance": float(opts.tolerance),
            }
        )
        diagnostics.update(_best_diagnostics(best_point))
        if opts.return_best_effort and best_point is not None:
            return _make_result(
                success=False,
                message=message,
                pressure=float(best_point["P"]),
                payload=best_point["payload"],
                species=species,
                vapor_species=vapor_labels,
                composition=composition,
                charge_residual=charge_residual,
                diagnostics=_json_like(diagnostics),
            )
        raise SolutionError(message, _json_like(diagnostics)) from exc
    low_p, low_eval, high_p, high_eval = bracket
    best_p = low_p
    best_eval = low_eval
    y_seed = initial_y_seed if initial_y_seed is not None else low_eval["y_vap"]
    low_log = float(np.log(low_p))
    high_log = float(np.log(high_p))
    for iteration in range(1, opts.max_iterations + 1):
        if abs(low_eval["objective"]) <= abs(high_eval["objective"]):
            best_p, best_eval = low_p, low_eval
        else:
            best_p, best_eval = high_p, high_eval
        mid_log = 0.5 * (low_log + high_log)
        mid_p = float(np.exp(mid_log))
        mid_eval = evaluate(mid_p, y_seed=y_seed)
        y_seed = mid_eval["y_vap"]
        if abs(mid_eval["objective"]) <= opts.tolerance:
            best_p, best_eval = mid_p, mid_eval
            break
        if np.sign(low_eval["objective"]) == np.sign(mid_eval["objective"]):
            low_p, low_eval = mid_p, mid_eval
            low_log = mid_log
        else:
            high_p, high_eval = mid_p, mid_eval
            high_log = mid_log
    else:
        message = "electrolyte bubble pressure did not converge"
        diagnostics = {
            "success": False,
            "message": message,
            "iterations": int(opts.max_iterations),
            "state_failure_count": int(state_failure_count),
            "pressure_bracket": [float(low_p), float(high_p)],
            "log_pressure_bracket": [float(low_log), float(high_log)],
            "vapor_history": _bounded_history(history),
            "tolerance": float(opts.tolerance),
        }
        diagnostics.update(_best_diagnostics(best_point))
        if opts.return_best_effort and best_point is not None:
            return _make_result(
                success=False,
                message=message,
                pressure=float(best_point["P"]),
                payload=best_point["payload"],
                species=species,
                vapor_species=vapor_labels,
                composition=composition,
                charge_residual=charge_residual,
                diagnostics=_json_like(diagnostics),
            )
        raise SolutionError(message, _json_like(diagnostics))

    diagnostics = {
        "success": True,
        "iterations": int(len(history)),
        "state_failure_count": int(state_failure_count),
        "pressure_bracket": [float(low_p), float(high_p)],
        "log_pressure_bracket": [float(low_log), float(high_log)],
        "vapor_history": _bounded_history(history),
        "vapor_species": list(vapor_labels),
        "volatile_species": list(volatile_labels),
        "nonvolatile_species": list(nonvolatile_labels),
        "tolerance": float(opts.tolerance),
        "log_pressure_solve": True,
        "used_initial_y_vap": initial_y_seed is not None,
        "message": "converged",
    }
    diagnostics.update(_best_diagnostics(best_point))
    return _make_result(
        success=True,
        message="converged",
        pressure=best_p,
        payload=best_eval,
        species=species,
        vapor_species=vapor_labels,
        composition=composition,
        charge_residual=charge_residual,
        diagnostics=diagnostics,
    )


def _normalize_options(options: ElectrolyteBubbleOptions | None) -> ElectrolyteBubbleOptions:
    if options is None:
        return ElectrolyteBubbleOptions()
    if not isinstance(options, ElectrolyteBubbleOptions):
        raise InputError("options must be an ElectrolyteBubbleOptions instance.")
    if options.max_iterations <= 0 or options.max_vapor_iterations <= 0:
        raise InputError("ElectrolyteBubbleOptions iteration counts must be positive.")
    if options.tolerance <= 0.0 or options.vapor_tolerance <= 0.0:
        raise InputError("ElectrolyteBubbleOptions tolerances must be positive.")
    if options.min_pressure <= 0.0 or options.max_pressure <= options.min_pressure:
        raise InputError("ElectrolyteBubbleOptions pressure bounds are invalid.")
    if options.initial_pressure <= 0.0:
        raise InputError("ElectrolyteBubbleOptions initial_pressure must be positive.")
    if options.pressure_factor <= 1.0:
        raise InputError("ElectrolyteBubbleOptions pressure_factor must be greater than 1.")
    if not isinstance(options.return_best_effort, bool):
        raise InputError("ElectrolyteBubbleOptions.return_best_effort must be a bool.")
    return options


def _positive_scalar(value: Any, name: str) -> float:
    try:
        out = float(value)
    except (TypeError, ValueError) as exc:
        raise InputError(f"{name} must be a positive finite number.") from exc
    if not np.isfinite(out) or out <= 0.0:
        raise InputError(f"{name} must be a positive finite number.")
    return out


def _normalize_liquid_composition(value: Any, ncomp: int, options: ElectrolyteBubbleOptions) -> np.ndarray:
    if value is None:
        raise InputError("electrolyte_bubble_pressure requires x_liq.")
    array = np.asarray(value, dtype=float).flatten()
    if array.size != int(ncomp):
        raise InputError("x_liq length must match mixture component count.")
    if np.any(~np.isfinite(array)):
        raise InputError("x_liq values must be finite.")
    if np.any(array < -options.min_composition):
        raise InputError("x_liq values must be non-negative.")
    array = np.clip(array, 0.0, None)
    total = float(np.sum(array))
    if total <= 0.0:
        raise InputError("x_liq must have a positive sum.")
    return array / total


def _mixture_charges(mixture: Any) -> np.ndarray:
    charges = np.asarray(mixture.parameters.get("z", []), dtype=float).flatten()
    if charges.size == 0:
        return np.zeros(int(mixture.ncomp), dtype=float)
    if charges.size != int(mixture.ncomp):
        raise InputError("mixture charge vector length must match component count.")
    return charges


def _normalize_species_subset(species: list[str], labels: Any) -> list[str]:
    if labels is None:
        return []
    if isinstance(labels, str):
        labels = [labels]
    out: list[str] = []
    for raw in labels:
        label = str(raw)
        if label not in species:
            raise InputError(f"Unknown species '{label}'. Available species={species}.")
        if label not in out:
            out.append(label)
    return out


def _normalize_vapor_seed(
    vapor_species: list[str], seed: Mapping[str, float] | None, options: ElectrolyteBubbleOptions
) -> np.ndarray | None:
    if seed is None:
        return None
    if not isinstance(seed, Mapping):
        raise InputError(
            "ElectrolyteBubbleOptions.initial_y_vap must be a mapping from vapor species to mole fraction."
        )
    unknown = sorted(set(str(label) for label in seed) - set(vapor_species))
    if unknown:
        raise InputError(f"initial_y_vap contains unknown vapor species: {unknown}.")
    missing = [label for label in vapor_species if label not in seed]
    if missing:
        raise InputError(f"initial_y_vap is missing vapor species: {missing}.")
    y = np.asarray([float(seed[label]) for label in vapor_species], dtype=float)
    if np.any(~np.isfinite(y)) or np.any(y < -options.min_composition):
        raise InputError("initial_y_vap values must be finite and non-negative.")
    y = np.clip(y, 0.0, None)
    total = float(np.sum(y))
    if total <= 0.0:
        raise InputError("initial_y_vap must have a positive sum.")
    return y / total


def _build_neutral_submixture(mixture: Any, indices: list[int]) -> Any:
    from .epcsaft import ePCSAFTMixture

    params = mixture.parameters
    ncomp = int(mixture.ncomp)
    out: dict[str, Any] = {}
    for key, value in params.items():
        if key in {"z", "dielc", "d_born", "f_solv", "elec_model"}:
            continue
        if isinstance(value, np.ndarray):
            if value.ndim == 1 and value.size == ncomp:
                out[key] = value[indices].copy()
            elif value.ndim == 2 and value.shape == (ncomp, ncomp):
                out[key] = value[np.ix_(indices, indices)].copy()
            else:
                out[key] = value.copy()
        elif isinstance(value, list) and len(value) == ncomp:
            out[key] = [value[i] for i in indices]
        else:
            out[key] = value
    return ePCSAFTMixture.from_params(out, species=[mixture.species[i] for i in indices])


def _bubble_objective(
    mixture: Any,
    vapor_mixture: Any,
    *,
    T: float,
    P: float,
    x_liq: np.ndarray,
    species: list[str],
    vapor_species: list[str],
    vapor_indices: list[int],
    options: ElectrolyteBubbleOptions,
    y_seed: np.ndarray | None,
) -> dict[str, Any]:
    liquid = mixture.state(T=T, P=P, x=x_liq, phase="liq")
    ln_phi_liq_full = np.asarray(liquid.fugacity_coefficient(), dtype=float)
    x_v = np.clip(x_liq[vapor_indices], options.min_composition, None)
    y = x_v / float(np.sum(x_v)) if y_seed is None else np.asarray(y_seed, dtype=float).flatten()
    y = np.clip(y, options.min_composition, None)
    y = y / float(np.sum(y))
    ln_phi_vap = None
    weights = None
    for _ in range(options.max_vapor_iterations):
        vapor = vapor_mixture.state(T=T, P=P, x=y, phase="vap")
        ln_phi_vap = np.asarray(vapor.fugacity_coefficient(), dtype=float)
        weights = x_v * np.exp(np.asarray(ln_phi_liq_full[vapor_indices], dtype=float) - ln_phi_vap)
        total = float(np.sum(weights))
        if not np.isfinite(total) or total <= 0.0:
            raise SolutionError("electrolyte bubble pressure produced invalid vapor weights")
        y_next = weights / total
        if float(np.max(np.abs(y_next - y))) <= options.vapor_tolerance:
            y = y_next
            break
        y = y_next
    assert ln_phi_vap is not None and weights is not None
    return {
        "objective": float(np.sum(weights) - 1.0),
        "y_vap": y,
        "ln_phi_liq_full": ln_phi_liq_full,
        "ln_phi_liq_vap": np.asarray(ln_phi_liq_full[vapor_indices], dtype=float),
        "ln_phi_vap": ln_phi_vap,
        "x_v": x_v,
        "vapor_species": list(vapor_species),
    }


def _find_pressure_bracket(
    evaluate: Any, options: ElectrolyteBubbleOptions, initial_y_seed: np.ndarray | None = None
) -> tuple[float, dict[str, Any], float, dict[str, Any]]:
    center = min(max(float(options.initial_pressure), float(options.min_pressure)), float(options.max_pressure))
    center_eval = evaluate(center, y_seed=initial_y_seed)
    if abs(center_eval["objective"]) <= options.tolerance:
        return center, center_eval, center, center_eval
    candidates: list[tuple[float, dict[str, Any]]] = [(center, center_eval)]
    min_log = float(np.log(options.min_pressure))
    max_log = float(np.log(options.max_pressure))
    step_log = float(np.log(options.pressure_factor))
    low_log = float(np.log(center))
    high_log = float(np.log(center))
    for _ in range(options.max_bracket_expansions):
        low_log = max(min_log, low_log - step_log)
        high_log = min(max_log, high_log + step_log)
        low = float(np.exp(low_log))
        high = float(np.exp(high_log))
        for pressure in (low, high):
            if any(abs(pressure - existing[0]) <= 1.0e-12 * max(1.0, pressure) for existing in candidates):
                continue
            try:
                candidates.append((pressure, evaluate(pressure)))
            except Exception:
                continue
        ordered = sorted(candidates, key=lambda item: item[0])
        for left, right in zip(ordered[:-1], ordered[1:]):
            if np.sign(left[1]["objective"]) == 0 or np.sign(right[1]["objective"]) == 0:
                return left[0], left[1], right[0], right[1]
            if np.sign(left[1]["objective"]) != np.sign(right[1]["objective"]):
                return left[0], left[1], right[0], right[1]
        if low <= options.min_pressure and high >= options.max_pressure:
            break
    diagnostics = {
        "success": False,
        "message": "electrolyte bubble pressure could not bracket a pressure root",
        "evaluated_pressures": [
            {"P": p, "logP": float(np.log(p)), "objective": item["objective"]} for p, item in candidates
        ],
        "pressure_bounds": [float(options.min_pressure), float(options.max_pressure)],
        "log_pressure_bounds": [min_log, max_log],
    }
    raise SolutionError("electrolyte bubble pressure could not bracket a pressure root", diagnostics)


def _fugacity_residual(payload: dict[str, Any], vapor_species: list[str]) -> dict[str, float]:
    y = np.asarray(payload["y_vap"], dtype=float)
    ln_phi_vap = np.asarray(payload["ln_phi_vap"], dtype=float)
    ln_phi_liq = np.asarray(payload["ln_phi_liq_vap"], dtype=float)
    x_v = np.asarray(payload["x_v"], dtype=float)
    residual = np.log(y) + ln_phi_vap - np.log(x_v) - ln_phi_liq
    return {label: float(value) for label, value in zip(vapor_species, residual)}


def _make_result(
    *,
    success: bool,
    message: str,
    pressure: float,
    payload: dict[str, Any],
    species: list[str],
    vapor_species: list[str],
    composition: np.ndarray,
    charge_residual: float,
    diagnostics: dict[str, Any],
) -> ElectrolyteBubbleResult:
    residual = _fugacity_residual(payload, vapor_species)
    residual_norm = float(max((abs(value) for value in residual.values()), default=0.0))
    y_vap = {label: float(value) for label, value in zip(vapor_species, payload["y_vap"])}
    return ElectrolyteBubbleResult(
        success=success,
        message=message,
        P=float(pressure),
        y_vap=y_vap,
        x_liq=composition,
        ln_phi_liq={label: float(value) for label, value in zip(species, payload["ln_phi_liq_full"])},
        ln_phi_vap={label: float(value) for label, value in zip(vapor_species, payload["ln_phi_vap"])},
        fugacity_residual=residual,
        fugacity_residual_norm=residual_norm,
        charge_residual=charge_residual,
        partial_pressures={label: float(pressure * y) for label, y in y_vap.items()},
        diagnostics=diagnostics,
    )


def _best_diagnostics(best_point: dict[str, Any] | None) -> dict[str, Any]:
    if best_point is None:
        return {
            "best_P": None,
            "best_objective": None,
            "best_y_vap": {},
            "best_partial_pressures": {},
            "best_fugacity_residual_norm": None,
        }
    return {
        "best_P": float(best_point["P"]),
        "best_objective": float(best_point["objective"]),
        "best_y_vap": dict(best_point["y_vap"]),
        "best_partial_pressures": dict(best_point["partial_pressures"]),
        "best_fugacity_residual_norm": float(best_point["fugacity_residual_norm"]),
    }


def _bounded_history(history: list[dict[str, Any]], limit: int = 50) -> list[dict[str, Any]]:
    return history[-limit:]


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
