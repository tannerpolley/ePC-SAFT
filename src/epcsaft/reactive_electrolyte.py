"""Sequential reactive-electrolyte bubble workflow contracts."""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import Any

from ._types import InputError, SolutionError
from .electrolyte_bubble import ElectrolyteBubbleOptions, electrolyte_bubble_pressure
from .reactive_speciation import ReactiveSpeciationOptions, ReactiveSpeciationResult, solve_reactive_speciation

_PHASE_HANDOFF_MASS_TOLERANCE = 1.0e-8
_PHASE_HANDOFF_CHARGE_TOLERANCE = 1.0e-8
_PHASE_HANDOFF_REACTION_TOLERANCE = 1.0e-5
_SWEEP_POINT_KEYS = frozenset({"T", "P", "totals", "initial_x", "options"})


@dataclass(frozen=True, slots=True)
class ReactiveElectrolyteBubbleOptions:
    """Controls for native reactive electrolyte bubble-pressure calculations."""

    speciation_options: ReactiveSpeciationOptions | None = None
    bubble_options: ElectrolyteBubbleOptions | None = None
    phase_handoff_mass_tolerance: float = _PHASE_HANDOFF_MASS_TOLERANCE
    phase_handoff_charge_tolerance: float = _PHASE_HANDOFF_CHARGE_TOLERANCE
    phase_handoff_reaction_tolerance: float = _PHASE_HANDOFF_REACTION_TOLERANCE
    error_mode: str = "raise"
    penalty_value: float = 1.0e6


@dataclass(frozen=True, slots=True)
class ReactiveElectrolyteBubbleResult:
    """Structured result shape reserved for native reactive electrolyte bubble calculations."""

    success: bool
    message: str
    x_liq: Mapping[str, float]
    activity_coefficients: Mapping[str, float]
    mass_balance_residuals: Mapping[str, float]
    charge_residual: float
    reaction_residuals: Sequence[float]
    named_reaction_residuals: Mapping[str, float]
    P_total: float
    y_vap: Mapping[str, float]
    partial_pressures: Mapping[str, float]
    fugacity_residual: Mapping[str, float]
    fugacity_residual_norm: float
    penalty_residuals: Sequence[float]
    diagnostics: Mapping[str, Any]

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-like result payload."""
        return {
            "success": bool(self.success),
            "message": str(self.message),
            "x_liq": {str(k): float(v) for k, v in self.x_liq.items()},
            "activity_coefficients": {str(k): float(v) for k, v in self.activity_coefficients.items()},
            "mass_balance_residuals": {str(k): float(v) for k, v in self.mass_balance_residuals.items()},
            "charge_residual": float(self.charge_residual),
            "reaction_residuals": [float(v) for v in self.reaction_residuals],
            "named_reaction_residuals": {str(k): float(v) for k, v in self.named_reaction_residuals.items()},
            "P_total": float(self.P_total),
            "y_vap": {str(k): float(v) for k, v in self.y_vap.items()},
            "partial_pressures": {str(k): float(v) for k, v in self.partial_pressures.items()},
            "fugacity_residual": {str(k): float(v) for k, v in self.fugacity_residual.items()},
            "fugacity_residual_norm": float(self.fugacity_residual_norm),
            "penalty_residuals": [float(v) for v in self.penalty_residuals],
            "diagnostics": dict(self.diagnostics),
        }


def solve_reactive_electrolyte_bubble(
    *,
    species: Sequence[str],
    mixture_factory: Any,
    T: float,
    P: float,
    balances: Mapping[str, Mapping[str, float]],
    totals: Mapping[str, float],
    reactions: Any,
    initial_x: Any,
    vapor_species: Any,
    volatile_species: Any = None,
    nonvolatile_species: Any = None,
    options: ReactiveElectrolyteBubbleOptions | None = None,
) -> ReactiveElectrolyteBubbleResult:
    """Run native chemical speciation followed by native Ipopt electrolyte bubble pressure."""
    if options is None:
        options = ReactiveElectrolyteBubbleOptions()
    if not isinstance(options, ReactiveElectrolyteBubbleOptions):
        raise InputError("options must be a ReactiveElectrolyteBubbleOptions instance.")
    speciation_options = options.speciation_options or ReactiveSpeciationOptions()
    bubble_options = options.bubble_options or ElectrolyteBubbleOptions()
    chemical = solve_reactive_speciation(
        species=species,
        mixture_factory=mixture_factory,
        T=T,
        P=P,
        balances=balances,
        totals=totals,
        reactions=reactions,
        initial_x=initial_x,
        options=speciation_options,
    )
    x_liq = {label: chemical.x[label] for label in species}
    mixture = mixture_factory([x_liq[label] for label in species], T, P)
    bubble_failure: SolutionError | None = None
    try:
        bubble = electrolyte_bubble_pressure(
            mixture,
            T=T,
            x_liq=[x_liq[label] for label in species],
            vapor_species=vapor_species if vapor_species is not None else volatile_species,
            volatile_species=volatile_species,
            nonvolatile_species=nonvolatile_species,
            options=bubble_options,
        )
    except SolutionError as exc:
        if options.error_mode == "raise":
            raise
        bubble_failure = exc
        bubble = _failed_bubble_result(
            message=str(getattr(exc, "message", str(exc))),
            diagnostics=dict(getattr(exc, "diagnostics", {}) or {}),
            species=species,
            x_liq=x_liq,
            vapor_species=vapor_species if vapor_species is not None else volatile_species,
        )
    handoff = _speciation_phase_handoff_diagnostics(chemical, options)
    diagnostics = {
        "speciation": chemical.to_dict(),
        "bubble": bubble.to_dict(),
        "speciation_strict_success": bool(chemical.success),
        "speciation_phase_handoff_success": handoff["success"],
        "speciation_phase_handoff": handoff,
        "bubble_success": bool(bubble.success),
        "native_entrypoint": "native_chemical_speciation_then_ipopt_electrolyte_bubble",
        "solver_language": "c++",
    }
    if bubble_failure is not None:
        diagnostics["bubble_failure"] = {
            "type": type(bubble_failure).__name__,
            "message": str(getattr(bubble_failure, "message", str(bubble_failure))),
            "diagnostics": dict(getattr(bubble_failure, "diagnostics", {}) or {}),
        }
    success = bool(diagnostics["speciation_phase_handoff_success"] and bubble.success)
    message = "converged" if success else _reactive_bubble_failure_message(chemical, bubble, diagnostics)
    return ReactiveElectrolyteBubbleResult(
        success=success,
        message=message,
        x_liq=x_liq,
        activity_coefficients=chemical.activity_coefficients,
        mass_balance_residuals=chemical.mass_balance_residuals,
        charge_residual=chemical.charge_residual,
        reaction_residuals=chemical.reaction_residuals,
        named_reaction_residuals=chemical.named_reaction_residuals,
        P_total=bubble.P,
        y_vap=bubble.y_vap,
        partial_pressures=bubble.partial_pressures,
        fugacity_residual=bubble.fugacity_residual,
        fugacity_residual_norm=bubble.fugacity_residual_norm,
        penalty_residuals=[],
        diagnostics=diagnostics,
    )


def solve_reactive_electrolyte_bubble_sweep(
    *,
    species: Sequence[str],
    mixture_factory: Any,
    points: Sequence[Mapping[str, Any]],
    balances: Mapping[str, Mapping[str, float]],
    reactions: Any,
    vapor_species: Any,
    volatile_species: Any = None,
    nonvolatile_species: Any = None,
    options: ReactiveElectrolyteBubbleOptions | None = None,
) -> list[ReactiveElectrolyteBubbleResult]:
    """Apply the reactive electrolyte bubble route contract across sweep points."""
    if options is None:
        options = ReactiveElectrolyteBubbleOptions()
    if not isinstance(options, ReactiveElectrolyteBubbleOptions):
        raise InputError("options must be a ReactiveElectrolyteBubbleOptions instance.")
    results: list[ReactiveElectrolyteBubbleResult] = []
    for point in points:
        unknown_keys = sorted(set(point) - _SWEEP_POINT_KEYS)
        if unknown_keys:
            raise InputError(f"Unsupported reactive electrolyte bubble sweep point key(s): {unknown_keys}.")
        if "T" not in point or "totals" not in point:
            raise InputError("Each reactive electrolyte bubble sweep point requires T and totals.")
        point_options = point.get("options", options)
        if not isinstance(point_options, ReactiveElectrolyteBubbleOptions):
            raise InputError("Each point options entry must be a ReactiveElectrolyteBubbleOptions instance.")
        pressure = float(point.get("P", 101325.0))
        result = solve_reactive_electrolyte_bubble(
            species=species,
            mixture_factory=mixture_factory,
            T=float(point["T"]),
            P=pressure,
            balances=balances,
            totals=point["totals"],
            reactions=reactions,
            initial_x=point.get("initial_x"),
            vapor_species=vapor_species,
            volatile_species=volatile_species,
            nonvolatile_species=nonvolatile_species,
            options=point_options,
        )
        results.append(result)
    return results


def _speciation_phase_handoff_diagnostics(
    result: ReactiveSpeciationResult,
    options: ReactiveElectrolyteBubbleOptions,
) -> dict[str, Any]:
    """Return residual-family diagnostics for speciation-to-phase handoff."""
    diagnostics = result.diagnostics
    mass_tolerance = _positive_tolerance(options.phase_handoff_mass_tolerance, "phase_handoff_mass_tolerance")
    charge_tolerance = _positive_tolerance(options.phase_handoff_charge_tolerance, "phase_handoff_charge_tolerance")
    reaction_tolerance = _positive_tolerance(
        options.phase_handoff_reaction_tolerance,
        "phase_handoff_reaction_tolerance",
    )
    mass_norm = float(diagnostics.get("mass_residual_norm", float("inf")))
    charge_norm = float(diagnostics.get("charge_residual_abs", abs(result.charge_residual)))
    reaction_norm = float(
        diagnostics.get(
            "reaction_residual_norm",
            max((abs(value) for value in result.reaction_residuals), default=0.0),
        )
    )
    residuals_finite = all(math.isfinite(value) for value in (mass_norm, charge_norm, reaction_norm))
    residuals_within_tolerance = (
        mass_norm <= mass_tolerance and charge_norm <= charge_tolerance and reaction_norm <= reaction_tolerance
    )
    if result.success:
        reason = "strict_success"
    elif not residuals_finite:
        reason = "nonfinite_residuals"
    elif residuals_within_tolerance:
        reason = "residuals_within_phase_handoff_tolerances"
    else:
        reason = "residuals_exceed_phase_handoff_tolerances"
    success = bool(result.success or (residuals_finite and residuals_within_tolerance))
    return {
        "success": success,
        "reason": reason,
        "native_success": bool(diagnostics.get("native_success", result.success)),
        "mass_residual_norm": mass_norm,
        "charge_residual_abs": charge_norm,
        "reaction_residual_norm": reaction_norm,
        "mass_tolerance": mass_tolerance,
        "charge_tolerance": charge_tolerance,
        "reaction_tolerance": reaction_tolerance,
    }


def _positive_tolerance(value: float, name: str) -> float:
    tolerance = float(value)
    if not math.isfinite(tolerance) or tolerance <= 0.0:
        raise InputError(f"{name} must be a finite positive value.")
    return tolerance


def _failed_bubble_result(
    *,
    message: str,
    diagnostics: Mapping[str, Any],
    species: Sequence[str],
    x_liq: Mapping[str, float],
    vapor_species: Any,
) -> Any:
    """Build a fixed-shape bubble result from a strict bubble failure."""
    from .electrolyte_bubble import ElectrolyteBubbleResult

    labels = _normalize_vapor_labels(vapor_species)
    pressure = float(diagnostics.get("best_P", diagnostics.get("P", 0.0)) or 0.0)
    y_vap = _mapping_for_labels(labels, diagnostics.get("best_y_vap"), default=0.0)
    partial_values = diagnostics.get("best_partial_pressures")
    partial = _mapping_for_labels(labels, partial_values, default=0.0)
    residual = _mapping_for_labels(labels, diagnostics.get("fugacity_residual"), default=float("nan"))
    if partial_values is None and labels:
        partial = {label: float(y_vap.get(label, 0.0)) * pressure for label in labels}
    return ElectrolyteBubbleResult(
        success=False,
        message=message,
        P=pressure,
        y_vap=y_vap,
        x_liq=[float(x_liq[label]) for label in species],
        ln_phi_liq={label: float("nan") for label in labels},
        ln_phi_vap={label: float("nan") for label in labels},
        fugacity_residual=residual,
        fugacity_residual_norm=float(diagnostics.get("best_fugacity_residual_norm", float("nan"))),
        charge_residual=float(diagnostics.get("charge_residual", float("nan"))),
        partial_pressures=partial,
        diagnostics=dict(diagnostics),
    )


def _normalize_vapor_labels(values: Any) -> list[str]:
    if values is None:
        return []
    if isinstance(values, str):
        return [values]
    return [str(value) for value in values]


def _mapping_for_labels(labels: Sequence[str], values: Any, *, default: float) -> dict[str, float]:
    if isinstance(values, Mapping):
        return {label: float(values.get(label, default)) for label in labels}
    if values is None:
        return {label: float(default) for label in labels}
    try:
        items = list(values)
    except TypeError:
        return {label: float(default) for label in labels}
    return {label: float(items[index]) if index < len(items) else float(default) for index, label in enumerate(labels)}


def _reactive_bubble_failure_message(
    chemical: ReactiveSpeciationResult,
    bubble: Any,
    diagnostics: Mapping[str, Any],
) -> str:
    if not bool(diagnostics["speciation_phase_handoff_success"]):
        if chemical.success:
            return "reactive electrolyte speciation failed phase-handoff checks"
        return "reactive electrolyte speciation did not meet phase-handoff tolerances"
    return str(bubble.message)
