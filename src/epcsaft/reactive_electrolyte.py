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
    state_failure_count: int
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
            "state_failure_count": int(self.state_failure_count),
            "penalty_residuals": [float(v) for v in self.penalty_residuals],
            "diagnostics": dict(self.diagnostics),
        }


def solve_reactive_electrolyte_bubble(
    *,
    species: Sequence[str],
    mixture_factory: Any,
    T: float,
    P_seed: float,
    balances: Mapping[str, Mapping[str, float]],
    totals: Mapping[str, float],
    reactions: Any,
    initial_x: Any,
    vapor_species: Any,
    volatile_species: Any = None,
    nonvolatile_species: Any = None,
    options: ReactiveElectrolyteBubbleOptions | None = None,
) -> ReactiveElectrolyteBubbleResult:
    """Run native chemical speciation followed by native electrolyte bubble pressure."""
    if options is None:
        options = ReactiveElectrolyteBubbleOptions()
    if not isinstance(options, ReactiveElectrolyteBubbleOptions):
        raise InputError("options must be a ReactiveElectrolyteBubbleOptions instance.")
    speciation_options = options.speciation_options or ReactiveSpeciationOptions()
    bubble_options = options.bubble_options or ElectrolyteBubbleOptions(initial_pressure=float(P_seed))
    if options.bubble_options is None:
        bubble_options = ElectrolyteBubbleOptions(
            initial_pressure=float(P_seed),
            return_best_effort=options.error_mode != "raise",
        )
    chemical = solve_reactive_speciation(
        species=species,
        mixture_factory=mixture_factory,
        T=T,
        P=P_seed,
        balances=balances,
        totals=totals,
        reactions=reactions,
        initial_x=initial_x,
        options=speciation_options,
    )
    x_liq = {label: chemical.x[label] for label in species}
    mixture = mixture_factory([x_liq[label] for label in species], T, P_seed)
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
        "native_entrypoint": "_solve_chemical_equilibrium_native_then__solve_electrolyte_bubble_native",
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
        state_failure_count=int(chemical.state_failure_count) + int(bubble.diagnostics.get("state_failure_count", 0)),
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
    continuation: str = "auto",
) -> list[ReactiveElectrolyteBubbleResult]:
    """Run a native-backed reactive electrolyte bubble sweep with pressure/y continuation."""
    continuation_mode = _normalize_continuation(continuation)
    if options is None:
        options = ReactiveElectrolyteBubbleOptions()
    if not isinstance(options, ReactiveElectrolyteBubbleOptions):
        raise InputError("options must be a ReactiveElectrolyteBubbleOptions instance.")
    results: list[ReactiveElectrolyteBubbleResult] = []
    last_pressure = None
    last_y = None
    for point in points:
        if "T" not in point or "totals" not in point:
            raise InputError("Each reactive electrolyte bubble sweep point requires T and totals.")
        point_options = point.get("options", options)
        if not isinstance(point_options, ReactiveElectrolyteBubbleOptions):
            raise InputError("Each point options entry must be a ReactiveElectrolyteBubbleOptions instance.")
        pressure_seed = float(point.get("P_seed", point.get("P", 101325.0)))
        if continuation_mode == "auto" and (last_pressure is not None or last_y is not None):
            point_options = _with_bubble_continuation(
                point_options,
                initial_pressure=float(last_pressure if last_pressure is not None else pressure_seed),
                initial_y_vap=last_y,
                fallback_pressure=pressure_seed,
            )
            pressure_seed = float(last_pressure if last_pressure is not None else pressure_seed)
        result = solve_reactive_electrolyte_bubble(
            species=species,
            mixture_factory=mixture_factory,
            T=float(point["T"]),
            P_seed=pressure_seed,
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
        if result.success:
            last_pressure = result.P_total
            last_y = result.y_vap
    return results


def _normalize_continuation(value: Any) -> str:
    if isinstance(value, bool):
        return "auto" if value else "none"
    token = str(value).strip().lower()
    if token in {"auto", "on", "true", "1"}:
        return "auto"
    if token in {"none", "off", "false", "0", "disabled"}:
        return "none"
    raise InputError("continuation must be 'auto' or 'none'.")


def _with_bubble_continuation(
    options: ReactiveElectrolyteBubbleOptions,
    *,
    initial_pressure: float,
    initial_y_vap: Mapping[str, float] | None,
    fallback_pressure: float,
) -> ReactiveElectrolyteBubbleOptions:
    base = options.bubble_options or ElectrolyteBubbleOptions(
        initial_pressure=float(fallback_pressure),
        return_best_effort=options.error_mode != "raise",
    )
    return ReactiveElectrolyteBubbleOptions(
        speciation_options=options.speciation_options,
        bubble_options=ElectrolyteBubbleOptions(
            initial_pressure=float(initial_pressure),
            min_pressure=base.min_pressure,
            max_pressure=base.max_pressure,
            max_iterations=base.max_iterations,
            max_vapor_iterations=base.max_vapor_iterations,
            max_bracket_expansions=base.max_bracket_expansions,
            tolerance=base.tolerance,
            vapor_tolerance=base.vapor_tolerance,
            pressure_factor=base.pressure_factor,
            min_composition=base.min_composition,
            charge_tolerance=base.charge_tolerance,
            return_best_effort=base.return_best_effort,
            initial_y_vap=initial_y_vap or base.initial_y_vap,
        ),
        error_mode=options.error_mode,
        penalty_value=options.penalty_value,
        phase_handoff_mass_tolerance=options.phase_handoff_mass_tolerance,
        phase_handoff_charge_tolerance=options.phase_handoff_charge_tolerance,
        phase_handoff_reaction_tolerance=options.phase_handoff_reaction_tolerance,
    )


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
    state_failure_count = int(diagnostics.get("state_failure_count", result.state_failure_count))
    residuals_finite = all(math.isfinite(value) for value in (mass_norm, charge_norm, reaction_norm))
    residuals_within_tolerance = (
        mass_norm <= mass_tolerance and charge_norm <= charge_tolerance and reaction_norm <= reaction_tolerance
    )
    if result.success:
        reason = "strict_success"
    elif state_failure_count > 0:
        reason = "state_failures"
    elif not residuals_finite:
        reason = "nonfinite_residuals"
    elif residuals_within_tolerance:
        reason = "residuals_within_phase_handoff_tolerances"
    else:
        reason = "residuals_exceed_phase_handoff_tolerances"
    success = bool(result.success or (state_failure_count == 0 and residuals_finite and residuals_within_tolerance))
    return {
        "success": success,
        "reason": reason,
        "native_success": bool(diagnostics.get("native_success", result.success)),
        "state_failure_count": state_failure_count,
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
