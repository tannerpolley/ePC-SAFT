"""Composable reactive electrolyte equilibrium workflows."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, replace
from typing import Any

import numpy as np

from ._types import InputError, SolutionError
from .electrolyte_bubble import ElectrolyteBubbleOptions
from .electrolyte_bubble import ElectrolyteBubbleResult
from .reactive_speciation import ReactiveSpeciationOptions
from .reactive_speciation import ReactiveSpeciationResult
from .reactive_speciation import solve_reactive_speciation


@dataclass(frozen=True, slots=True)
class ReactiveElectrolyteBubbleOptions:
    """Controls for composed reactive-speciation plus electrolyte-bubble workflows."""

    speciation_options: ReactiveSpeciationOptions | None = None
    bubble_options: ElectrolyteBubbleOptions | None = None
    error_mode: str = "raise"
    penalty_value: float = 1.0e6


@dataclass(frozen=True, slots=True)
class ReactiveElectrolyteBubbleResult:
    """Structured result for a composed reactive electrolyte bubble calculation."""

    success: bool
    message: str
    x_liq: dict[str, float]
    activity_coefficients: dict[str, float]
    mass_balance_residuals: dict[str, float]
    charge_residual: float
    reaction_residuals: list[float]
    named_reaction_residuals: dict[str, float]
    P_total: float
    y_vap: dict[str, float]
    partial_pressures: dict[str, float]
    fugacity_residual: dict[str, float]
    fugacity_residual_norm: float
    state_failure_count: int
    penalty_residuals: list[float]
    diagnostics: dict[str, Any]

    def __post_init__(self) -> None:
        object.__setattr__(self, "success", bool(self.success))
        object.__setattr__(self, "message", str(self.message))
        object.__setattr__(self, "x_liq", {str(k): float(v) for k, v in self.x_liq.items()})
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
        object.__setattr__(self, "P_total", float(self.P_total))
        object.__setattr__(self, "y_vap", {str(k): float(v) for k, v in self.y_vap.items()})
        object.__setattr__(self, "partial_pressures", {str(k): float(v) for k, v in self.partial_pressures.items()})
        object.__setattr__(self, "fugacity_residual", {str(k): float(v) for k, v in self.fugacity_residual.items()})
        object.__setattr__(self, "fugacity_residual_norm", float(self.fugacity_residual_norm))
        object.__setattr__(self, "state_failure_count", int(self.state_failure_count))
        object.__setattr__(self, "penalty_residuals", [float(v) for v in self.penalty_residuals])
        object.__setattr__(self, "diagnostics", dict(self.diagnostics))

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-like result payload."""
        return {
            "success": self.success,
            "message": self.message,
            "x_liq": dict(self.x_liq),
            "activity_coefficients": dict(self.activity_coefficients),
            "mass_balance_residuals": dict(self.mass_balance_residuals),
            "charge_residual": self.charge_residual,
            "reaction_residuals": list(self.reaction_residuals),
            "named_reaction_residuals": dict(self.named_reaction_residuals),
            "P_total": self.P_total,
            "y_vap": dict(self.y_vap),
            "partial_pressures": dict(self.partial_pressures),
            "fugacity_residual": dict(self.fugacity_residual),
            "fugacity_residual_norm": self.fugacity_residual_norm,
            "state_failure_count": self.state_failure_count,
            "penalty_residuals": list(self.penalty_residuals),
            "diagnostics": _json_like(self.diagnostics),
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
    """Solve reactive speciation, then bubble pressure for the true-species liquid."""
    opts = _normalize_options(options)
    labels = [str(label) for label in species]
    pressure_seed = _positive_scalar(P_seed, "P_seed")
    strict = opts.error_mode == "raise"
    spec_options = opts.speciation_options or ReactiveSpeciationOptions()
    bubble_options = opts.bubble_options or ElectrolyteBubbleOptions(initial_pressure=pressure_seed)
    if not strict:
        spec_options = replace(spec_options, return_best_effort=True)
        bubble_options = replace(bubble_options, return_best_effort=True)
    if bubble_options.initial_pressure != pressure_seed:
        pressure_seed = float(bubble_options.initial_pressure)

    try:
        speciation = solve_reactive_speciation(
            species=labels,
            mixture_factory=mixture_factory,
            T=T,
            P=pressure_seed,
            balances=balances,
            totals=totals,
            reactions=reactions,
            initial_x=initial_x,
            options=spec_options,
        )
    except SolutionError as exc:
        if strict:
            raise
        diagnostics = {
            "success": False,
            "message": exc.message,
            "failure_stage": "reactive_speciation",
            "reactive_speciation_diagnostics": _json_like(exc.diagnostics or {}),
        }
        return _failed_result(
            labels=labels,
            message=exc.message,
            pressure=pressure_seed,
            penalty_value=opts.penalty_value,
            diagnostics=diagnostics,
        )

    x_array = np.asarray([speciation.x[label] for label in labels], dtype=float)
    mixture = mixture_factory(x_array, T, pressure_seed)
    try:
        bubble = mixture.equilibrium(
            kind="electrolyte_bubble_pressure",
            T=T,
            x_liq=x_array,
            volatile_species=vapor_species if volatile_species is None else volatile_species,
            vapor_species=vapor_species,
            nonvolatile_species=nonvolatile_species,
            options=bubble_options,
        )
    except SolutionError as exc:
        if strict:
            raise
        diagnostics = {
            "success": False,
            "message": exc.message,
            "failure_stage": "electrolyte_bubble_pressure",
            "reactive_speciation": speciation.to_dict(),
            "electrolyte_bubble_diagnostics": _json_like(exc.diagnostics or {}),
        }
        return _result_from_parts(
            success=False,
            message=exc.message,
            species=labels,
            speciation=speciation,
            bubble=None,
            pressure=pressure_seed,
            penalty_value=opts.penalty_value,
            diagnostics=diagnostics,
        )

    success = bool(speciation.success and bubble.success)
    message = "converged" if success else "reactive electrolyte bubble workflow returned a diagnostic result"
    diagnostics = {
        "success": success,
        "message": message,
        "failure_stage": "" if success else "nonconverged_subsolve",
        "reactive_speciation": speciation.to_dict(),
        "electrolyte_bubble_pressure": bubble.to_dict(),
        "error_mode": opts.error_mode,
    }
    return _result_from_parts(
        success=success,
        message=message,
        species=labels,
        speciation=speciation,
        bubble=bubble,
        pressure=bubble.P,
        penalty_value=opts.penalty_value,
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
    """Solve an ordered reactive electrolyte bubble sweep with continuation seeds."""
    if continuation not in {"auto", "none"}:
        raise InputError("continuation must be 'auto' or 'none'.")
    if not points:
        raise InputError("points must contain at least one sweep point.")
    opts = _normalize_options(options)
    labels = [str(label) for label in species]
    results: list[ReactiveElectrolyteBubbleResult] = []
    previous: ReactiveElectrolyteBubbleResult | None = None
    for index, point in enumerate(points):
        payload = dict(point)
        if "T" not in payload:
            raise InputError("Each reactive electrolyte sweep point must provide T.")
        if "totals" not in payload:
            raise InputError("Each reactive electrolyte sweep point must provide totals.")
        initial_x = payload.pop("initial_x", None)
        if continuation == "auto" and previous is not None and previous.x_liq:
            initial_x = [previous.x_liq[label] for label in labels]
        if initial_x is None:
            raise InputError("Each reactive electrolyte sweep point must provide initial_x or a continuation seed.")
        point_options = payload.pop("options", opts)
        point_options = _normalize_options(point_options)
        if continuation == "auto" and previous is not None and previous.P_total > 0.0:
            point_options = _options_with_continuation(point_options, previous)
        result = solve_reactive_electrolyte_bubble(
            species=labels,
            mixture_factory=mixture_factory,
            T=payload.pop("T"),
            P_seed=payload.pop("P_seed", previous.P_total if previous is not None else 1.0e5),
            balances=payload.pop("balances", balances),
            totals=payload.pop("totals"),
            reactions=payload.pop("reactions", reactions),
            initial_x=initial_x,
            vapor_species=payload.pop("vapor_species", vapor_species),
            volatile_species=payload.pop("volatile_species", volatile_species),
            nonvolatile_species=payload.pop("nonvolatile_species", nonvolatile_species),
            options=point_options,
        )
        diagnostics = dict(result.diagnostics)
        diagnostics["sweep_index"] = int(index)
        diagnostics["continuation_used"] = bool(continuation == "auto" and previous is not None)
        result = replace(result, diagnostics=diagnostics)
        results.append(result)
        if result.success:
            previous = result
    return results


def _normalize_options(options: ReactiveElectrolyteBubbleOptions | None) -> ReactiveElectrolyteBubbleOptions:
    if options is None:
        return ReactiveElectrolyteBubbleOptions()
    if not isinstance(options, ReactiveElectrolyteBubbleOptions):
        raise InputError("options must be a ReactiveElectrolyteBubbleOptions instance.")
    if options.error_mode not in {"raise", "result"}:
        raise InputError("ReactiveElectrolyteBubbleOptions.error_mode must be 'raise' or 'result'.")
    if options.penalty_value <= 0.0 or not np.isfinite(options.penalty_value):
        raise InputError("ReactiveElectrolyteBubbleOptions.penalty_value must be positive and finite.")
    if options.speciation_options is not None and not isinstance(options.speciation_options, ReactiveSpeciationOptions):
        raise InputError("ReactiveElectrolyteBubbleOptions.speciation_options must be ReactiveSpeciationOptions.")
    if options.bubble_options is not None and not isinstance(options.bubble_options, ElectrolyteBubbleOptions):
        raise InputError("ReactiveElectrolyteBubbleOptions.bubble_options must be ElectrolyteBubbleOptions.")
    return options


def _options_with_continuation(
    options: ReactiveElectrolyteBubbleOptions, previous: ReactiveElectrolyteBubbleResult
) -> ReactiveElectrolyteBubbleOptions:
    bubble = options.bubble_options or ElectrolyteBubbleOptions()
    bubble = replace(
        bubble,
        initial_pressure=float(previous.P_total),
        initial_y_vap=dict(previous.y_vap) if previous.y_vap else bubble.initial_y_vap,
    )
    return replace(options, bubble_options=bubble)


def _result_from_parts(
    *,
    success: bool,
    message: str,
    species: list[str],
    speciation: ReactiveSpeciationResult,
    bubble: ElectrolyteBubbleResult | None,
    pressure: float,
    penalty_value: float,
    diagnostics: dict[str, Any],
) -> ReactiveElectrolyteBubbleResult:
    penalty = [] if success else [float(penalty_value)]
    return ReactiveElectrolyteBubbleResult(
        success=success,
        message=message,
        x_liq={label: float(speciation.x.get(label, 0.0)) for label in species},
        activity_coefficients={label: float(speciation.activity_coefficients.get(label, 0.0)) for label in species},
        mass_balance_residuals=dict(speciation.mass_balance_residuals),
        charge_residual=speciation.charge_residual,
        reaction_residuals=list(speciation.reaction_residuals),
        named_reaction_residuals=dict(speciation.named_reaction_residuals),
        P_total=float(pressure),
        y_vap={} if bubble is None else dict(bubble.y_vap),
        partial_pressures={} if bubble is None else dict(bubble.partial_pressures),
        fugacity_residual={} if bubble is None else dict(bubble.fugacity_residual),
        fugacity_residual_norm=0.0 if bubble is None else bubble.fugacity_residual_norm,
        state_failure_count=int(
            speciation.state_failure_count + (0 if bubble is None else bubble.diagnostics.get("state_failure_count", 0))
        ),
        penalty_residuals=penalty,
        diagnostics=diagnostics,
    )


def _failed_result(
    *,
    labels: list[str],
    message: str,
    pressure: float,
    penalty_value: float,
    diagnostics: dict[str, Any],
) -> ReactiveElectrolyteBubbleResult:
    return ReactiveElectrolyteBubbleResult(
        success=False,
        message=message,
        x_liq={label: 0.0 for label in labels},
        activity_coefficients={label: 0.0 for label in labels},
        mass_balance_residuals={},
        charge_residual=0.0,
        reaction_residuals=[],
        named_reaction_residuals={},
        P_total=float(pressure),
        y_vap={},
        partial_pressures={},
        fugacity_residual={},
        fugacity_residual_norm=0.0,
        state_failure_count=0,
        penalty_residuals=[float(penalty_value)],
        diagnostics=diagnostics,
    )


def _positive_scalar(value: Any, name: str) -> float:
    try:
        out = float(value)
    except (TypeError, ValueError) as exc:
        raise InputError(f"{name} must be a positive finite number.") from exc
    if not np.isfinite(out) or out <= 0.0:
        raise InputError(f"{name} must be a positive finite number.")
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
