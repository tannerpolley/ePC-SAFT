"""Native-only reactive electrolyte bubble workflow contracts."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import Any

from ._types import InputError
from ._types import SolutionError
from .electrolyte_bubble import ElectrolyteBubbleOptions
from .electrolyte_bubble import electrolyte_bubble_pressure
from .reactive_speciation import ReactiveSpeciationOptions
from .reactive_speciation import solve_reactive_speciation


@dataclass(frozen=True, slots=True)
class ReactiveElectrolyteBubbleOptions:
    """Controls reserved for the future native reactive electrolyte bubble backend."""

    speciation_options: ReactiveSpeciationOptions | None = None
    bubble_options: ElectrolyteBubbleOptions | None = None
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
    except SolutionError:
        if options.error_mode == "raise":
            raise
        raise
    diagnostics = {
        "speciation": chemical.to_dict(),
        "bubble": bubble.to_dict(),
        "native_entrypoint": "_solve_chemical_equilibrium_native_then__solve_electrolyte_bubble_native",
        "solver_language": "c++",
    }
    return ReactiveElectrolyteBubbleResult(
        success=bool(chemical.success and bubble.success),
        message="converged" if chemical.success and bubble.success else bubble.message,
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
    _ = continuation
    if options is None:
        options = ReactiveElectrolyteBubbleOptions()
    results: list[ReactiveElectrolyteBubbleResult] = []
    last_pressure = None
    last_y = None
    for point in points:
        if "T" not in point or "totals" not in point:
            raise InputError("Each reactive electrolyte bubble sweep point requires T and totals.")
        point_options = options
        if options.bubble_options is not None and (last_pressure is not None or last_y is not None):
            base = options.bubble_options
            point_options = ReactiveElectrolyteBubbleOptions(
                speciation_options=options.speciation_options,
                bubble_options=ElectrolyteBubbleOptions(
                    initial_pressure=float(last_pressure if last_pressure is not None else base.initial_pressure),
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
                    initial_y_vap=last_y or base.initial_y_vap,
                ),
                error_mode=options.error_mode,
                penalty_value=options.penalty_value,
            )
        result = solve_reactive_electrolyte_bubble(
            species=species,
            mixture_factory=mixture_factory,
            T=float(point["T"]),
            P_seed=float(point.get("P_seed", last_pressure or point.get("P", 101325.0))),
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
