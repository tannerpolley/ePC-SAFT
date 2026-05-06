"""Native-only reactive electrolyte bubble workflow contracts."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import Any

from ._types import InputError
from .electrolyte_bubble import ElectrolyteBubbleOptions
from .reactive_speciation import ReactiveSpeciationOptions


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
    """Reject reactive electrolyte bubble solves until a native backend exists."""
    _ = (
        species,
        mixture_factory,
        T,
        P_seed,
        balances,
        totals,
        reactions,
        initial_x,
        vapor_species,
        volatile_species,
        nonvolatile_species,
        options,
    )
    raise InputError(
        "reactive electrolyte bubble pressure is disabled until a native C++ backend is implemented; "
        "Python-side equilibrium orchestration is not exposed by this package."
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
    """Reject reactive electrolyte bubble sweeps until a native backend exists."""
    _ = (
        species,
        mixture_factory,
        points,
        balances,
        reactions,
        vapor_species,
        volatile_species,
        nonvolatile_species,
        options,
        continuation,
    )
    raise InputError(
        "reactive electrolyte bubble sweeps are disabled until a native C++ backend is implemented; "
        "Python-side equilibrium orchestration is not exposed by this package."
    )
