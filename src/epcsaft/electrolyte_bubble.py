"""Native-only electrolyte bubble-pressure public contracts."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import Any

from ._types import InputError


@dataclass(frozen=True, slots=True)
class ElectrolyteBubbleOptions:
    """Numerical controls for the future native electrolyte bubble-pressure backend."""

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
    """Structured result shape reserved for native electrolyte bubble-pressure calculations."""

    success: bool
    message: str
    P: float
    y_vap: Mapping[str, float]
    x_liq: Sequence[float]
    ln_phi_liq: Mapping[str, float]
    ln_phi_vap: Mapping[str, float]
    fugacity_residual: Mapping[str, float]
    fugacity_residual_norm: float
    charge_residual: float
    partial_pressures: Mapping[str, float]
    diagnostics: Mapping[str, Any]

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-like result payload."""
        return {
            "success": bool(self.success),
            "message": str(self.message),
            "P": float(self.P),
            "y_vap": {str(k): float(v) for k, v in self.y_vap.items()},
            "x_liq": [float(v) for v in self.x_liq],
            "ln_phi_liq": {str(k): float(v) for k, v in self.ln_phi_liq.items()},
            "ln_phi_vap": {str(k): float(v) for k, v in self.ln_phi_vap.items()},
            "fugacity_residual": {str(k): float(v) for k, v in self.fugacity_residual.items()},
            "fugacity_residual_norm": float(self.fugacity_residual_norm),
            "charge_residual": float(self.charge_residual),
            "partial_pressures": {str(k): float(v) for k, v in self.partial_pressures.items()},
            "diagnostics": dict(self.diagnostics),
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
    """Reject electrolyte bubble-pressure solves until a native backend exists."""
    _ = (mixture, T, x_liq, z, volatile_species, vapor_species, nonvolatile_species, options)
    raise InputError(
        "electrolyte_bubble_pressure is disabled until a native C++ backend is implemented; "
        "Python-side equilibrium solvers are not exposed by this package."
    )
