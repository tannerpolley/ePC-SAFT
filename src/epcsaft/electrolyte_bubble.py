"""Native-only electrolyte bubble-pressure public contracts."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import Any

import numpy as np

from ._types import InputError


@dataclass(frozen=True, slots=True)
class ElectrolyteBubbleOptions:
    """Route controls reserved for native Ipopt fixed-liquid electrolyte bubble pressure."""

    initial_pressure: float = 1.0e5
    max_iterations: int = 80
    tolerance: float = 1.0e-6
    min_composition: float = 1.0e-14
    charge_tolerance: float = 1.0e-8
    initial_y_vap: Mapping[str, float] | None = None


@dataclass(frozen=True, slots=True)
class ElectrolyteBubbleResult:
    """Structured result returned by native electrolyte bubble-pressure calculations."""

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
    """Require the native Ipopt electrolyte bubble route builder before solving."""
    if options is None:
        options = ElectrolyteBubbleOptions()
    if not isinstance(options, ElectrolyteBubbleOptions):
        raise InputError("electrolyte_bubble_pressure options must be an ElectrolyteBubbleOptions instance.")
    if x_liq is None:
        if z is None:
            raise InputError("electrolyte_bubble_pressure requires x_liq or z.")
        x_liq = z
    species = list(getattr(mixture, "species", []))
    if not species:
        raise InputError("electrolyte_bubble_pressure requires mixture species labels.")
    vapor_labels = _normalize_species_labels(
        vapor_species if vapor_species is not None else volatile_species,
        field_name="vapor_species",
    )
    if not vapor_labels:
        raise InputError("electrolyte_bubble_pressure requires vapor_species or volatile_species.")
    if nonvolatile_species is not None:
        nonvolatile = set(_normalize_species_labels(nonvolatile_species, field_name="nonvolatile_species"))
        overlap = nonvolatile.intersection(vapor_labels)
        if overlap:
            raise InputError("vapor_species and nonvolatile_species overlap: " + ", ".join(sorted(overlap)))
    x_values = np.asarray(x_liq, dtype=float).flatten()
    if x_values.size != len(species):
        raise InputError("x_liq length must match mixture species count.")
    raise InputError(
        "electrolyte_bubble_pressure requires the native Ipopt equilibrium route builder; "
        "no package-owned alternate electrolyte bubble-pressure solver is available for this public route."
    )


def _normalize_species_labels(values: Any, *, field_name: str) -> tuple[str, ...]:
    if values is None:
        return ()
    if isinstance(values, str):
        return (values,)
    try:
        return tuple(str(value) for value in values)
    except TypeError as exc:
        raise InputError(f"{field_name} must be a string or sequence of strings.") from exc
