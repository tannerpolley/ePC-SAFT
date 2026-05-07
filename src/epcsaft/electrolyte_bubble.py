"""Native-only electrolyte bubble-pressure public contracts."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import Any

import numpy as np

from ._types import InputError
from ._types import SolutionError


@dataclass(frozen=True, slots=True)
class ElectrolyteBubbleOptions:
    """Numerical controls for native fixed-liquid electrolyte bubble pressure."""

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
    """Solve fixed-liquid electrolyte bubble pressure through the native backend."""
    from . import _core

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
    request = {
        "T": float(T),
        "x_liq": x_values.tolist(),
        "species": species,
        "vapor_species": list(vapor_labels),
        "options": _options_to_native_dict(options, vapor_species=vapor_labels),
    }
    try:
        payload = _core._solve_electrolyte_bubble_native(mixture._native, request)
    except _core.NativeValueError as exc:
        raise InputError(str(exc)) from exc
    result = _result_from_native_payload(
        payload,
        species=species,
        vapor_species=vapor_labels,
        x_liq=x_values,
    )
    if not result.success and not options.return_best_effort:
        raise SolutionError(result.message, result.diagnostics)
    return result


def _normalize_species_labels(values: Any, *, field_name: str) -> tuple[str, ...]:
    if values is None:
        return ()
    if isinstance(values, str):
        return (values,)
    try:
        return tuple(str(value) for value in values)
    except TypeError as exc:
        raise InputError(f"{field_name} must be a string or sequence of strings.") from exc


def _options_to_native_dict(options: ElectrolyteBubbleOptions, *, vapor_species: Sequence[str]) -> dict[str, Any]:
    initial_y = None
    if options.initial_y_vap is not None:
        if isinstance(options.initial_y_vap, Mapping):
            missing = [label for label in vapor_species if label not in options.initial_y_vap]
            if missing:
                raise InputError("initial_y_vap is missing vapor species: " + ", ".join(missing))
            initial_y = [float(options.initial_y_vap[label]) for label in vapor_species]
        else:
            initial_y = [float(value) for value in options.initial_y_vap]
    return {
        "initial_pressure": float(options.initial_pressure),
        "min_pressure": float(options.min_pressure),
        "max_pressure": float(options.max_pressure),
        "max_iterations": int(options.max_iterations),
        "max_vapor_iterations": int(options.max_vapor_iterations),
        "max_bracket_expansions": int(options.max_bracket_expansions),
        "tolerance": float(options.tolerance),
        "vapor_tolerance": float(options.vapor_tolerance),
        "pressure_factor": float(options.pressure_factor),
        "min_composition": float(options.min_composition),
        "charge_tolerance": float(options.charge_tolerance),
        "return_best_effort": bool(options.return_best_effort),
        "initial_y_vap": initial_y,
    }


def _phase_by_label(payload: Mapping[str, Any], label: str) -> Mapping[str, Any] | None:
    for phase in payload.get("phases", ()):
        if phase.get("label") == label:
            return phase
    return None


def _result_from_native_payload(
    payload: Mapping[str, Any],
    *,
    species: Sequence[str],
    vapor_species: Sequence[str],
    x_liq: np.ndarray,
) -> ElectrolyteBubbleResult:
    diagnostics = dict(payload.get("diagnostics") or {})
    success = bool(diagnostics.get("success", False))
    message = str(diagnostics.get("message", "converged" if success else "electrolyte bubble pressure failed"))
    vapor_indices = [species.index(label) for label in vapor_species]
    liquid = _phase_by_label(payload, "liq")
    vapor = _phase_by_label(payload, "vap")
    liquid_ln_phi = np.asarray((liquid or {}).get("ln_fugacity_coefficient", []), dtype=float)
    vapor_ln_phi = np.asarray((vapor or {}).get("ln_fugacity_coefficient", []), dtype=float)
    y_values = np.asarray(diagnostics.get("best_y_vap", []), dtype=float)
    if y_values.size != len(vapor_species) and vapor is not None:
        vapor_comp = np.asarray(vapor.get("composition", []), dtype=float)
        if vapor_comp.size == len(species):
            y_values = vapor_comp[vapor_indices]
            total = float(np.sum(y_values))
            if total > 0.0:
                y_values = y_values / total
    if y_values.size != len(vapor_species):
        y_values = np.zeros(len(vapor_species), dtype=float)
    partial = np.asarray(diagnostics.get("best_partial_pressures", []), dtype=float)
    if partial.size != len(vapor_species):
        pressure = float(diagnostics.get("best_P", 0.0))
        partial = y_values * pressure
    residual = np.asarray(diagnostics.get("fugacity_residual", []), dtype=float)
    if residual.size != len(vapor_species):
        residual = np.full(len(vapor_species), float("nan"))
    p_value = float(diagnostics.get("best_P", 0.0))
    ln_phi_liq = {
        label: float(liquid_ln_phi[index]) if liquid_ln_phi.size == len(species) else float("nan")
        for label, index in zip(vapor_species, vapor_indices)
    }
    if vapor_ln_phi.size == len(species):
        ln_phi_vap = {label: float(vapor_ln_phi[index]) for label, index in zip(vapor_species, vapor_indices)}
    elif vapor_ln_phi.size == len(vapor_species):
        ln_phi_vap = {label: float(vapor_ln_phi[pos]) for pos, label in enumerate(vapor_species)}
    else:
        ln_phi_vap = {label: float("nan") for label in vapor_species}
    return ElectrolyteBubbleResult(
        success=success,
        message=message,
        P=p_value,
        y_vap={label: float(value) for label, value in zip(vapor_species, y_values)},
        x_liq=[float(value) for value in x_liq],
        ln_phi_liq=ln_phi_liq,
        ln_phi_vap=ln_phi_vap,
        fugacity_residual={label: float(value) for label, value in zip(vapor_species, residual)},
        fugacity_residual_norm=float(diagnostics.get("best_fugacity_residual_norm", float("nan"))),
        charge_residual=float(diagnostics.get("charge_residual", float("nan"))),
        partial_pressures={label: float(value) for label, value in zip(vapor_species, partial)},
        diagnostics=diagnostics,
    )
