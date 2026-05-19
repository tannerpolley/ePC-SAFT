"""Native-only electrolyte bubble-pressure public contracts."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from numbers import Integral, Real
from typing import Any

import numpy as np

from ._types import InputError, SolutionError

_CANONICAL_PRESSURE_SCALE = 1.0e5


@dataclass(frozen=True, slots=True)
class ElectrolyteBubbleOptions:
    """Route controls reserved for native Ipopt fixed-liquid electrolyte bubble pressure."""

    max_iterations: int = 80
    tolerance: float = 1.0e-6
    min_composition: float = 1.0e-14
    charge_tolerance: float = 1.0e-8
    hessian_mode: str = "auto"
    ipopt_iteration_history_limit: int = 20
    ipopt_linear_solver: str = "auto"
    ipopt_acceptable_tolerance: float | None = None
    ipopt_constraint_violation_tolerance: float | None = None
    ipopt_dual_infeasibility_tolerance: float | None = None
    ipopt_complementarity_tolerance: float | None = None
    continuation_state: Mapping[str, Any] | None = None


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
    """Solve fixed-liquid electrolyte bubble pressure through the native Ipopt route."""
    from .equilibrium import (
        _ipopt_continuation_context,
        _prepare_ipopt_continuation_state,
        _stamp_ipopt_continuation_state,
    )
    from .equilibrium_core.native_results import native_route_diagnostics

    if options is None:
        options = ElectrolyteBubbleOptions()
    if not isinstance(options, ElectrolyteBubbleOptions):
        raise InputError("electrolyte_bubble_pressure options must be an ElectrolyteBubbleOptions instance.")
    hessian_mode = str(options.hessian_mode).strip().lower().replace("_", "-")
    if hessian_mode not in {"auto", "exact", "limited-memory"}:
        raise InputError("ElectrolyteBubbleOptions.hessian_mode must be 'auto', 'exact', or 'limited-memory'.")
    if isinstance(options.ipopt_iteration_history_limit, bool) or not isinstance(
        options.ipopt_iteration_history_limit, Integral
    ):
        raise InputError(
            "ElectrolyteBubbleOptions.ipopt_iteration_history_limit must be an integer greater than or equal to zero."
        )
    iteration_history_limit = int(options.ipopt_iteration_history_limit)
    if iteration_history_limit < 0:
        raise InputError(
            "ElectrolyteBubbleOptions.ipopt_iteration_history_limit must be an integer greater than or equal to zero."
        )
    ipopt_linear_solver = str(options.ipopt_linear_solver).strip().lower()
    if not ipopt_linear_solver:
        raise InputError("ElectrolyteBubbleOptions.ipopt_linear_solver must be a non-empty string.")
    ipopt_acceptable_tolerance = _optional_positive_float_option(
        options.ipopt_acceptable_tolerance,
        "ipopt_acceptable_tolerance",
    )
    ipopt_constraint_violation_tolerance = _optional_positive_float_option(
        options.ipopt_constraint_violation_tolerance,
        "ipopt_constraint_violation_tolerance",
    )
    ipopt_dual_infeasibility_tolerance = _optional_positive_float_option(
        options.ipopt_dual_infeasibility_tolerance,
        "ipopt_dual_infeasibility_tolerance",
    )
    ipopt_complementarity_tolerance = _optional_positive_float_option(
        options.ipopt_complementarity_tolerance,
        "ipopt_complementarity_tolerance",
    )
    if options.continuation_state is not None and not isinstance(options.continuation_state, Mapping):
        raise InputError("ElectrolyteBubbleOptions.continuation_state must be a mapping when provided.")
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
    if not np.all(np.isfinite(x_values)) or np.any(x_values <= 0.0):
        raise InputError("x_liq values must be positive and finite.")
    x_values = x_values / float(np.sum(x_values))
    unknown_vapor = sorted(set(vapor_labels).difference(species))
    if unknown_vapor:
        raise InputError("Unknown vapor species: " + ", ".join(unknown_vapor))
    continuation_context = _ipopt_continuation_context(
        route_kind="electrolyte_bubble_pressure",
        mixture=mixture,
        fixed_specs={
            "fixed": ["T", "x_liq"],
            "vapor_species": list(vapor_labels),
        },
    )
    prepared_continuation_state = _prepare_ipopt_continuation_state(
        ElectrolyteBubbleOptions(
            max_iterations=int(options.max_iterations),
            tolerance=float(options.tolerance),
            min_composition=float(options.min_composition),
            charge_tolerance=float(options.charge_tolerance),
            hessian_mode=hessian_mode,
            ipopt_iteration_history_limit=iteration_history_limit,
            ipopt_linear_solver=ipopt_linear_solver,
            ipopt_acceptable_tolerance=ipopt_acceptable_tolerance,
            ipopt_constraint_violation_tolerance=ipopt_constraint_violation_tolerance,
            ipopt_dual_infeasibility_tolerance=ipopt_dual_infeasibility_tolerance,
            ipopt_complementarity_tolerance=ipopt_complementarity_tolerance,
            continuation_state=options.continuation_state,
        ),
        continuation_context=continuation_context,
    )

    from . import _core

    route = _core._native_electrolyte_bubble_p_eos_route_result(
        mixture._native,
        float(T),
        x_values.tolist(),
        int(options.max_iterations),
        float(options.tolerance),
        hessian_mode,
        iteration_history_limit,
        float(options.tolerance),
        max(_CANONICAL_PRESSURE_SCALE * float(options.tolerance), float(options.tolerance)),
        float(options.charge_tolerance),
        float(options.tolerance),
        max(10.0 * float(options.min_composition), 1.0e-8),
        prepared_continuation_state,
        linear_solver=ipopt_linear_solver,
        acceptable_tolerance=_resolved_optional_tolerance(
            ipopt_acceptable_tolerance,
            max(100.0 * float(options.tolerance), 1.0e-10),
        ),
        constraint_violation_tolerance=_resolved_optional_tolerance(
            ipopt_constraint_violation_tolerance,
            float(options.tolerance),
        ),
        dual_infeasibility_tolerance=_resolved_optional_tolerance(
            ipopt_dual_infeasibility_tolerance,
            float(options.tolerance),
        ),
        complementarity_tolerance=_resolved_optional_tolerance(
            ipopt_complementarity_tolerance,
            float(options.tolerance),
        ),
    )
    if str(route.get("status", "")) == "ipopt_dependency_required":
        _raise_native_ipopt_electrolyte_bubble_required()
    if not bool(route.get("accepted", False)):
        diagnostics = native_route_diagnostics(route)
        _stamp_ipopt_continuation_state(diagnostics, route, continuation_context=continuation_context)
        raise SolutionError("Native electrolyte bubble-pressure route was rejected.", diagnostics)
    return _accepted_native_electrolyte_bubble_result(
        mixture,
        T=float(T),
        x_liq=x_values,
        vapor_labels=tuple(vapor_labels),
        route=route,
        continuation_context=continuation_context,
    )


def _raise_native_ipopt_electrolyte_bubble_required() -> None:
    raise InputError("electrolyte_bubble_pressure requires the native Ipopt equilibrium route builder.")


def _accepted_native_electrolyte_bubble_result(
    mixture: Any,
    *,
    T: float,
    x_liq: np.ndarray,
    vapor_labels: tuple[str, ...],
    route: Mapping[str, Any],
    continuation_context: Mapping[str, Any] | None = None,
) -> ElectrolyteBubbleResult:
    from .equilibrium import _stamp_ipopt_continuation_state
    from .equilibrium_core.native_results import native_route_diagnostics

    species = list(getattr(mixture, "species", []))
    phase_amounts = np.asarray(route.get("phase_amounts"), dtype=float)
    if phase_amounts.ndim != 2 or phase_amounts.shape != (2, len(species)):
        raise SolutionError("Native electrolyte bubble-pressure route returned invalid phase amounts.")
    phase_volumes = np.asarray(route.get("phase_volumes"), dtype=float)
    if phase_volumes.shape != (2,):
        raise SolutionError("Native electrolyte bubble-pressure route returned invalid phase volumes.")
    variables = np.asarray(route.get("variables"), dtype=float).flatten()
    if variables.size == 0:
        raise SolutionError("Native electrolyte bubble-pressure route returned no solver variables.")
    pressure = float(variables[-1])
    if not np.isfinite(pressure) or pressure <= 0.0:
        raise SolutionError("Native electrolyte bubble-pressure route returned invalid pressure.")

    liquid = phase_amounts[0] / float(np.sum(phase_amounts[0]))
    vapor = phase_amounts[1] / float(np.sum(phase_amounts[1]))
    liquid_state = mixture.state(T=T, P=pressure, x=liquid, phase="liq")
    vapor_state = mixture.state(T=T, P=pressure, x=vapor, phase="vap")
    ln_phi_liq = np.asarray(liquid_state.fugacity_coefficient(), dtype=float)
    ln_phi_vap = np.asarray(vapor_state.fugacity_coefficient(), dtype=float)

    vapor_indices = [species.index(label) for label in vapor_labels]
    y_vap = {species[index]: float(vapor[index]) for index in vapor_indices}
    ln_liq = {species[index]: float(ln_phi_liq[index]) for index in vapor_indices}
    ln_vap = {species[index]: float(ln_phi_vap[index]) for index in vapor_indices}
    partial_pressures = {species[index]: float(pressure * vapor[index]) for index in vapor_indices}
    fugacity_residual = {
        species[index]: float(np.log(liquid[index]) + ln_phi_liq[index] - np.log(vapor[index]) - ln_phi_vap[index])
        for index in vapor_indices
    }
    fugacity_residual_norm = max((abs(value) for value in fugacity_residual.values()), default=0.0)
    diagnostics = native_route_diagnostics(route)
    diagnostics["backend"] = "ipopt"
    if continuation_context is not None:
        _stamp_ipopt_continuation_state(diagnostics, route, continuation_context=continuation_context)
    return ElectrolyteBubbleResult(
        success=True,
        message="converged",
        P=pressure,
        y_vap=y_vap,
        x_liq=[float(value) for value in x_liq],
        ln_phi_liq=ln_liq,
        ln_phi_vap=ln_vap,
        fugacity_residual=fugacity_residual,
        fugacity_residual_norm=float(fugacity_residual_norm),
        charge_residual=float(diagnostics.get("charge_balance_norm", 0.0)),
        partial_pressures=partial_pressures,
        diagnostics=diagnostics,
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


def _optional_positive_float_option(value: Any, label: str) -> float | None:
    if value is None:
        return None
    if isinstance(value, bool) or not isinstance(value, Real):
        raise InputError(f"ElectrolyteBubbleOptions.{label} must be a finite real number when provided.")
    out = float(value)
    if not np.isfinite(out) or out <= 0.0:
        raise InputError(f"ElectrolyteBubbleOptions.{label} must be positive when provided.")
    return out


def _resolved_optional_tolerance(value: float | None, default: float) -> float:
    if value is None:
        return float(default)
    return float(value)
