"""Native-backed phase-equilibrium helpers and Python input adapters."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from numbers import Integral, Real
from typing import Any, Literal

import numpy as np

from ._types import InputError, SolutionError
from .equilibrium_core.electrolyte_basis import (
    electrolyte_feed_from_molality_inputs,
    electrolyte_formula_basis,
    normalize_salt_molality,
)
from .equilibrium_core.native_requests import neutral_two_phase_eos_tolerances
from .equilibrium_core.native_results import (
    native_route_solved_pressure,
    native_route_summed_phase_amounts,
    neutral_two_phase_payload_to_result,
)

_ASCANI_2022_REFERENCE = {
    "authors": "Ascani, Sadowski, and Held",
    "year": 2022,
    "title": "Calculation of Multiphase Equilibria Containing Mixed Solvents and Mixed Electrolytes",
    "doi": "10.1021/acs.jced.1c00866",
}


def _raise_native_ipopt_equilibrium_required(route: str) -> None:
    raise InputError(f"{route} requires a native Ipopt equilibrium NLP route.")


def _raise_native_ipopt_reactive_phase_required(route: str) -> None:
    raise InputError(f"{route} requires a native Ipopt reactive phase-equilibrium NLP route.")


def _raise_native_ipopt_lle_required(route: str) -> None:
    raise InputError(f"{route} requires a native Ipopt equilibrium NLP route.")


def _raise_native_ipopt_tp_flash_required() -> None:
    raise InputError("tp_flash requires a native Ipopt equilibrium NLP route.")


def _raise_native_ipopt_stability_required(route: str) -> None:
    raise InputError(f"{route} requires a native Ipopt equilibrium stability NLP route.")


@dataclass(frozen=True, slots=True)
class EquilibriumOptions:
    """Numerical controls for equilibrium solvers."""

    max_iterations: int = 180
    tolerance: float = 1.0e-6
    min_composition: float = 1.0e-12
    jacobian_backend: Literal["auto", "analytic", "cppad"] = "auto"
    solver_backend: Literal["auto", "ipopt"] = "auto"
    timeout_seconds: float | None = None


@dataclass(frozen=True, slots=True)
class EquilibriumProblem:
    """Base class for typed equilibrium problem objects."""

    def solve(self, mixture):
        """Solve this problem with a mixture instance."""
        raise InputError("EquilibriumProblem is abstract; instantiate a concrete equilibrium problem class.")


@dataclass(frozen=True, slots=True)
class TPFlash(EquilibriumProblem):
    """Neutral TP flash problem."""

    T: float
    P: float
    z: Any
    options: EquilibriumOptions | None = None

    def solve(self, mixture):
        return mixture.flash_tp(T=self.T, P=self.P, z=self.z, options=self.options)


@dataclass(frozen=True, slots=True)
class StabilityAnalysis(EquilibriumProblem):
    """Neutral or electrolyte tangent-plane-distance stability problem."""

    T: float
    P: float
    z: Any | None = None
    options: EquilibriumOptions | None = None
    parent_phase: str | None = None
    trial_phases: Any | None = None
    solvent_feed: Mapping[str, float] | None = None
    salt_molality: Mapping[str, float] | None = None
    electrolyte: bool = False

    def solve(self, mixture):
        if self.electrolyte or self.solvent_feed is not None or self.salt_molality is not None:
            return mixture.electrolyte_stability_tp(
                T=self.T,
                P=self.P,
                z=self.z,
                solvent_feed=self.solvent_feed,
                salt_molality=self.salt_molality,
                options=self.options,
            )
        return mixture.stability_tp(
            T=self.T,
            P=self.P,
            z=self.z,
            options=self.options,
            parent_phase=self.parent_phase,
            trial_phases=self.trial_phases,
        )


@dataclass(frozen=True, slots=True)
class BubblePoint(EquilibriumProblem):
    """Neutral bubble-point problem."""

    T: float
    x: Any
    options: EquilibriumOptions | None = None

    def solve(self, mixture):
        return mixture.bubble_p(T=self.T, x=self.x, options=self.options)


@dataclass(frozen=True, slots=True)
class DewPoint(EquilibriumProblem):
    """Neutral dew-point problem."""

    y: Any
    T: float | None = None
    P: float | None = None
    options: EquilibriumOptions | None = None

    def solve(self, mixture):
        if (self.T is None) == (self.P is None):
            raise InputError("DewPoint requires exactly one of T or P.")
        if self.T is not None:
            return mixture.dew_p(T=self.T, y=self.y, options=self.options)
        return mixture.dew_t(P=self.P, y=self.y, options=self.options)


@dataclass(frozen=True, slots=True)
class LLEProblem(EquilibriumProblem):
    """Neutral liquid-liquid equilibrium problem."""

    T: float
    P: float
    z: Any
    options: EquilibriumOptions | None = None

    def solve(self, mixture):
        return mixture.lle_tp(
            T=self.T,
            P=self.P,
            z=self.z,
            options=self.options,
        )


@dataclass(frozen=True, slots=True)
class ElectrolyteLLEProblem(EquilibriumProblem):
    """Charge-constrained electrolyte LLE problem."""

    T: float
    P: float
    z: Any | None = None
    solvent_feed: Mapping[str, float] | None = None
    salt_molality: Mapping[str, float] | None = None
    options: EquilibriumOptions | None = None

    def solve(self, mixture):
        return mixture.electrolyte_lle_tp(
            T=self.T,
            P=self.P,
            z=self.z,
            solvent_feed=self.solvent_feed,
            salt_molality=self.salt_molality,
            options=self.options,
        )


@dataclass(frozen=True, slots=True)
class ElectrolyteBubblePoint(EquilibriumProblem):
    """Fixed-liquid electrolyte bubble-point problem."""

    T: float
    x_liq: Any | None = None
    z: Any | None = None
    vapor_species: Any | None = None
    volatile_species: Any | None = None
    nonvolatile_species: Any | None = None
    options: Any | None = None

    def solve(self, mixture):
        return mixture.electrolyte_bubble_p(
            T=self.T,
            x_liq=self.x_liq,
            z=self.z,
            vapor_species=self.vapor_species,
            volatile_species=self.volatile_species,
            nonvolatile_species=self.nonvolatile_species,
            options=self.options,
        )


@dataclass(frozen=True, slots=True)
class ReactiveSpeciationProblem(EquilibriumProblem):
    """Homogeneous reactive speciation problem."""

    T: float
    P: float
    balances: Any
    totals: Mapping[str, float]
    reactions: Any
    initial_x: Any | None = None
    z: Any | None = None
    options: Any | None = None

    def solve(self, mixture):
        return mixture.chemical_equilibrium(
            T=self.T,
            P=self.P,
            balances=self.balances,
            totals=self.totals,
            reactions=self.reactions,
            initial_x=self.initial_x,
            z=self.z,
            options=self.options,
        )


@dataclass(frozen=True, slots=True)
class ReactiveElectrolyteBubbleProblem(ReactiveSpeciationProblem):
    """Reactive speciation followed by fixed-liquid electrolyte bubble pressure."""

    vapor_species: Any | None = None
    volatile_species: Any | None = None
    nonvolatile_species: Any | None = None

    def solve(self, mixture):
        return mixture.reactive_electrolyte_bubble_p(
            T=self.T,
            P=self.P,
            balances=self.balances,
            totals=self.totals,
            reactions=self.reactions,
            initial_x=self.initial_x,
            z=self.z,
            vapor_species=self.vapor_species,
            volatile_species=self.volatile_species,
            nonvolatile_species=self.nonvolatile_species,
            options=self.options,
        )


@dataclass(frozen=True, slots=True)
class ReactivePhaseEquilibriumProblem(ReactiveSpeciationProblem):
    """Coupled native reactive phase-equilibrium problem."""

    phase_kind: str = "auto"
    phase_options: Any | None = None
    phase_kwargs: Mapping[str, Any] | None = None

    def solve(self, mixture):
        return reactive_phase_equilibrium(
            mixture,
            T=self.T,
            P=self.P,
            balances=self.balances,
            totals=self.totals,
            reactions=self.reactions,
            initial_x=self.initial_x,
            z=self.z,
            phase_kind=self.phase_kind,
            options=self.options,
            phase_options=self.phase_options,
            phase_kwargs=self.phase_kwargs,
        )


@dataclass(frozen=True, slots=True, init=False)
class EquilibriumPhase:
    """One phase returned by an equilibrium calculation."""

    label: str
    composition: np.ndarray
    density: float
    temperature: float
    pressure: float
    phase_fraction: float
    ln_fugacity_coefficient: np.ndarray | None = None
    diagnostics: dict[str, Any] | None = None

    def __init__(
        self,
        label: str,
        composition: Any,
        density: float,
        temperature: float,
        pressure: float,
        phase_fraction: float,
        ln_fugacity_coefficient: Any = None,
        fugacity_coefficient: Any = None,
        diagnostics: dict[str, Any] | None = None,
    ) -> None:
        """Create a phase payload.

        ``ln_fugacity_coefficient`` is the explicit natural-log fugacity
        coefficient field. ``fugacity_coefficient`` accepts coefficient-form
        phi values and is converted to ``ln(phi)`` when the log field is not
        supplied.
        """
        if ln_fugacity_coefficient is None:
            if fugacity_coefficient is None:
                ln_fugacity_coefficient = None
            else:
                ln_fugacity_coefficient = np.log(np.asarray(fugacity_coefficient, dtype=float))
        object.__setattr__(self, "label", str(label))
        object.__setattr__(self, "composition", np.asarray(composition, dtype=float))
        object.__setattr__(self, "density", float(density))
        object.__setattr__(self, "temperature", float(temperature))
        object.__setattr__(self, "pressure", float(pressure))
        object.__setattr__(self, "phase_fraction", float(phase_fraction))
        if ln_fugacity_coefficient is not None:
            ln_fugacity_coefficient = np.asarray(ln_fugacity_coefficient, dtype=float)
        object.__setattr__(self, "ln_fugacity_coefficient", ln_fugacity_coefficient)
        object.__setattr__(self, "diagnostics", None if diagnostics is None else dict(diagnostics))

    @property
    def fugacity_coefficient(self) -> np.ndarray | None:
        """Return coefficient-form fugacity coefficients."""
        if self.ln_fugacity_coefficient is None:
            return None
        return np.exp(self.ln_fugacity_coefficient)

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-like phase payload."""
        ln_fugacity = None if self.ln_fugacity_coefficient is None else self.ln_fugacity_coefficient.tolist()
        return {
            "label": self.label,
            "composition": self.composition.tolist(),
            "density": self.density,
            "temperature": self.temperature,
            "pressure": self.pressure,
            "phase_fraction": self.phase_fraction,
            "ln_fugacity_coefficient": ln_fugacity,
            "fugacity_coefficient": None if self.fugacity_coefficient is None else self.fugacity_coefficient.tolist(),
            "diagnostics": _json_like(self.diagnostics),
        }


@dataclass(frozen=True, slots=True)
class EquilibriumResult:
    """Structured result returned by an equilibrium calculation."""

    backend: str
    problem_kind: str
    phases: tuple[EquilibriumPhase, ...]
    stable: bool
    split_detected: bool
    diagnostics: dict[str, Any]

    def __post_init__(self) -> None:
        object.__setattr__(self, "backend", str(self.backend))
        object.__setattr__(self, "problem_kind", str(self.problem_kind))
        object.__setattr__(self, "phases", tuple(self.phases))
        object.__setattr__(self, "stable", bool(self.stable))
        object.__setattr__(self, "split_detected", bool(self.split_detected))
        object.__setattr__(self, "diagnostics", dict(self.diagnostics))

    @property
    def phase_labels(self) -> list[str]:
        """Return phase labels in result order."""
        return [phase.label for phase in self.phases]

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-like result payload."""
        return {
            "backend": self.backend,
            "problem_kind": self.problem_kind,
            "phase_labels": self.phase_labels,
            "phases": [phase.to_dict() for phase in self.phases],
            "stable": self.stable,
            "split_detected": self.split_detected,
            "diagnostics": _json_like(self.diagnostics),
        }


@dataclass(frozen=True, slots=True)
class StabilityTrial:
    """One tangent-plane-distance trial calculation."""

    parent_phase: str
    trial_phase: str
    seed_name: str
    composition: np.ndarray
    tpd: float
    iterations: int
    converged: bool
    unstable: bool
    diagnostics: dict[str, Any]

    def __post_init__(self) -> None:
        object.__setattr__(self, "parent_phase", str(self.parent_phase))
        object.__setattr__(self, "trial_phase", str(self.trial_phase))
        object.__setattr__(self, "seed_name", str(self.seed_name))
        object.__setattr__(self, "composition", np.asarray(self.composition, dtype=float))
        object.__setattr__(self, "tpd", float(self.tpd))
        object.__setattr__(self, "iterations", int(self.iterations))
        object.__setattr__(self, "converged", bool(self.converged))
        object.__setattr__(self, "unstable", bool(self.unstable))
        object.__setattr__(self, "diagnostics", dict(self.diagnostics))

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-like stability-trial payload."""
        return {
            "parent_phase": self.parent_phase,
            "trial_phase": self.trial_phase,
            "seed_name": self.seed_name,
            "composition": self.composition.tolist(),
            "tpd": self.tpd,
            "iterations": self.iterations,
            "converged": self.converged,
            "unstable": self.unstable,
            "diagnostics": _json_like(self.diagnostics),
        }


@dataclass(frozen=True, slots=True)
class StabilityResult:
    """Structured result returned by neutral TPD stability analysis."""

    backend: str
    problem_kind: str
    stable: bool
    min_tpd: float
    parent_phase: str
    trial_phase: str
    trial_composition: np.ndarray
    trials: tuple[StabilityTrial, ...]
    diagnostics: dict[str, Any]

    def __post_init__(self) -> None:
        object.__setattr__(self, "backend", str(self.backend))
        object.__setattr__(self, "problem_kind", str(self.problem_kind))
        object.__setattr__(self, "stable", bool(self.stable))
        object.__setattr__(self, "min_tpd", float(self.min_tpd))
        object.__setattr__(self, "parent_phase", str(self.parent_phase))
        object.__setattr__(self, "trial_phase", str(self.trial_phase))
        object.__setattr__(self, "trial_composition", np.asarray(self.trial_composition, dtype=float))
        object.__setattr__(self, "trials", tuple(self.trials))
        object.__setattr__(self, "diagnostics", dict(self.diagnostics))

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-like stability-result payload."""
        return {
            "backend": self.backend,
            "problem_kind": self.problem_kind,
            "stable": self.stable,
            "min_tpd": self.min_tpd,
            "parent_phase": self.parent_phase,
            "trial_phase": self.trial_phase,
            "trial_composition": self.trial_composition.tolist(),
            "trials": [trial.to_dict() for trial in self.trials],
            "diagnostics": _json_like(self.diagnostics),
        }


def electrolyte_lle_flash(
    mixture: Any,
    *,
    T: float,
    P: float,
    z: Any = None,
    solvent_feed: Any = None,
    salt_molality: Any = None,
    options: EquilibriumOptions | None = None,
) -> EquilibriumResult:
    """Run the public electrolyte LLE route."""
    return electrolyte_lle_flash_native(
        mixture,
        T=T,
        P=P,
        z=z,
        solvent_feed=solvent_feed,
        salt_molality=salt_molality,
        options=options,
    )


def electrolyte_feed_from_molality(
    mixture: Any,
    *,
    solvent_feed: Any,
    salt_molality: Any,
    basis_mass_kg: float = 1.0,
) -> np.ndarray:
    """Convert mixed-solvent salt molality input into a species mole-fraction feed."""
    species = list(mixture.species)
    charges = _mixture_charges(mixture)
    mw = np.asarray(mixture.parameters.get("MW", []), dtype=float).flatten()
    if mw.size != len(species):
        raise InputError("mixture parameters must include one MW value per species.")
    return electrolyte_feed_from_molality_inputs(
        species,
        charges,
        mw,
        solvent_feed=solvent_feed,
        salt_molality=salt_molality,
        basis_mass_kg=basis_mass_kg,
    )


def _normalize_options(options: EquilibriumOptions | Mapping[str, Any] | None) -> EquilibriumOptions:
    if options is None:
        return EquilibriumOptions()
    if isinstance(options, Mapping):
        raw = dict(options)
        allowed = {
            "max_iterations",
            "tolerance",
            "min_composition",
            "jacobian_backend",
            "solver_backend",
            "timeout_seconds",
        }
        unknown = sorted(set(raw) - allowed)
        if unknown:
            raise InputError("Unknown equilibrium option key(s): {}.".format(", ".join(unknown)))
        options = EquilibriumOptions(**raw)
    if not isinstance(options, EquilibriumOptions):
        raise InputError("options must be an EquilibriumOptions instance.")
    if isinstance(options.max_iterations, bool) or not isinstance(options.max_iterations, Integral):
        raise InputError("options.max_iterations must be an integer greater than zero.")
    max_iterations = int(options.max_iterations)
    if max_iterations <= 0:
        raise InputError("options.max_iterations must be an integer greater than zero.")
    tolerance = _finite_float_option(options.tolerance, "tolerance")
    if tolerance <= 0.0:
        raise InputError("options.tolerance must be positive.")
    min_composition = _finite_float_option(options.min_composition, "min_composition")
    if min_composition <= 0.0:
        raise InputError("options.min_composition must be positive.")
    jacobian_backend = str(options.jacobian_backend).strip().lower()
    if jacobian_backend not in {"auto", "analytic", "cppad"}:
        raise InputError("options.jacobian_backend must be 'auto', 'analytic', or 'cppad'.")
    solver_backend = str(options.solver_backend).strip().lower()
    if solver_backend not in {"auto", "ipopt"}:
        raise InputError("options.solver_backend must be 'auto' or 'ipopt'.")
    timeout_seconds = _optional_positive_float_option(options.timeout_seconds, "timeout_seconds")
    return EquilibriumOptions(
        max_iterations=max_iterations,
        tolerance=tolerance,
        min_composition=min_composition,
        jacobian_backend=jacobian_backend,  # type: ignore[arg-type]
        solver_backend=solver_backend,  # type: ignore[arg-type]
        timeout_seconds=timeout_seconds,
    )


def _finite_float_option(value: Any, label: str) -> float:
    if isinstance(value, bool) or not isinstance(value, Real):
        raise InputError(f"options.{label} must be a finite real number.")
    out = float(value)
    if not np.isfinite(out):
        raise InputError(f"options.{label} must be finite.")
    return out


def _optional_positive_float_option(value: Any, label: str) -> float | None:
    if value is None:
        return None
    out = _finite_float_option(value, label)
    if out <= 0.0:
        raise InputError(f"options.{label} must be positive when provided.")
    return out


def _native_timeout_seconds(options: EquilibriumOptions) -> float:
    return 0.0 if options.timeout_seconds is None else float(options.timeout_seconds)


def _normalize_feed(z: Any, ncomp: int, min_composition: float, kind: str) -> np.ndarray:
    if z is None:
        raise InputError(f"z is required for kind='{kind}'.")
    feed = np.asarray(z, dtype=float).flatten()
    if feed.size != int(ncomp):
        raise InputError(f"Feed composition length ({feed.size}) must match mixture component count ({ncomp}).")
    if not np.all(np.isfinite(feed)):
        raise InputError("Feed composition z must contain only finite values.")
    if np.any(feed < 0.0):
        raise InputError("Feed composition z must be non-negative.")
    total = float(np.sum(feed))
    if total <= 0.0:
        raise InputError("Feed composition z must have a positive sum.")
    feed = feed / total
    if np.any(feed < min_composition):
        raise InputError(f"{kind} requires each feed composition entry to be >= min_composition.")
    return feed


def _positive_scalar(value: Any, label: str, kind: str) -> float:
    if value is None:
        raise InputError(f"{label} is required for kind='{kind}'.")
    out = float(value)
    if not np.isfinite(out) or out <= 0.0:
        raise InputError(f"{label} must be a positive finite scalar.")
    return out


def _reject_ion_containing_mixture(mixture: Any) -> None:
    z = np.asarray(mixture.parameters.get("z", []), dtype=float).flatten()
    if z.size and np.any(np.abs(z) > 1.0e-12):
        raise InputError("Neutral equilibrium does not support ion-containing mixtures.")


def _require_ion_containing_mixture(mixture: Any, kind: str) -> None:
    charges = _mixture_charges(mixture)
    if not np.any(np.abs(charges) > 1.0e-12):
        raise InputError(f"{kind} requires an ion-containing mixture.")


def _mixture_charges(mixture: Any) -> np.ndarray:
    charges = np.asarray(mixture.parameters.get("z", []), dtype=float).flatten()
    if charges.size != int(mixture.ncomp):
        raise InputError("mixture parameters must include one charge value per species in params['z'].")
    return charges


def _require_charge_neutral(composition: np.ndarray, charges: np.ndarray, label: str) -> None:
    charge = float(np.dot(np.asarray(composition, dtype=float), np.asarray(charges, dtype=float)))
    if abs(charge) > 1.0e-10:
        raise InputError(f"{label} must be charge neutral; charge balance is {charge}.")


def _normalize_electrolyte_feed(
    mixture: Any,
    *,
    z: Any,
    solvent_feed: Any,
    salt_molality: Any,
    options: EquilibriumOptions,
) -> tuple[np.ndarray, dict[str, Any]]:
    if z is not None and (solvent_feed is not None or salt_molality is not None):
        raise InputError("Use either direct mole-fraction z or solvent_feed plus salt_molality, not both.")
    if z is not None:
        feed = _normalize_feed(z, mixture.ncomp, options.min_composition, "electrolyte_lle")
        return feed, {"composition_basis": "mole_fraction"}
    if solvent_feed is None or salt_molality is None:
        raise InputError("electrolyte_lle requires z or solvent_feed plus salt_molality.")
    feed = electrolyte_feed_from_molality(mixture, solvent_feed=solvent_feed, salt_molality=salt_molality)
    return feed, {
        "composition_basis": "molality",
        "salt_molality": dict(normalize_salt_molality(salt_molality)),
        "solvent_feed": _json_like(solvent_feed),
    }


def _accepted_native_neutral_two_phase_result(
    mixture: Any,
    *,
    T: float,
    P: float,
    feed: np.ndarray,
    route: Mapping[str, Any],
    tolerances: tuple[float, float, float, float],
    route_label: str,
    problem_kind: str,
    phase_labels: tuple[str, str],
    route_family: str = "neutral",
) -> EquilibriumResult:
    from . import _core

    material_tolerance, pressure_tolerance, chemical_potential_tolerance, phase_distance_tolerance = tolerances
    if not bool(route.get("accepted", False)):
        postsolve = route.get("postsolve", {})
        diagnostics = dict(postsolve) if isinstance(postsolve, Mapping) else {}
        if route_status := route.get("status"):
            diagnostics["route_status"] = route_status
        if solver_status := route.get("solver_status"):
            diagnostics["solver_status"] = solver_status
        raise SolutionError(f"Native {route_family} {route_label} route was rejected.", diagnostics)

    result_payload = _core._native_neutral_two_phase_eos_result(
        mixture._native,
        T,
        P,
        route["phase_amounts"],
        route["phase_volumes"],
        feed.tolist(),
        material_tolerance,
        pressure_tolerance,
        chemical_potential_tolerance,
        phase_distance_tolerance,
    )
    return neutral_two_phase_payload_to_result(
        result_payload,
        problem_kind=problem_kind,
        phase_labels=phase_labels,
    )


def _native_neutral_fixed_temperature_pressure(
    mixture: Any,
    *,
    T: float,
    composition: np.ndarray,
    options: EquilibriumOptions,
    route_label: str,
    route_binding: str,
    problem_kind: str,
) -> EquilibriumResult:
    from . import _core

    route_tolerances = (
        options.tolerance,
        max(1.0e5 * options.tolerance, options.tolerance),
        options.tolerance,
        max(10.0 * options.min_composition, 1.0e-8),
    )
    route = getattr(_core, route_binding)(
        mixture._native,
        T,
        composition.tolist(),
        options.max_iterations,
        options.tolerance,
        _native_timeout_seconds(options),
        *route_tolerances,
    )
    if str(route.get("status", "")) == "ipopt_dependency_required":
        _raise_native_ipopt_equilibrium_required(route_label)

    pressure = native_route_solved_pressure(route, route_label) if bool(route.get("accepted", False)) else 1.0
    feed = (
        native_route_summed_phase_amounts(route, mixture.ncomp, route_label)
        if bool(route.get("accepted", False))
        else composition
    )
    return _accepted_native_neutral_two_phase_result(
        mixture,
        T=T,
        P=pressure,
        feed=feed,
        route=route,
        tolerances=neutral_two_phase_eos_tolerances(pressure, options),
        route_label=route_label,
        problem_kind=problem_kind,
        phase_labels=("liq", "vap"),
    )


def _native_neutral_tp_flash(
    mixture: Any,
    *,
    T: float,
    P: float,
    feed: np.ndarray,
    options: EquilibriumOptions,
) -> EquilibriumResult:
    from . import _core

    tolerances = neutral_two_phase_eos_tolerances(P, options)
    route = _core._native_neutral_tp_flash_eos_route_result(
        mixture._native,
        T,
        P,
        feed.tolist(),
        options.max_iterations,
        options.tolerance,
        _native_timeout_seconds(options),
        *tolerances,
    )
    if str(route.get("status", "")) == "ipopt_dependency_required":
        _raise_native_ipopt_tp_flash_required()
    return _accepted_native_neutral_two_phase_result(
        mixture,
        T=T,
        P=P,
        feed=feed,
        route=route,
        tolerances=tolerances,
        route_label="TP flash",
        problem_kind="neutral_tp_flash",
        phase_labels=("phase_0", "phase_1"),
    )


def _native_neutral_lle_flash(
    mixture: Any,
    *,
    T: float,
    P: float,
    feed: np.ndarray,
    options: EquilibriumOptions,
) -> EquilibriumResult:
    from . import _core

    tolerances = neutral_two_phase_eos_tolerances(P, options)
    route = _core._native_neutral_lle_eos_route_result(
        mixture._native,
        T,
        P,
        feed.tolist(),
        options.max_iterations,
        options.tolerance,
        _native_timeout_seconds(options),
        *tolerances,
    )
    if str(route.get("status", "")) == "ipopt_dependency_required":
        _raise_native_ipopt_lle_required("lle_flash")
    return _accepted_native_neutral_two_phase_result(
        mixture,
        T=T,
        P=P,
        feed=feed,
        route=route,
        tolerances=tolerances,
        route_label="LLE",
        problem_kind="neutral_lle",
        phase_labels=("liq1", "liq2"),
    )


def _normalize_parent_phases(parent_phase: Any) -> tuple[str, ...]:
    if parent_phase is None:
        return ("liq", "vap")
    return (_normalize_phase_token(parent_phase, "parent_phase"),)


def _normalize_trial_phases(trial_phases: Any) -> tuple[str, ...]:
    if trial_phases is None:
        return ("liq", "vap")
    if isinstance(trial_phases, str):
        return (_normalize_phase_token(trial_phases, "trial_phases"),)
    try:
        tokens = tuple(_normalize_phase_token(item, "trial_phases") for item in trial_phases)
    except TypeError as exc:
        raise InputError("trial_phases must be None, a phase string, or an iterable of phase strings.") from exc
    if not tokens:
        raise InputError("trial_phases must contain at least one phase.")
    return tokens


def _normalize_phase_token(value: Any, label: str) -> str:
    token = str(value).strip().lower()
    if token not in {"liq", "vap"}:
        raise InputError(f"{label} must be None, 'liq', or 'vap'.")
    return token


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


def reactive_phase_equilibrium(
    mixture: Any,
    *,
    T: float,
    P: float,
    balances: Mapping[str, Mapping[str, float]],
    totals: Mapping[str, float],
    reactions: Any,
    initial_x: Any = None,
    z: Any = None,
    phase_kind: str = "auto",
    options: Any = None,
    phase_options: EquilibriumOptions | Mapping[str, Any] | None = None,
    phase_kwargs: Mapping[str, Any] | None = None,
) -> EquilibriumResult:
    """Validate a reactive phase-equilibrium request and require the native Ipopt route."""
    from .reactive_speciation import (
        _normalize_balances as _normalize_reactive_balances,
    )
    from .reactive_speciation import (
        _normalize_options as _normalize_reactive_options,
    )
    from .reactive_speciation import (
        _normalize_reactions,
    )

    species = [str(label) for label in getattr(mixture, "species", [])]
    if not species:
        raise InputError("reactive phase equilibrium requires mixture species.")
    route = _normalize_reactive_phase_route(mixture, phase_kind, phase_kwargs)
    extra_phase_kwargs = dict(phase_kwargs or {})
    _reject_reactive_phase_kwargs(extra_phase_kwargs, route)
    _, solver_options = _reactive_phase_option_pair(
        options=options,
        phase_options=phase_options,
        normalize_reactive_options=_normalize_reactive_options,
    )
    _positive_scalar(T, "T", "reactive_phase_equilibrium")
    _positive_scalar(P, "P", "reactive_phase_equilibrium")
    if route == "electrolyte_lle":
        _require_ion_containing_mixture(mixture, "reactive_electrolyte_lle")
        charges = _mixture_charges(mixture)
        feed, feed_diagnostics = _normalize_electrolyte_feed(
            mixture,
            z=z if z is not None else initial_x,
            solvent_feed=extra_phase_kwargs.get("solvent_feed"),
            salt_molality=extra_phase_kwargs.get("salt_molality"),
            options=solver_options,
        )
        _require_charge_neutral(feed, charges, "reactive_electrolyte_lle feed")
        electrolyte_formula_basis(mixture.species, charges, feed, salt_labels=tuple(feed_diagnostics.get("salt_molality", {})))
    else:
        if extra_phase_kwargs.get("solvent_feed") is not None or extra_phase_kwargs.get("salt_molality") is not None:
            raise InputError("solvent_feed and salt_molality require reactive_electrolyte_lle.")
        _reject_ion_containing_mixture(mixture)
        feed_source = z if z is not None else initial_x
        feed = _normalize_feed(feed_source, int(mixture.ncomp), solver_options.min_composition, "reactive_lle")

    _normalize_reactive_balances(species, balances, totals)
    reaction_defs = _normalize_reactions(species, reactions)
    _reaction_phase_stoichiometry_matrix(species, reaction_defs, route)
    _raise_native_ipopt_reactive_phase_required(route)


def _normalize_reactive_phase_route(
    mixture: Any,
    phase_kind: Any,
    phase_kwargs: Mapping[str, Any] | None,
) -> str:
    token = str("auto" if phase_kind is None else phase_kind).strip().lower()
    aliases = {
        "reactive_lle": "lle_flash",
        "reactive_lle_flash": "lle_flash",
        "lle_tp": "lle_flash",
        "reactive_electrolyte_lle": "electrolyte_lle",
        "reactive_electrolyte_lle_flash": "electrolyte_lle",
        "electrolyte_lle_tp": "electrolyte_lle",
    }
    token = aliases.get(token, token)
    if token == "auto":
        kwargs = dict(phase_kwargs or {})
        if kwargs.get("solvent_feed") is not None or kwargs.get("salt_molality") is not None:
            return "electrolyte_lle"
        charges = np.asarray(getattr(mixture, "parameters", {}).get("z", []), dtype=float).flatten()
        if charges.size == int(getattr(mixture, "ncomp", 0)) and np.any(np.abs(charges) > 1.0e-12):
            return "electrolyte_lle"
        return "lle_flash"
    if token not in {"lle_flash", "electrolyte_lle"}:
        raise InputError(
            "ReactivePhaseEquilibriumProblem production solves currently support phase_kind='lle_flash' "
            "or phase_kind='electrolyte_lle'. Use reactive_staged_equilibrium for explicit staged workflows."
        )
    return token


def _reject_reactive_phase_kwargs(phase_kwargs: Mapping[str, Any], route: str) -> None:
    allowed = {"solvent_feed", "salt_molality"}
    unsupported = sorted(key for key, value in phase_kwargs.items() if value is not None and key not in allowed)
    if unsupported:
        raise InputError(
            "reactive phase equilibrium does not support phase_kwargs key(s): {}.".format(", ".join(unsupported))
        )
    if route == "lle_flash" and (
        phase_kwargs.get("solvent_feed") is not None or phase_kwargs.get("salt_molality") is not None
    ):
        raise InputError("solvent_feed and salt_molality require reactive_electrolyte_lle.")


def _reactive_phase_option_pair(
    *,
    options: Any,
    phase_options: EquilibriumOptions | Mapping[str, Any] | None,
    normalize_reactive_options: Any,
) -> tuple[Any, EquilibriumOptions]:
    if isinstance(options, EquilibriumOptions) or isinstance(options, Mapping):
        if phase_options is not None:
            raise InputError("Use options or phase_options for reactive phase solver controls, not both.")
        return normalize_reactive_options(None), _normalize_options(options)
    reactive_options = normalize_reactive_options(options)
    if phase_options is not None:
        return reactive_options, _normalize_options(phase_options)
    return reactive_options, _equilibrium_options_from_reactive_options(reactive_options)


def _equilibrium_options_from_reactive_options(options: Any) -> EquilibriumOptions:
    return EquilibriumOptions(
        max_iterations=int(options.max_iterations),
        tolerance=float(options.tolerance),
        min_composition=float(options.min_mole_fraction),
        jacobian_backend=str(options.jacobian_backend),
        solver_backend="auto",
    )


def _reaction_phase_stoichiometry_matrix(
    species: list[str],
    reactions: list[Any],
    route: str,
) -> tuple[np.ndarray | None, str]:
    has_phase_terms = [reaction.phase_stoichiometry is not None for reaction in reactions]
    if not any(has_phase_terms):
        return None, "per_phase_same_stoichiometry"
    if not all(has_phase_terms):
        raise InputError("All reactions must use phase_stoichiometry when any reaction uses phase-tagged terms.")
    aliases = {
        "phase1": 0,
        "liq1": 0,
        "liquid1": 0,
        "phase_1": 0,
        "aqueous": 0,
        "aq": 0,
        "water": 0,
        "phase2": 1,
        "liq2": 1,
        "liquid2": 1,
        "phase_2": 1,
        "organic": 1,
        "org": 1,
        "solvent": 1,
    }
    if route == "lle_flash":
        aliases.update({"reactant_liquid": 0, "extract_liquid": 1})
    matrix = np.zeros((len(reactions), 2, len(species)), dtype=float)
    for reaction_index, reaction in enumerate(reactions):
        assert reaction.phase_stoichiometry is not None
        seen_phases: set[int] = set()
        for phase_label, coeffs in reaction.phase_stoichiometry.items():
            phase_key = str(phase_label).strip().lower()
            if phase_key not in aliases:
                supported = "', '".join(sorted(aliases))
                raise InputError(
                    f"Unknown phase label '{phase_label}' in reaction phase_stoichiometry; "
                    f"supported labels include '{supported}'."
                )
            phase_index = aliases[phase_key]
            seen_phases.add(phase_index)
            for label, coefficient in coeffs.items():
                matrix[reaction_index, phase_index, species.index(str(label))] = float(coefficient)
        if seen_phases != {0, 1}:
            raise InputError("phase-tagged reactions must include terms for both liquid phases.")
    return matrix, "phase_tagged_cross_phase"


def bubble_p(mixture: Any, *, T: float, x: Any, options: EquilibriumOptions | None = None) -> EquilibriumResult:
    """Solve a neutral bubble pressure at fixed liquid composition and temperature."""
    opts = _normalize_options(options)
    composition = _normalize_feed(x, mixture.ncomp, opts.min_composition, "bubble_p")
    _reject_ion_containing_mixture(mixture)
    temperature = _positive_scalar(T, "T", "bubble_p")
    return _native_neutral_fixed_temperature_pressure(
        mixture,
        T=temperature,
        composition=composition,
        options=opts,
        route_label="bubble_p",
        route_binding="_native_neutral_bubble_p_eos_route_result",
        problem_kind="neutral_bubble_p",
    )


def dew_p(mixture: Any, *, T: float, y: Any, options: EquilibriumOptions | None = None) -> EquilibriumResult:
    """Solve a neutral dew pressure at fixed vapor composition and temperature."""
    opts = _normalize_options(options)
    composition = _normalize_feed(y, mixture.ncomp, opts.min_composition, "dew_p")
    _reject_ion_containing_mixture(mixture)
    temperature = _positive_scalar(T, "T", "dew_p")
    return _native_neutral_fixed_temperature_pressure(
        mixture,
        T=temperature,
        composition=composition,
        options=opts,
        route_label="dew_p",
        route_binding="_native_neutral_dew_p_eos_route_result",
        problem_kind="neutral_dew_p",
    )


def bubble_t(mixture: Any, *, P: float, x: Any, options: EquilibriumOptions | None = None) -> EquilibriumResult:
    """Solve a neutral bubble temperature at fixed liquid composition and pressure."""
    opts = _normalize_options(options)
    _normalize_feed(x, mixture.ncomp, opts.min_composition, "bubble_t")
    _reject_ion_containing_mixture(mixture)
    _positive_scalar(P, "P", "bubble_t")
    _raise_native_ipopt_equilibrium_required("bubble_t")


def dew_t(mixture: Any, *, P: float, y: Any, options: EquilibriumOptions | None = None) -> EquilibriumResult:
    """Solve a neutral dew temperature at fixed vapor composition and pressure."""
    opts = _normalize_options(options)
    _normalize_feed(y, mixture.ncomp, opts.min_composition, "dew_t")
    _reject_ion_containing_mixture(mixture)
    _positive_scalar(P, "P", "dew_t")
    _raise_native_ipopt_equilibrium_required("dew_t")


def tp_flash(
    mixture: Any, *, T: float, P: float, z: Any, options: EquilibriumOptions | None = None
) -> EquilibriumResult:
    """Validate a neutral TP flash request and require the native Ipopt route."""
    opts = _normalize_options(options)
    feed = _normalize_feed(z, mixture.ncomp, opts.min_composition, "tp_flash")
    _reject_ion_containing_mixture(mixture)
    temperature = _positive_scalar(T, "T", "tp_flash")
    pressure = _positive_scalar(P, "P", "tp_flash")
    return _native_neutral_tp_flash(mixture, T=temperature, P=pressure, feed=feed, options=opts)


def lle_flash(
    mixture: Any,
    *,
    T: float,
    P: float,
    z: Any,
    options: EquilibriumOptions | None = None,
) -> EquilibriumResult:
    """Validate a neutral LLE flash request and require the native Ipopt route."""
    opts = _normalize_options(options)
    feed = _normalize_feed(z, mixture.ncomp, opts.min_composition, "lle_flash")
    _reject_ion_containing_mixture(mixture)
    temperature = _positive_scalar(T, "T", "lle_flash")
    pressure = _positive_scalar(P, "P", "lle_flash")
    return _native_neutral_lle_flash(mixture, T=temperature, P=pressure, feed=feed, options=opts)


def neutral_stability(
    mixture: Any,
    *,
    T: float,
    P: float,
    z: Any,
    options: EquilibriumOptions | None = None,
    parent_phase: Any = None,
    trial_phases: Any = None,
) -> StabilityResult:
    """Validate a neutral stability request and require the native Ipopt route."""
    opts = _normalize_options(options)
    _normalize_feed(z, mixture.ncomp, opts.min_composition, "stability")
    _reject_ion_containing_mixture(mixture)
    _positive_scalar(T, "T", "stability")
    _positive_scalar(P, "P", "stability")
    if parent_phase is not None:
        _normalize_parent_phases(parent_phase)
    if trial_phases is not None:
        _normalize_trial_phases(trial_phases)
    _raise_native_ipopt_stability_required("stability")


def electrolyte_stability(
    mixture: Any,
    *,
    T: float,
    P: float,
    z: Any = None,
    solvent_feed: Any = None,
    salt_molality: Any = None,
    options: EquilibriumOptions | None = None,
) -> StabilityResult:
    """Validate an electrolyte stability request and require the native Ipopt route."""
    opts = _normalize_options(options)
    _require_ion_containing_mixture(mixture, "electrolyte_stability")
    _positive_scalar(T, "T", "electrolyte_stability")
    _positive_scalar(P, "P", "electrolyte_stability")
    feed, feed_diagnostics = _normalize_electrolyte_feed(
        mixture,
        z=z,
        solvent_feed=solvent_feed,
        salt_molality=salt_molality,
        options=opts,
    )
    charges = _mixture_charges(mixture)
    _require_charge_neutral(feed, charges, "electrolyte_stability feed")
    electrolyte_formula_basis(mixture.species, charges, feed, salt_labels=tuple(feed_diagnostics.get("salt_molality", {})))
    _raise_native_ipopt_stability_required("electrolyte_stability")


def electrolyte_lle_flash_native(
    mixture: Any,
    *,
    T: float,
    P: float,
    z: Any = None,
    solvent_feed: Any = None,
    salt_molality: Any = None,
    options: EquilibriumOptions | None = None,
) -> EquilibriumResult:
    """Validate an electrolyte LLE request and require the native Ipopt route."""
    opts = _normalize_options(options)
    _require_ion_containing_mixture(mixture, "electrolyte_lle")
    _positive_scalar(T, "T", "electrolyte_lle")
    _positive_scalar(P, "P", "electrolyte_lle")
    feed, feed_diagnostics = _normalize_electrolyte_feed(
        mixture,
        z=z,
        solvent_feed=solvent_feed,
        salt_molality=salt_molality,
        options=opts,
    )
    charges = _mixture_charges(mixture)
    _require_charge_neutral(feed, charges, "electrolyte_lle feed")
    electrolyte_formula_basis(mixture.species, charges, feed, salt_labels=tuple(feed_diagnostics.get("salt_molality", {})))
    from . import _core

    route_tolerances = neutral_two_phase_eos_tolerances(P, opts)
    material_tolerance, pressure_tolerance, chemical_potential_tolerance, phase_distance_tolerance = route_tolerances
    charge_tolerance = min(opts.tolerance, 1.0e-8)
    route = _core._native_electrolyte_lle_eos_route_result(
        mixture._native,
        T,
        P,
        feed.tolist(),
        opts.max_iterations,
        opts.tolerance,
        _native_timeout_seconds(opts),
        material_tolerance,
        pressure_tolerance,
        charge_tolerance,
        chemical_potential_tolerance,
        phase_distance_tolerance,
    )
    if str(route.get("status", "")) == "ipopt_dependency_required":
        _raise_native_ipopt_lle_required("electrolyte_lle")
    return _accepted_native_neutral_two_phase_result(
        mixture,
        T=T,
        P=P,
        feed=feed,
        route=route,
        tolerances=route_tolerances,
        route_label="LLE",
        problem_kind="electrolyte_lle",
        phase_labels=("aq", "org"),
        route_family="electrolyte",
    )
