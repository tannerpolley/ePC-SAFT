"""Native-backed phase-equilibrium helpers and Python input adapters."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, replace
from numbers import Integral, Real
from typing import Any, Literal

import numpy as np

from ._types import InputError, SolutionError
from .equilibrium_core.electrolyte_basis import (
    electrolyte_feed_from_molality_inputs,
    electrolyte_formula_basis,
    normalize_salt_molality,
)
from .equilibrium_core.native_requests import (
    build_reactive_two_phase_eos_native_request,
    neutral_two_phase_eos_tolerances,
)
from .equilibrium_core.native_results import (
    native_route_solved_pressure,
    native_route_solved_temperature,
    native_route_summed_phase_amounts,
    native_route_diagnostics,
    neutral_two_phase_payload_to_result,
    raise_native_route_rejected,
)

_ASCANI_2022_REFERENCE = {
    "authors": "Ascani, Sadowski, and Held",
    "year": 2022,
    "title": "Calculation of Multiphase Equilibria Containing Mixed Solvents and Mixed Electrolytes",
    "doi": "10.1021/acs.jced.1c00866",
}

_ELECTROLYTE_LLE_TRACE_FLOOR = 1.0e-10
_ELECTROLYTE_LLE_TPDF_TOLERANCE = 1.0e-8


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


def _classify_problem_route(mixture: Any, kind: str) -> dict[str, str]:
    from .equilibrium_core.classify import classify_equilibrium_route

    return classify_equilibrium_route(mixture, kind)


def _with_route_diagnostics(result: Any, route: Mapping[str, str]) -> Any:
    diagnostics = dict(result.diagnostics)
    diagnostics["equilibrium_route"] = route["route"]
    diagnostics["route_reason"] = route["reason"]
    return replace(result, diagnostics=diagnostics)


def _raise_with_route_diagnostics(exc: SolutionError, route: Mapping[str, str]) -> None:
    diagnostics = dict(getattr(exc, "diagnostics", None) or {})
    diagnostics["equilibrium_route"] = route["route"]
    diagnostics["route_reason"] = route["reason"]
    raise SolutionError(exc.message, diagnostics) from exc


def _solve_problem_route(mixture: Any, kind: str, solve) -> Any:
    route = _classify_problem_route(mixture, kind)
    try:
        result = solve()
    except SolutionError as exc:
        _raise_with_route_diagnostics(exc, route)
    return _with_route_diagnostics(result, route)


def _reject_phase_controls_for_kind(
    kind: str,
    *,
    x_liq: Any = None,
    volatile_species: Any = None,
    vapor_species: Any = None,
    nonvolatile_species: Any = None,
) -> None:
    if x_liq is not None or volatile_species is not None or vapor_species is not None or nonvolatile_species is not None:
        raise InputError("x_liq and vapor species controls are only supported for kind='electrolyte_bubble_pressure'.")


def _reject_parent_controls_for_phase_equilibrium(parent_phase: Any, trial_phases: Any) -> None:
    if parent_phase is not None or trial_phases is not None:
        raise InputError("parent_phase and trial_phases are only supported for kind='stability'.")


def equilibrium_problem_from_request(
    *,
    kind: str,
    T: Any = None,
    P: Any = None,
    z: Any = None,
    x_liq: Any = None,
    volatile_species: Any = None,
    vapor_species: Any = None,
    nonvolatile_species: Any = None,
    solvent_feed: Any = None,
    salt_molality: Any = None,
    options: Any = None,
    parent_phase: Any = None,
    trial_phases: Any = None,
) -> EquilibriumProblem:
    """Convert a public string equilibrium request into a typed problem object."""
    token = str(kind).strip().lower()
    if token in {
        "bubble_p",
        "neutral_bubble_p",
        "bubble_t",
        "neutral_bubble_t",
        "dew_p",
        "neutral_dew_p",
        "dew_t",
        "neutral_dew_t",
    }:
        if solvent_feed is not None or salt_molality is not None:
            raise InputError("solvent_feed and salt_molality are only supported for kind='electrolyte_lle'.")
        if volatile_species is not None or vapor_species is not None or nonvolatile_species is not None:
            raise InputError("vapor species controls are only supported for kind='electrolyte_bubble_pressure'.")
        _reject_parent_controls_for_phase_equilibrium(parent_phase, trial_phases)
        if token in {"bubble_p", "neutral_bubble_p"}:
            if P is not None:
                raise InputError("P is solved by kind='bubble_p'.")
            composition = x_liq if x_liq is not None else z
            if composition is None:
                raise InputError("kind='bubble_p' requires x_liq or z as the liquid composition.")
            return BubblePoint(T=T, x=composition, options=options)
        if token in {"bubble_t", "neutral_bubble_t"}:
            if T is not None:
                raise InputError("T is solved by kind='bubble_t'.")
            composition = x_liq if x_liq is not None else z
            if composition is None:
                raise InputError("kind='bubble_t' requires x_liq or z as the liquid composition.")
            return BubblePoint(T=None, P=P, x=composition, options=options)
        if x_liq is not None:
            raise InputError("x_liq is only supported for bubble-point routes.")
        if z is None:
            raise InputError(f"kind='{kind}' requires z as the vapor composition.")
        if token in {"dew_p", "neutral_dew_p"}:
            if P is not None:
                raise InputError("P is solved by kind='dew_p'.")
            return DewPoint(T=T, y=z, options=options)
        if T is not None:
            raise InputError("T is solved by kind='dew_t'.")
        return DewPoint(P=P, y=z, options=options)

    if token in {"electrolyte_bubble_pressure", "electrolyte_bubble", "bubble_pressure"}:
        if P is not None:
            raise InputError("P is solved by kind='electrolyte_bubble_pressure'.")
        if solvent_feed is not None or salt_molality is not None:
            raise InputError("solvent_feed and salt_molality are not supported for kind='electrolyte_bubble_pressure'.")
        _reject_parent_controls_for_phase_equilibrium(parent_phase, trial_phases)
        return ElectrolyteBubblePoint(
            T=T,
            x_liq=x_liq,
            z=z,
            vapor_species=vapor_species,
            volatile_species=volatile_species,
            nonvolatile_species=nonvolatile_species,
            options=options,
        )

    if token == "tp_flash":
        if solvent_feed is not None or salt_molality is not None:
            raise InputError("solvent_feed and salt_molality are only supported for kind='electrolyte_lle'.")
        _reject_phase_controls_for_kind(
            token,
            x_liq=x_liq,
            volatile_species=volatile_species,
            vapor_species=vapor_species,
            nonvolatile_species=nonvolatile_species,
        )
        _reject_parent_controls_for_phase_equilibrium(parent_phase, trial_phases)
        return TPFlash(T=T, P=P, z=z, options=options)

    if token == "lle_flash":
        if solvent_feed is not None or salt_molality is not None:
            raise InputError("solvent_feed and salt_molality are only supported for kind='electrolyte_lle'.")
        _reject_phase_controls_for_kind(
            token,
            x_liq=x_liq,
            volatile_species=volatile_species,
            vapor_species=vapor_species,
            nonvolatile_species=nonvolatile_species,
        )
        _reject_parent_controls_for_phase_equilibrium(parent_phase, trial_phases)
        return LLEProblem(T=T, P=P, z=z, options=options)

    if token in {"electrolyte_lle", "electrolyte_lle_flash"}:
        _reject_phase_controls_for_kind(
            token,
            x_liq=x_liq,
            volatile_species=volatile_species,
            vapor_species=vapor_species,
            nonvolatile_species=nonvolatile_species,
        )
        _reject_parent_controls_for_phase_equilibrium(parent_phase, trial_phases)
        return ElectrolyteLLEProblem(
            T=T,
            P=P,
            z=z,
            solvent_feed=solvent_feed,
            salt_molality=salt_molality,
            options=options,
        )

    if token == "electrolyte_stability":
        _reject_phase_controls_for_kind(
            token,
            x_liq=x_liq,
            volatile_species=volatile_species,
            vapor_species=vapor_species,
            nonvolatile_species=nonvolatile_species,
        )
        if parent_phase is not None or trial_phases is not None:
            raise InputError("parent_phase and trial_phases are not supported for kind='electrolyte_stability'.")
        return StabilityAnalysis(
            T=T,
            P=P,
            z=z,
            options=options,
            solvent_feed=solvent_feed,
            salt_molality=salt_molality,
            electrolyte=True,
        )

    if token == "stability":
        if solvent_feed is not None or salt_molality is not None:
            raise InputError("solvent_feed and salt_molality are only supported for kind='electrolyte_lle'.")
        _reject_phase_controls_for_kind(
            token,
            x_liq=x_liq,
            volatile_species=volatile_species,
            vapor_species=vapor_species,
            nonvolatile_species=nonvolatile_species,
        )
        return StabilityAnalysis(
            T=T,
            P=P,
            z=z,
            options=options,
            parent_phase=parent_phase,
            trial_phases=trial_phases,
        )

    raise InputError(
        "Only kind='tp_flash', kind='auto', kind='bubble_p', kind='bubble_t', kind='dew_p', kind='dew_t', kind='lle_flash', kind='electrolyte_lle', kind='electrolyte_bubble_pressure', kind='electrolyte_stability', kind='stability', kind='chemical_equilibrium', kind='reactive_staged_equilibrium', kind='reactive_lle', kind='reactive_electrolyte_lle', kind='reactive_stability', or kind='reactive_electrolyte_bubble_pressure' is supported by equilibrium."
    )


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
        return _solve_problem_route(
            mixture,
            "tp_flash",
            lambda: mixture.flash_tp(T=self.T, P=self.P, z=self.z, options=self.options),
        )


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
            return _solve_problem_route(
                mixture,
                "electrolyte_stability",
                lambda: mixture.electrolyte_stability_tp(
                    T=self.T,
                    P=self.P,
                    z=self.z,
                    solvent_feed=self.solvent_feed,
                    salt_molality=self.salt_molality,
                    options=self.options,
                ),
            )
        return _solve_problem_route(
            mixture,
            "stability",
            lambda: mixture.stability_tp(
                T=self.T,
                P=self.P,
                z=self.z,
                options=self.options,
                parent_phase=self.parent_phase,
                trial_phases=self.trial_phases,
            ),
        )


@dataclass(frozen=True, slots=True)
class BubblePoint(EquilibriumProblem):
    """Neutral bubble-point problem."""

    T: float | None
    x: Any
    P: float | None = None
    options: EquilibriumOptions | None = None

    def solve(self, mixture):
        if (self.T is None) == (self.P is None):
            raise InputError("BubblePoint requires exactly one of T or P.")
        if self.T is not None:
            return _solve_problem_route(
                mixture,
                "bubble_p",
                lambda: mixture.bubble_p(T=self.T, x=self.x, options=self.options),
            )
        return _solve_problem_route(
            mixture,
            "bubble_t",
            lambda: mixture.bubble_t(P=self.P, x=self.x, options=self.options),
        )


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
            return _solve_problem_route(
                mixture,
                "dew_p",
                lambda: mixture.dew_p(T=self.T, y=self.y, options=self.options),
            )
        return _solve_problem_route(
            mixture,
            "dew_t",
            lambda: mixture.dew_t(P=self.P, y=self.y, options=self.options),
        )


@dataclass(frozen=True, slots=True)
class LLEProblem(EquilibriumProblem):
    """Neutral liquid-liquid equilibrium problem."""

    T: float
    P: float
    z: Any
    options: EquilibriumOptions | None = None

    def solve(self, mixture):
        return _solve_problem_route(
            mixture,
            "lle_flash",
            lambda: mixture.lle_tp(
                T=self.T,
                P=self.P,
                z=self.z,
                options=self.options,
            ),
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
        return _solve_problem_route(
            mixture,
            "electrolyte_lle",
            lambda: mixture.electrolyte_lle_tp(
                T=self.T,
                P=self.P,
                z=self.z,
                solvent_feed=self.solvent_feed,
                salt_molality=self.salt_molality,
                options=self.options,
            ),
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
        return _solve_problem_route(
            mixture,
            "electrolyte_bubble_pressure",
            lambda: mixture.electrolyte_bubble_p(
                T=self.T,
                x_liq=self.x_liq,
                z=self.z,
                vapor_species=self.vapor_species,
                volatile_species=self.volatile_species,
                nonvolatile_species=self.nonvolatile_species,
                options=self.options,
            ),
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
    feed_diagnostics: Mapping[str, Any] | None = None,
    electrolyte_basis_diagnostics: Mapping[str, Any] | None = None,
    options: EquilibriumOptions | None = None,
) -> EquilibriumResult:
    from . import _core

    material_tolerance, pressure_tolerance, chemical_potential_tolerance, phase_distance_tolerance = tolerances
    if not bool(route.get("accepted", False)):
        raise_native_route_rejected(route, f"Native {route_family} {route_label} route was rejected.")

    if route_family == "electrolyte":
        postsolve = route.get("postsolve", {})
        if not isinstance(postsolve, Mapping) or not bool(postsolve.get("accepted", False)):
            raise_native_route_rejected(route, f"Native electrolyte {route_label} postsolve was rejected.")
        phase_compositions = np.asarray(postsolve.get("phase_compositions", ()), dtype=float)
        phase_amount_totals = np.asarray(postsolve.get("phase_amount_totals", ()), dtype=float).reshape(-1)
        phase_volumes = np.asarray(postsolve.get("phase_volumes", ()), dtype=float).reshape(-1)
        phase_ln_phi_payload = postsolve.get("phase_ln_fugacity_coefficients", ())
        phase_ln_phi = np.asarray(phase_ln_phi_payload, dtype=float) if phase_ln_phi_payload else np.asarray(())
        if (
            phase_compositions.ndim != 2
            or phase_compositions.shape[0] != 2
            or phase_compositions.shape[1] != int(mixture.ncomp)
            or phase_amount_totals.size != 2
            or phase_volumes.size != 2
            or not np.all(np.isfinite(phase_compositions))
            or not np.all(np.isfinite(phase_amount_totals))
            or not np.all(np.isfinite(phase_volumes))
            or np.any(phase_amount_totals <= 0.0)
            or np.any(phase_volumes <= 0.0)
        ):
            raise SolutionError("Native electrolyte LLE route accepted without a valid retained phase split.")
        total_amount = float(np.sum(phase_amount_totals))
        phase_order = [0, 1]
        if phase_labels == ("aq", "org"):
            charges = np.abs(_mixture_charges(mixture))
            if charges.size == phase_compositions.shape[1] and np.any(charges > 0.0):
                ion_fractions = phase_compositions @ (charges > 0.0).astype(float)
                aqueous_index = int(np.argmax(ion_fractions))
                phase_order = [aqueous_index, 1 - aqueous_index]
        phases = []
        for output_index, label in enumerate(phase_labels):
            index = phase_order[output_index]
            density = float(phase_amount_totals[index] / phase_volumes[index])
            ln_phi = None
            if phase_ln_phi.ndim == 2 and phase_ln_phi.shape == phase_compositions.shape:
                ln_phi = phase_ln_phi[index]
            phases.append(
                EquilibriumPhase(
                    label,
                    composition=phase_compositions[index],
                    density=density,
                    temperature=T,
                    pressure=P,
                    phase_fraction=float(phase_amount_totals[index] / total_amount),
                    ln_fugacity_coefficient=ln_phi,
                    diagnostics={
                        "amount_total": float(phase_amount_totals[index]),
                        "volume": float(phase_volumes[index]),
                    },
                )
            )
        density_errors = []
        for phase in phases:
            if float(phase.density) < 1000.0:
                raise SolutionError(
                    "Native electrolyte LLE route accepted a non-liquid density root.",
                    {"phase": phase.label, "density_mol_m3": float(phase.density)},
                )
            state = mixture.state(T, phase.composition, P=P, phase="liq", rho_guess=phase.density)
            recomputed_density = float(state.molar_density())
            rel_error = abs(recomputed_density - float(phase.density)) / max(abs(float(phase.density)), 1.0)
            density_errors.append(
                {
                    "phase": phase.label,
                    "reported_density_mol_m3": float(phase.density),
                    "recomputed_density_mol_m3": recomputed_density,
                    "relative_error": rel_error,
                }
            )
            if rel_error > 1.0e-8:
                raise SolutionError(
                    "Native electrolyte LLE liquid density recomputation failed.",
                    {"phase": phase.label, "relative_error": rel_error},
                )
        if float(postsolve.get("gibbs_delta", 0.0)) >= -1.0e-8:
            raise SolutionError(
                "Native electrolyte LLE route did not lower the retained Gibbs objective.",
                {"gibbs_delta": float(postsolve.get("gibbs_delta", 0.0))},
            )
        diagnostics = native_route_diagnostics(
            route,
            defaults={
                "solver_backend": "ipopt",
                "problem_name": "electrolyte_lle_eos",
                "adapter_kind": "native_tnlp_adapter",
                "exact_gradient_required": True,
                "exact_jacobian_required": True,
                "gradient_approximation": "exact",
                "jacobian_approximation": "exact",
                "hessian_approximation": "limited-memory",
            },
        )
        diagnostics["feed_basis"] = _json_like(feed_diagnostics or {})
        basis_payload = dict(electrolyte_basis_diagnostics or {})
        diagnostics["electrolyte_basis"] = _json_like(basis_payload)
        if basis_payload:
            diagnostics["variable_model"] = str(basis_payload.get("variable_model", ""))
            diagnostics["salt_pairs"] = _json_like(basis_payload.get("salt_pairs", ()))
            diagnostics["basis_rank"] = int(basis_payload.get("basis_rank", 0))
            diagnostics["formula_feed"] = _json_like(basis_payload.get("formula_feed", ()))
        diagnostics["species_labels"] = list(getattr(mixture, "species", ()))
        diagnostics["charge_vector"] = _mixture_charges(mixture).tolist()
        diagnostics["density_recompute_relative_errors"] = density_errors
        _add_electrolyte_lle_residual_diagnostics(
            diagnostics,
            phases=tuple(phases),
            basis=basis_payload,
            trace_floor=_ELECTROLYTE_LLE_TRACE_FLOOR,
        )
        if float(diagnostics["material_balance_norm"]) > 1.0e-8:
            raise SolutionError("Native electrolyte LLE material-balance gate failed.", diagnostics)
        if float(diagnostics["charge_balance_norm"]) > 1.0e-8:
            raise SolutionError("Native electrolyte LLE charge-balance gate failed.", diagnostics)
        diagnostics["neutral_fugacity_residual_tolerance"] = 1.0e-7
        diagnostics["salt_pair_fugacity_residual_tolerance"] = 1.0e-7
        if float(diagnostics["neutral_fugacity_residual_norm"]) > 1.0e-7:
            raise SolutionError("Native electrolyte LLE neutral fugacity gate failed.", diagnostics)
        if float(diagnostics["salt_pair_fugacity_residual_norm"]) > 1.0e-7:
            raise SolutionError("Native electrolyte LLE salt-pair fugacity gate failed.", diagnostics)
        tpdf_certificate = _electrolyte_lle_tpdf_certificate(
            mixture,
            T=T,
            P=P,
            feed=feed,
            phases=tuple(phases),
            options=_normalize_options(options),
            stability_tolerance=_ELECTROLYTE_LLE_TPDF_TOLERANCE,
        )
        diagnostics["tpdf_stability"] = tpdf_certificate
        if not bool(tpdf_certificate.get("accepted", False)):
            raise SolutionError("Native electrolyte LLE hard TPDF stability certification failed.", diagnostics)
        return EquilibriumResult(
            backend="native_equilibrium_nlp",
            problem_kind=problem_kind,
            phases=tuple(phases),
            stable=False,
            split_detected=float(postsolve.get("phase_distance", 0.0)) >= phase_distance_tolerance,
            diagnostics=diagnostics,
        )

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


def _accepted_native_reactive_two_phase_result(
    mixture: Any,
    *,
    T: float,
    P: float,
    route: Mapping[str, Any],
    tolerances: tuple[float, float, float, float],
    route_label: str,
    problem_kind: str,
    phase_labels: tuple[str, str],
) -> EquilibriumResult:
    from . import _core

    if not bool(route.get("accepted", False)):
        raise_native_route_rejected(route, f"Native reactive {route_label} route was rejected.")

    material_tolerance, pressure_tolerance, chemical_potential_tolerance, phase_distance_tolerance = tolerances
    postsolve = route.get("postsolve", {})
    if isinstance(postsolve, Mapping) and str(postsolve.get("density_backend", "")) == "liquid_pressure_root":
        phase_compositions = np.asarray(postsolve.get("phase_compositions", ()), dtype=float)
        phase_amount_totals = np.asarray(postsolve.get("phase_amount_totals", ()), dtype=float).reshape(-1)
        phase_volumes = np.asarray(postsolve.get("phase_volumes", ()), dtype=float).reshape(-1)
        phase_ln_phi_payload = postsolve.get("phase_ln_fugacity_coefficients", ())
        phase_ln_phi = np.asarray(phase_ln_phi_payload, dtype=float) if phase_ln_phi_payload else np.asarray(())
        if (
            phase_compositions.ndim != 2
            or phase_compositions.shape[0] != 2
            or phase_compositions.shape[1] != int(mixture.ncomp)
            or phase_amount_totals.size != 2
            or phase_volumes.size != 2
            or not np.all(np.isfinite(phase_compositions))
            or not np.all(np.isfinite(phase_amount_totals))
            or not np.all(np.isfinite(phase_volumes))
            or np.any(phase_amount_totals <= 0.0)
            or np.any(phase_volumes <= 0.0)
        ):
            raise SolutionError(f"Native reactive {route_label} route accepted without a valid liquid-root phase split.")
        total_amount = float(np.sum(phase_amount_totals))
        phases = []
        for index, label in enumerate(phase_labels):
            ln_phi = None
            if phase_ln_phi.ndim == 2 and phase_ln_phi.shape == phase_compositions.shape:
                ln_phi = phase_ln_phi[index]
            phases.append(
                EquilibriumPhase(
                    label,
                    composition=phase_compositions[index],
                    density=float(phase_amount_totals[index] / phase_volumes[index]),
                    temperature=T,
                    pressure=P,
                    phase_fraction=float(phase_amount_totals[index] / total_amount),
                    ln_fugacity_coefficient=ln_phi,
                    diagnostics={
                        "amount_total": float(phase_amount_totals[index]),
                        "volume": float(phase_volumes[index]),
                    },
                )
            )
        diagnostics = native_route_diagnostics(
            route,
            route_status_key="reactive_route_status",
            defaults={
                "solver_backend": "ipopt",
                "problem_name": "reactive_liquid_root_eos",
            },
        )
        return EquilibriumResult(
            backend="native_equilibrium_nlp",
            problem_kind=problem_kind,
            phases=tuple(phases),
            stable=False,
            split_detected=float(postsolve.get("phase_distance", 0.0)) >= phase_distance_tolerance,
            diagnostics=diagnostics,
        )

    phase_amounts = np.asarray(route["phase_amounts"], dtype=float)
    if phase_amounts.ndim != 2 or phase_amounts.shape[1] != int(mixture.ncomp):
        raise SolutionError(f"Native reactive {route_label} route phase amounts had an invalid shape.")
    reacted_feed = np.sum(phase_amounts, axis=0)
    if not np.all(np.isfinite(reacted_feed)) or np.any(reacted_feed <= 0.0):
        raise SolutionError(f"Native reactive {route_label} route did not return positive species totals.")

    result_payload = _core._native_neutral_two_phase_eos_result(
        mixture._native,
        T,
        P,
        route["phase_amounts"],
        route["phase_volumes"],
        reacted_feed.tolist(),
        material_tolerance,
        pressure_tolerance,
        chemical_potential_tolerance,
        phase_distance_tolerance,
    )
    result = neutral_two_phase_payload_to_result(
        result_payload,
        problem_kind=problem_kind,
        phase_labels=phase_labels,
    )
    diagnostics = dict(result.diagnostics)
    diagnostics["reactive_postsolve"] = _json_like(route.get("postsolve", {}))
    diagnostics["reactive_route_status"] = str(route.get("status", ""))
    diagnostics["reaction_count"] = int(route.get("reaction_count", 0))
    diagnostics["balance_row_count"] = int(route.get("balance_row_count", 0))
    return EquilibriumResult(
        backend=result.backend,
        problem_kind=result.problem_kind,
        phases=result.phases,
        stable=result.stable,
        split_detected=result.split_detected,
        diagnostics=diagnostics,
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


def _native_neutral_fixed_pressure_temperature(
    mixture: Any,
    *,
    P: float,
    composition: np.ndarray,
    options: EquilibriumOptions,
    route_label: str,
    route_binding: str,
    problem_kind: str,
) -> EquilibriumResult:
    from . import _core

    route_tolerances = (
        options.tolerance,
        max(P * options.tolerance, options.tolerance),
        options.tolerance,
        max(10.0 * options.min_composition, 1.0e-8),
    )
    route = getattr(_core, route_binding)(
        mixture._native,
        P,
        composition.tolist(),
        options.max_iterations,
        options.tolerance,
        _native_timeout_seconds(options),
        *route_tolerances,
    )
    if str(route.get("status", "")) == "ipopt_dependency_required":
        _raise_native_ipopt_equilibrium_required(route_label)

    temperature = native_route_solved_temperature(route, route_label) if bool(route.get("accepted", False)) else 1.0
    feed = (
        native_route_summed_phase_amounts(route, mixture.ncomp, route_label)
        if bool(route.get("accepted", False))
        else composition
    )
    return _accepted_native_neutral_two_phase_result(
        mixture,
        T=temperature,
        P=P,
        feed=feed,
        route=route,
        tolerances=neutral_two_phase_eos_tolerances(P, options),
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


def _native_reactive_two_phase_flash(
    mixture: Any,
    *,
    T: float,
    P: float,
    feed: np.ndarray,
    balance_matrix: np.ndarray,
    total_vector: np.ndarray,
    species: list[str],
    reactions: list[Any],
    reaction_phase_stoichiometry: np.ndarray | None,
    options: EquilibriumOptions,
    route_binding: str,
    required_route: str,
    problem_kind: str,
    phase_labels: tuple[str, str],
    phase_models: Mapping[str, Any] | None = None,
) -> EquilibriumResult:
    from . import _core

    request = build_reactive_two_phase_eos_native_request(
        T=T,
        P=P,
        feed=feed,
        balance_matrix=balance_matrix,
        total_vector=total_vector,
        species=species,
        reactions=reactions,
        reaction_phase_stoichiometry=reaction_phase_stoichiometry,
    )
    material_tolerance, _, chemical_potential_tolerance, phase_distance_tolerance = neutral_two_phase_eos_tolerances(
        P,
        options,
    )
    if phase_models is not None:
        if route_binding != "_native_reactive_electrolyte_lle_eos_route_result":
            raise InputError("phase_models requires reactive_electrolyte_lle.")
        route = _core._native_reactive_electrolyte_lle_phase_model_eos_route_result(
            mixture._native,
            phase_models["aq"]._native,
            phase_models["org"]._native,
            phase_models["aq_indices"],
            phase_models["org_indices"],
            request["T"],
            request["P"],
            request["feed_amounts"],
            request["balance_rows"],
            request["balance_matrix"],
            request["total_vector"],
            request["reaction_rows"],
            request["reaction_stoichiometry"],
            request["log_equilibrium_constants"],
            options.max_iterations,
            options.tolerance,
            _native_timeout_seconds(options),
            material_tolerance,
            chemical_potential_tolerance,
            chemical_potential_tolerance,
            phase_distance_tolerance,
            options.min_composition,
            request["reaction_standard_states"],
            request["reaction_phase_stoichiometry"],
        )
    else:
        route = getattr(_core, route_binding)(
            mixture._native,
            request["T"],
            request["P"],
            request["feed_amounts"],
            request["balance_rows"],
            request["balance_matrix"],
            request["total_vector"],
            request["reaction_rows"],
            request["reaction_stoichiometry"],
            request["log_equilibrium_constants"],
            options.max_iterations,
            options.tolerance,
            _native_timeout_seconds(options),
            material_tolerance,
            chemical_potential_tolerance,
            chemical_potential_tolerance,
            phase_distance_tolerance,
            options.min_composition,
            request["reaction_standard_states"],
            request["reaction_phase_stoichiometry"],
        )
    if str(route.get("status", "")) == "ipopt_dependency_required":
        _raise_native_ipopt_reactive_phase_required(required_route)
    return _accepted_native_reactive_two_phase_result(
        mixture,
        T=T,
        P=P,
        route=route,
        tolerances=neutral_two_phase_eos_tolerances(P, options),
        route_label="LLE",
        problem_kind=problem_kind,
        phase_labels=phase_labels,
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


def _native_stability_trial_from_route(route: Mapping[str, Any], stability_tolerance: float) -> StabilityTrial:
    composition = np.asarray(route.get("trial_composition", route.get("variables", [])), dtype=float).flatten()
    if composition.size == 0 or not np.all(np.isfinite(composition)):
        raise SolutionError("Native stability route did not return a finite trial composition.", dict(route))
    min_tpd = float(route.get("min_tpd", route.get("objective", np.nan)))
    if not np.isfinite(min_tpd):
        raise SolutionError("Native stability route did not return a finite TPD objective.", dict(route))
    diagnostics = {
        "route_status": str(route.get("status", "")),
        "solver_status": str(route.get("solver_status", "")),
        "application_status": str(route.get("application_status", "")),
        "derivative_backend": str(route.get("derivative_backend", "")),
        "constraints": _json_like(route.get("constraints", [])),
    }
    return StabilityTrial(
        parent_phase=str(route.get("parent_phase", "")),
        trial_phase=str(route.get("trial_phase", "")),
        seed_name=str(route.get("seed_name", "canonical_shifted_feed")),
        composition=composition,
        tpd=min_tpd,
        iterations=0,
        converged=bool(route.get("solver_accepted", False)),
        unstable=min_tpd < -abs(stability_tolerance),
        diagnostics=diagnostics,
    )


def _native_neutral_stability(
    mixture: Any,
    *,
    T: float,
    P: float,
    feed: np.ndarray,
    options: EquilibriumOptions,
    parent_phases: tuple[str, ...],
    trial_phases: tuple[str, ...],
) -> StabilityResult:
    from . import _core

    route_payloads: list[Mapping[str, Any]] = []
    for parent in parent_phases:
        for trial in trial_phases:
            route = _core._native_neutral_stability_tpd_route_result(
                mixture._native,
                T,
                P,
                feed.tolist(),
                parent,
                trial,
                options.max_iterations,
                options.tolerance,
                _native_timeout_seconds(options),
                options.tolerance,
            )
            if str(route.get("status", "")) == "ipopt_dependency_required":
                _raise_native_ipopt_stability_required("stability")
            route_payloads.append(route)

    rejected = [route for route in route_payloads if not bool(route.get("accepted", False))]
    if rejected:
        diagnostics = {
            "route_statuses": [str(route.get("status", "")) for route in route_payloads],
            "solver_statuses": [str(route.get("solver_status", "")) for route in route_payloads],
            "parent_trial_routes": [
                [str(route.get("parent_phase", "")), str(route.get("trial_phase", ""))] for route in route_payloads
            ],
        }
        raise SolutionError("Native neutral stability route was rejected.", diagnostics)

    trials = tuple(_native_stability_trial_from_route(route, options.tolerance) for route in route_payloads)
    if not trials:
        raise SolutionError("Native neutral stability produced no trial routes.")
    tpd_values = np.asarray([trial.tpd for trial in trials], dtype=float)
    selected_index = int(np.argmin(tpd_values))
    selected = trials[selected_index]
    stable = selected.tpd >= -abs(options.tolerance)
    return StabilityResult(
        backend="native_equilibrium_nlp",
        problem_kind="neutral_stability",
        stable=stable,
        min_tpd=selected.tpd,
        parent_phase=selected.parent_phase,
        trial_phase=selected.trial_phase,
        trial_composition=selected.composition,
        trials=trials,
        diagnostics={
            "backend": "ipopt",
            "derivative_backend": "cppad_implicit",
            "route_count": len(trials),
            "selected_trial_index": selected_index,
            "stability_tolerance": options.tolerance,
            "parent_phases": list(parent_phases),
            "trial_phases": list(trial_phases),
        },
    )


def _native_electrolyte_stability(
    mixture: Any,
    *,
    T: float,
    P: float,
    feed: np.ndarray,
    feed_diagnostics: Mapping[str, Any],
    options: EquilibriumOptions,
) -> StabilityResult:
    from . import _core

    route = _core._native_electrolyte_stability_tpd_route_result(
        mixture._native,
        T,
        P,
        feed.tolist(),
        options.max_iterations,
        options.tolerance,
        _native_timeout_seconds(options),
        options.tolerance,
    )
    if str(route.get("status", "")) == "ipopt_dependency_required":
        _raise_native_ipopt_stability_required("electrolyte_stability")
    if not bool(route.get("accepted", False)):
        diagnostics = native_route_diagnostics(route)
        diagnostics["parent_phase"] = str(route.get("parent_phase", ""))
        diagnostics["trial_phase"] = str(route.get("trial_phase", ""))
        raise SolutionError("Native electrolyte stability route was rejected.", diagnostics)

    trial = _native_stability_trial_from_route(route, options.tolerance)
    stable = trial.tpd >= -abs(options.tolerance)
    return StabilityResult(
        backend="native_equilibrium_nlp",
        problem_kind="electrolyte_stability",
        stable=stable,
        min_tpd=trial.tpd,
        parent_phase=trial.parent_phase,
        trial_phase=trial.trial_phase,
        trial_composition=trial.composition,
        trials=(trial,),
        diagnostics={
            "backend": "ipopt",
            "derivative_backend": "cppad_implicit",
            "route_count": 1,
            "selected_trial_index": 0,
            "stability_tolerance": options.tolerance,
            "feed_basis": _json_like(feed_diagnostics),
        },
    )


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


def _phase_log_fugacity_vector(phase: EquilibriumPhase) -> np.ndarray:
    x = np.asarray(phase.composition, dtype=float)
    ln_phi = np.asarray(phase.ln_fugacity_coefficient, dtype=float)
    if ln_phi.shape != x.shape:
        raise SolutionError("Electrolyte LLE phase is missing retained ln fugacity coefficients.")
    return np.log(np.maximum(x, 1.0e-300)) + ln_phi


def _add_electrolyte_lle_residual_diagnostics(
    diagnostics: dict[str, Any],
    *,
    phases: tuple[EquilibriumPhase, EquilibriumPhase],
    basis: Mapping[str, Any],
    trace_floor: float,
) -> None:
    if len(phases) != 2:
        raise SolutionError("Electrolyte LLE residual diagnostics require exactly two phases.")
    left, right = phases
    left_ln_f = _phase_log_fugacity_vector(left)
    right_ln_f = _phase_log_fugacity_vector(right)
    neutral_indices = tuple(int(index) for index in basis.get("neutral_indices", ()))
    salt_pairs = tuple(dict(pair) for pair in basis.get("salt_pairs", ()))
    species = list(diagnostics.get("species_labels", ()))

    neutral_raw: dict[str, float] = {}
    neutral_weighted: dict[str, float] = {}
    for index in neutral_indices:
        label = species[index] if 0 <= index < len(species) else f"species_{index}"
        raw = float(left_ln_f[index] - right_ln_f[index])
        weight = _trace_residual_weight(left.composition[index], right.composition[index], trace_floor=trace_floor)
        neutral_raw[label] = raw
        neutral_weighted[label] = raw * weight

    salt_raw: dict[str, float] = {}
    salt_weighted: dict[str, float] = {}
    for pair in salt_pairs:
        cation = int(pair["cation"])
        anion = int(pair["anion"])
        c_stoich = float(pair.get("cation_stoich", 1.0))
        a_stoich = float(pair.get("anion_stoich", 1.0))
        left_value = c_stoich * left_ln_f[cation] + a_stoich * left_ln_f[anion]
        right_value = c_stoich * right_ln_f[cation] + a_stoich * right_ln_f[anion]
        raw = float(left_value - right_value)
        weight = _trace_residual_weight(
            left.composition[cation],
            right.composition[cation],
            left.composition[anion],
            right.composition[anion],
            trace_floor=trace_floor,
        )
        label = str(pair.get("label", f"{cation}:{anion}"))
        salt_raw[label] = raw
        salt_weighted[label] = raw * weight

    diagnostics["trace_floor"] = float(trace_floor)
    diagnostics["neutral_log_fugacity_residuals_raw"] = neutral_raw
    diagnostics["neutral_log_fugacity_residuals_weighted"] = neutral_weighted
    diagnostics["salt_pair_log_fugacity_residuals_raw"] = salt_raw
    diagnostics["salt_pair_log_fugacity_residuals_weighted"] = salt_weighted
    neutral_norm = max((abs(value) for value in neutral_weighted.values()), default=0.0)
    salt_norm = max((abs(value) for value in salt_weighted.values()), default=0.0)
    diagnostics["neutral_fugacity_residual_norm"] = float(neutral_norm)
    diagnostics["salt_pair_fugacity_residual_norm"] = float(salt_norm)
    diagnostics["mean_ionic_fugacity_residual_norm"] = float(salt_norm)


def _trace_residual_weight(*values: float, trace_floor: float) -> float:
    if trace_floor <= 0.0:
        return 1.0
    minimum = min(abs(float(value)) for value in values)
    return min(1.0, minimum / trace_floor)


def _electrolyte_lle_tpdf_certificate(
    mixture: Any,
    *,
    T: float,
    P: float,
    feed: np.ndarray,
    phases: tuple[EquilibriumPhase, EquilibriumPhase],
    options: EquilibriumOptions,
    stability_tolerance: float,
) -> dict[str, Any]:
    from . import _core

    charges = _mixture_charges(mixture)
    seeds = _electrolyte_lle_tpdf_seed_set(feed=feed, phases=phases, charges=charges)
    trial_records: list[dict[str, Any]] = []
    parent_specs: list[tuple[str, np.ndarray, str]] = [("feed", np.asarray(feed, dtype=float), "unstable")]
    parent_specs.extend((phase.label, np.asarray(phase.composition, dtype=float), "stable") for phase in phases)
    max_iterations = min(int(options.max_iterations), 120)
    for parent_label, parent_composition, expected in parent_specs:
        for seed_name, trial_seed in seeds:
            restart_count = 0
            try:
                route = _core._native_electrolyte_stability_tpd_route_result(
                    mixture._native,
                    T,
                    P,
                    parent_composition.tolist(),
                    max_iterations,
                    options.tolerance,
                    _native_timeout_seconds(options),
                    stability_tolerance,
                    np.asarray(trial_seed, dtype=float).tolist(),
                )
            except Exception as exc:
                trial_records.append(
                    {
                        "parent_phase": parent_label,
                        "seed_name": seed_name,
                        "status": "solver_exception",
                        "min_tpd": None,
                        "tolerance": stability_tolerance,
                        "failure_reason": str(exc),
                    }
                )
                continue
            if not bool(route.get("accepted", False)):
                restart = np.asarray(route.get("variables", ()), dtype=float).reshape(-1)
                if restart.size == feed.size and np.all(np.isfinite(restart)) and np.all(restart > 0.0):
                    try:
                        route = _core._native_electrolyte_stability_tpd_route_result(
                            mixture._native,
                            T,
                            P,
                            parent_composition.tolist(),
                            max_iterations,
                            options.tolerance,
                            _native_timeout_seconds(options),
                            stability_tolerance,
                            restart.tolist(),
                        )
                        restart_count = 1
                    except Exception:
                        pass
            min_tpd = float(route.get("min_tpd", route.get("objective", np.nan)))
            accepted = bool(route.get("accepted", False)) and np.isfinite(min_tpd)
            trial_records.append(
                {
                    "parent_phase": parent_label,
                    "expected": expected,
                    "seed_name": seed_name,
                    "status": str(route.get("status", "")),
                    "solver_status": str(route.get("solver_status", "")),
                    "application_status": str(route.get("application_status", "")),
                    "min_tpd": None if not np.isfinite(min_tpd) else min_tpd,
                    "tolerance": stability_tolerance,
                    "trial_composition": _json_like(route.get("trial_composition", ())),
                    "initial_composition": _json_like(route.get("initial_composition", trial_seed)),
                    "restart_count": restart_count,
                    "failure_reason": "" if accepted else "solver_rejected",
                }
            )
            if parent_label == "feed" and accepted and min_tpd < -abs(stability_tolerance):
                break
    feed_trials = [row for row in trial_records if row["parent_phase"] == "feed"]
    phase_trials = [row for row in trial_records if row["parent_phase"] != "feed"]
    accepted_feed_trials = [
        row for row in feed_trials if row.get("status") == "accepted" and row.get("min_tpd") is not None
    ]
    accepted_phase_trials = [
        row for row in phase_trials if row.get("status") == "accepted" and row.get("min_tpd") is not None
    ]
    rejected_feed_trials = len(accepted_feed_trials) != len(feed_trials)
    rejected_phase_trials = len(accepted_phase_trials) != len(phase_trials)
    feed_unstable = any(float(row["min_tpd"]) < -abs(stability_tolerance) for row in accepted_feed_trials)
    final_phases_stable = bool(phase_trials) and not rejected_phase_trials and all(
        float(row["min_tpd"]) >= -abs(stability_tolerance) for row in accepted_phase_trials
    )
    if (not feed_unstable and rejected_feed_trials) or rejected_phase_trials:
        return {
            "accepted": False,
            "status": "blocked_solver",
            "tolerance": stability_tolerance,
            "feed_unstable": feed_unstable,
            "final_phases_stable": final_phases_stable,
            "trials": trial_records,
        }
    accepted = feed_unstable and final_phases_stable
    failed_trial = None
    if not feed_unstable:
        failed_trial = "feed_not_unstable"
    elif not final_phases_stable:
        failed_trial = "final_phase_unstable"
    return {
        "accepted": accepted,
        "status": "accepted" if accepted else "failed_gate",
        "tolerance": stability_tolerance,
        "feed_unstable": feed_unstable,
        "final_phases_stable": final_phases_stable,
        "failure_reason": "" if accepted else failed_trial,
        "trials": trial_records,
    }


def _electrolyte_lle_tpdf_seed_set(
    *,
    feed: np.ndarray,
    phases: tuple[EquilibriumPhase, EquilibriumPhase],
    charges: np.ndarray,
) -> tuple[tuple[str, np.ndarray], ...]:
    feed = np.asarray(feed, dtype=float)
    charged_total = float(np.sum(feed[np.abs(charges) > 1.0e-12]))
    seeds: list[tuple[str, np.ndarray]] = [
        ("feed_like", feed),
    ]
    for phase in phases:
        seeds.append((f"returned_{phase.label}_phase", np.asarray(phase.composition, dtype=float)))
    seeds.extend(
        [
            (
                "water_rich",
                _charge_neutral_electrolyte_seed(
                    feed=feed,
                    charges=charges,
                    neutral_weights=(0.985, 0.015),
                    ion_fraction=max(charged_total, 1.0e-6),
                ),
            ),
            (
                "organic_rich",
                _charge_neutral_electrolyte_seed(
                    feed=feed,
                    charges=charges,
                    neutral_weights=(0.05, 0.95),
                    ion_fraction=max(min(charged_total * 0.05, 1.0e-3), 1.0e-6),
                ),
            ),
            (
                "salt_rich",
                _charge_neutral_electrolyte_seed(
                    feed=feed,
                    charges=charges,
                    neutral_weights=None,
                    ion_fraction=min(max(charged_total * 5.0, 0.05), 0.2),
                ),
            ),
        ]
    )
    return tuple(seeds)


def _charge_neutral_electrolyte_seed(
    *,
    feed: np.ndarray,
    charges: np.ndarray,
    neutral_weights: tuple[float, float] | None,
    ion_fraction: float,
) -> np.ndarray:
    feed = np.asarray(feed, dtype=float)
    charges = np.asarray(charges, dtype=float)
    neutral_indices = np.flatnonzero(np.abs(charges) <= 1.0e-12)
    charged_indices = np.flatnonzero(np.abs(charges) > 1.0e-12)
    if neutral_indices.size < 2 or charged_indices.size == 0:
        return feed.copy()
    ion_fraction = float(np.clip(ion_fraction, 1.0e-8, 0.5))
    out = np.zeros_like(feed, dtype=float)
    charged_template = feed[charged_indices]
    charged_sum = float(np.sum(charged_template))
    if charged_sum <= 0.0:
        return feed.copy()
    out[charged_indices] = charged_template / charged_sum * ion_fraction
    if neutral_weights is None:
        neutral_template = feed[neutral_indices]
    else:
        neutral_template = np.ones(neutral_indices.size, dtype=float) * _ELECTROLYTE_LLE_TRACE_FLOOR
        neutral_template[0] = float(neutral_weights[0])
        neutral_template[1] = float(neutral_weights[1])
    neutral_template = np.asarray(neutral_template, dtype=float)
    neutral_template = np.maximum(neutral_template, _ELECTROLYTE_LLE_TRACE_FLOOR)
    neutral_template /= float(np.sum(neutral_template))
    out[neutral_indices] = neutral_template * (1.0 - ion_fraction)
    return out / float(np.sum(out))


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
    phase_models: Mapping[str, Any] | None = None,
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
    extra_phase_kwargs = dict(phase_kwargs or {})
    if phase_models is not None:
        extra_phase_kwargs["phase_models"] = phase_models
    route = _normalize_reactive_phase_route(mixture, phase_kind, extra_phase_kwargs)
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

    balance_matrix, total_vector, _ = _normalize_reactive_balances(species, balances, totals)
    reaction_defs = _normalize_reactions(species, reactions)
    reaction_phase_stoichiometry, _ = _reaction_phase_stoichiometry_matrix(species, reaction_defs, route)
    _require_reactive_phase_standard_states(reaction_defs, route)
    phase_model_payload = _normalize_reactive_phase_models(
        species,
        extra_phase_kwargs.get("phase_models"),
        route,
    )
    if route == "lle_flash":
        route_binding = "_native_reactive_lle_eos_route_result"
        required_route = "reactive_lle"
        problem_kind = "reactive_lle"
        phase_labels = ("liq1", "liq2")
    else:
        route_binding = "_native_reactive_electrolyte_lle_eos_route_result"
        required_route = "reactive_electrolyte_lle"
        problem_kind = "reactive_electrolyte_lle"
        phase_labels = ("aq", "org")
    return _native_reactive_two_phase_flash(
        mixture,
        T=T,
        P=P,
        feed=feed,
        balance_matrix=balance_matrix,
        total_vector=total_vector,
        species=species,
        reactions=reaction_defs,
        reaction_phase_stoichiometry=reaction_phase_stoichiometry,
        options=solver_options,
        route_binding=route_binding,
        required_route=required_route,
        problem_kind=problem_kind,
        phase_labels=phase_labels,
        phase_models=phase_model_payload,
    )


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
        if (
            kwargs.get("solvent_feed") is not None
            or kwargs.get("salt_molality") is not None
            or kwargs.get("phase_models") is not None
        ):
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
    allowed = {"solvent_feed", "salt_molality", "phase_models"}
    unsupported = sorted(key for key, value in phase_kwargs.items() if value is not None and key not in allowed)
    if unsupported:
        raise InputError(
            "reactive phase equilibrium does not support phase_kwargs key(s): {}.".format(", ".join(unsupported))
        )
    if route == "lle_flash" and (
        phase_kwargs.get("solvent_feed") is not None or phase_kwargs.get("salt_molality") is not None
    ):
        raise InputError("solvent_feed and salt_molality require reactive_electrolyte_lle.")
    if route == "lle_flash" and phase_kwargs.get("phase_models") is not None:
        raise InputError("phase_models requires reactive_electrolyte_lle.")


def _normalize_reactive_phase_models(
    species: list[str],
    phase_models: Any,
    route: str,
) -> dict[str, Any] | None:
    if phase_models is None:
        return None
    if route != "electrolyte_lle":
        raise InputError("phase_models requires reactive_electrolyte_lle.")
    if not isinstance(phase_models, Mapping):
        raise InputError("phase_models must be a mapping with 'aq' and 'org' mixtures.")
    try:
        aq_model = phase_models["aq"]
        org_model = phase_models["org"]
    except KeyError as exc:
        raise InputError("phase_models must provide 'aq' and 'org' mixtures.") from exc
    species_index = {label: index for index, label in enumerate(species)}

    def indices_for(model: Any, label: str) -> list[int]:
        model_species = [str(item) for item in getattr(model, "species", [])]
        if not model_species:
            raise InputError(f"phase_models['{label}'] must expose species labels.")
        missing = [item for item in model_species if item not in species_index]
        if missing:
            raise InputError(
                f"phase_models['{label}'] species are not present in the global reactive mixture: {missing}."
            )
        if len(set(model_species)) != len(model_species):
            raise InputError(f"phase_models['{label}'] species labels must be unique.")
        if not hasattr(model, "_native"):
            raise InputError(f"phase_models['{label}'] must be an ePCSAFTMixture.")
        return [species_index[item] for item in model_species]

    return {
        "aq": aq_model,
        "org": org_model,
        "aq_indices": indices_for(aq_model, "aq"),
        "org_indices": indices_for(org_model, "org"),
    }


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


def _require_reactive_phase_standard_states(reactions: list[Any], route: str) -> None:
    supported = {"mole_fraction_activity", "thermodynamic_activity"}
    unsupported = sorted(
        {str(reaction.standard_state) for reaction in reactions if str(reaction.standard_state) not in supported}
    )
    if unsupported:
        allowed = "', '".join(sorted(supported))
        observed = "', '".join(unsupported)
        raise InputError(f"{route} requires reaction standard_state in '{allowed}'; got '{observed}'.")


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
    composition = _normalize_feed(x, mixture.ncomp, opts.min_composition, "bubble_t")
    _reject_ion_containing_mixture(mixture)
    pressure = _positive_scalar(P, "P", "bubble_t")
    return _native_neutral_fixed_pressure_temperature(
        mixture,
        P=pressure,
        composition=composition,
        options=opts,
        route_label="bubble_t",
        route_binding="_native_neutral_bubble_t_eos_route_result",
        problem_kind="neutral_bubble_t",
    )


def dew_t(mixture: Any, *, P: float, y: Any, options: EquilibriumOptions | None = None) -> EquilibriumResult:
    """Solve a neutral dew temperature at fixed vapor composition and pressure."""
    opts = _normalize_options(options)
    composition = _normalize_feed(y, mixture.ncomp, opts.min_composition, "dew_t")
    _reject_ion_containing_mixture(mixture)
    pressure = _positive_scalar(P, "P", "dew_t")
    return _native_neutral_fixed_pressure_temperature(
        mixture,
        P=pressure,
        composition=composition,
        options=opts,
        route_label="dew_t",
        route_binding="_native_neutral_dew_t_eos_route_result",
        problem_kind="neutral_dew_t",
    )


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
    feed = _normalize_feed(z, mixture.ncomp, opts.min_composition, "stability")
    _reject_ion_containing_mixture(mixture)
    temperature = _positive_scalar(T, "T", "stability")
    pressure = _positive_scalar(P, "P", "stability")
    parents = _normalize_parent_phases(parent_phase)
    trials = _normalize_trial_phases(trial_phases)
    return _native_neutral_stability(
        mixture,
        T=temperature,
        P=pressure,
        feed=feed,
        options=opts,
        parent_phases=parents,
        trial_phases=trials,
    )


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
    temperature = _positive_scalar(T, "T", "electrolyte_stability")
    pressure = _positive_scalar(P, "P", "electrolyte_stability")
    feed, feed_diagnostics = _normalize_electrolyte_feed(
        mixture,
        z=z,
        solvent_feed=solvent_feed,
        salt_molality=salt_molality,
        options=opts,
    )
    charges = _mixture_charges(mixture)
    _require_charge_neutral(feed, charges, "electrolyte_stability feed")
    electrolyte_basis_diagnostics = electrolyte_formula_basis(
        mixture.species,
        charges,
        feed,
        salt_labels=tuple(feed_diagnostics.get("salt_molality", {})),
    )
    return _native_electrolyte_stability(
        mixture,
        T=temperature,
        P=pressure,
        feed=feed,
        feed_diagnostics=feed_diagnostics,
        options=opts,
    )


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
    electrolyte_basis_diagnostics = electrolyte_formula_basis(
        mixture.species,
        charges,
        feed,
        salt_labels=tuple(feed_diagnostics.get("salt_molality", {})),
    )
    from . import _core

    def run_route(route_options: EquilibriumOptions) -> tuple[Mapping[str, Any], tuple[float, float, float, float]]:
        route_tolerances = neutral_two_phase_eos_tolerances(P, route_options)
        material_tolerance, pressure_tolerance, chemical_potential_tolerance, phase_distance_tolerance = route_tolerances
        charge_tolerance = min(route_options.tolerance, 1.0e-8)
        route_result = _core._native_electrolyte_lle_eos_route_result(
            mixture._native,
            T,
            P,
            feed.tolist(),
            route_options.max_iterations,
            route_options.tolerance,
            _native_timeout_seconds(route_options),
            material_tolerance,
            pressure_tolerance,
            charge_tolerance,
            chemical_potential_tolerance,
            phase_distance_tolerance,
        )
        return route_result, route_tolerances

    route_options = opts
    route, route_tolerances = run_route(route_options)
    postsolve = route.get("postsolve", {})
    if (
        bool(route.get("accepted", False))
        and isinstance(postsolve, Mapping)
        and float(postsolve.get("ln_fugacity_consistency_norm", 0.0)) > 1.0e-7
        and opts.tolerance > 1.0e-8
    ):
        route_options = EquilibriumOptions(
            max_iterations=max(opts.max_iterations, 500),
            tolerance=1.0e-8,
            min_composition=opts.min_composition,
            jacobian_backend=opts.jacobian_backend,
            solver_backend=opts.solver_backend,
            timeout_seconds=opts.timeout_seconds,
        )
        route, route_tolerances = run_route(route_options)
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
        feed_diagnostics=feed_diagnostics,
        electrolyte_basis_diagnostics=electrolyte_basis_diagnostics,
        options=route_options,
    )
