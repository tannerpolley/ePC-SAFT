"""Native-backed phase-equilibrium helpers and Python input adapters."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from numbers import Integral, Real
from typing import Any, Literal

import numpy as np

from ._types import InputError, SolutionError
from .equilibrium_core.electrolyte_basis import build_electrolyte_basis

_ASCANI_2022_REFERENCE = {
    "authors": "Ascani, Sadowski, and Held",
    "year": 2022,
    "title": "Calculation of Multiphase Equilibria Containing Mixed Solvents and Mixed Electrolytes",
    "doi": "10.1021/acs.jced.1c00866",
}


def _raise_native_ipopt_not_routed(route: str) -> None:
    raise InputError(
        f"solver_backend='ipopt' was requested, but the native Ipopt adapter for {route} is not wired to public "
        "equilibrium routes yet. Use the current native route until the Ipopt NLP adapter is implemented."
    )


def _raise_native_ipopt_equilibrium_required(route: str) -> None:
    raise InputError(
        f"{route} requires a native Ipopt equilibrium NLP route. The previous Python scalar solve route was removed "
        "by the solver gate."
    )


@dataclass(frozen=True, slots=True)
class EquilibriumOptions:
    """Numerical controls for equilibrium solvers."""

    max_iterations: int = 180
    tolerance: float = 1.0e-6
    min_composition: float = 1.0e-12
    include_phase_diagnostics: bool = False
    stability_precheck: bool = True
    ignored_legacy_options: tuple[str, ...] = ()
    density_diagnostics: Literal["auto", "off", "full"] = "auto"
    experimental_coupled_density_lle: bool = False
    jacobian_backend: Literal["auto", "analytic", "cppad"] = "auto"
    solver_backend: Literal["auto", "ipopt"] = "auto"
    timeout_seconds: float | None = None
    max_seed_attempts: int | None = None
    max_density_failures: int | None = None
    max_total_objective_evaluations: int | None = None


@dataclass(frozen=True, slots=True)
class EquilibriumProblem:
    """Base class for typed equilibrium problem objects."""

    def solve(self, mixture):
        """Solve this problem with a mixture instance."""
        raise NotImplementedError("EquilibriumProblem subclasses define solve(mixture).")


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
    initial_phases: Any | None = None

    def solve(self, mixture):
        return mixture.lle_tp(
            T=self.T,
            P=self.P,
            z=self.z,
            options=self.options,
            initial_phases=self.initial_phases,
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
    initial_phases: Any | None = None

    def solve(self, mixture):
        return mixture.electrolyte_lle_tp(
            T=self.T,
            P=self.P,
            z=self.z,
            solvent_feed=self.solvent_feed,
            salt_molality=self.salt_molality,
            initial_phases=self.initial_phases,
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
    """Homogeneous activity-coupled reactive speciation problem."""

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

    P_seed: float | None = None
    vapor_species: Any | None = None
    volatile_species: Any | None = None
    nonvolatile_species: Any | None = None

    def solve(self, mixture):
        return mixture.reactive_electrolyte_bubble_p(
            T=self.T,
            P_seed=self.P if self.P_seed is None else self.P_seed,
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
        coefficient field. ``fugacity_coefficient`` remains accepted as a
        compatibility alias for older callers and carries the same ln(phi)
        values, not coefficient-form phi values.
        """
        if ln_fugacity_coefficient is None:
            ln_fugacity_coefficient = fugacity_coefficient
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
        """Backward-compatible alias for natural-log fugacity coefficients."""
        return self.ln_fugacity_coefficient

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
            "fugacity_coefficient": ln_fugacity,
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
    initial_phases: Any = None,
    options: EquilibriumOptions | None = None,
) -> EquilibriumResult:
    """Run native C++ electrolyte LLE through the public compatibility name."""
    return electrolyte_lle_flash_native(
        mixture,
        T=T,
        P=P,
        z=z,
        solvent_feed=solvent_feed,
        salt_molality=salt_molality,
        initial_phases=initial_phases,
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
    basis_mass_kg = _positive_scalar(basis_mass_kg, "basis_mass_kg", "electrolyte_feed_from_molality")
    solvent_x = _normalize_solvent_feed(species, charges, solvent_feed)
    salt_items = _normalize_salt_molality(salt_molality)
    neutral_mw = float(np.sum(solvent_x * mw))
    if not np.isfinite(neutral_mw) or neutral_mw <= 0.0:
        raise InputError("solvent_feed produced an invalid salt-free solvent molecular weight.")
    moles = np.zeros(len(species), dtype=float)
    neutral_moles = basis_mass_kg / neutral_mw
    moles += solvent_x * neutral_moles
    for salt_label, molality in salt_items.items():
        cation_i, anion_i = _species_pair_for_salt(species, charges, salt_label)
        z_cat = round(abs(float(charges[cation_i])))
        z_an = round(abs(float(charges[anion_i])))
        gcd_z = int(np.gcd(z_cat, z_an))
        nu_cation = z_an // gcd_z
        nu_anion = z_cat // gcd_z
        salt_moles = float(molality) * basis_mass_kg
        moles[cation_i] += salt_moles * nu_cation
        moles[anion_i] += salt_moles * nu_anion
    total = float(np.sum(moles))
    if total <= 0.0:
        raise InputError("Computed electrolyte feed has non-positive total moles.")
    feed = moles / total
    _require_charge_neutral(feed, charges, "molality-derived electrolyte feed")
    return feed


def _normalize_options(options: EquilibriumOptions | Mapping[str, Any] | None) -> EquilibriumOptions:
    if options is None:
        return EquilibriumOptions()
    if isinstance(options, Mapping):
        raw = dict(options)
        ignored_keys = []
        translated: dict[str, Any] = {}
        legacy_map = {
            "max_nfev": "max_iterations",
            "solver_tol": "tolerance",
        }
        for source, target in legacy_map.items():
            if source in raw:
                translated[target] = raw.pop(source)
        ignored_legacy = {"tpdf_global_trials", "tpdf_local_trials", "charge_weight", "seed_x", "force_seed_solve"}
        for key in sorted(ignored_legacy):
            if key in raw:
                raw.pop(key)
                ignored_keys.append(key)
        allowed = {
            "max_iterations",
            "tolerance",
            "min_composition",
            "include_phase_diagnostics",
            "stability_precheck",
            "density_diagnostics",
            "experimental_coupled_density_lle",
            "jacobian_backend",
            "solver_backend",
            "timeout_seconds",
            "max_seed_attempts",
            "max_density_failures",
            "max_total_objective_evaluations",
        }
        unknown = sorted(set(raw) - allowed)
        if unknown:
            raise InputError("Unknown equilibrium option key(s): {}.".format(", ".join(unknown)))
        translated.update(raw)
        translated["ignored_legacy_options"] = tuple(ignored_keys)
        options = EquilibriumOptions(**translated)
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
    if not isinstance(options.include_phase_diagnostics, bool):
        raise InputError("options.include_phase_diagnostics must be boolean.")
    if not isinstance(options.stability_precheck, bool):
        raise InputError("options.stability_precheck must be boolean.")
    ignored_legacy_options = tuple(str(item) for item in options.ignored_legacy_options)
    density_diagnostics = str(options.density_diagnostics).strip().lower()
    if density_diagnostics not in {"auto", "off", "full"}:
        raise InputError("options.density_diagnostics must be 'auto', 'off', or 'full'.")
    if not isinstance(options.experimental_coupled_density_lle, bool):
        raise InputError("options.experimental_coupled_density_lle must be boolean.")
    jacobian_backend = str(options.jacobian_backend).strip().lower()
    if jacobian_backend not in {"auto", "analytic", "cppad"}:
        raise InputError("options.jacobian_backend must be 'auto', 'analytic', or 'cppad'.")
    solver_backend = str(options.solver_backend).strip().lower()
    if solver_backend not in {"auto", "ipopt"}:
        raise InputError("options.solver_backend must be 'auto' or 'ipopt'.")
    timeout_seconds = _optional_positive_float_option(options.timeout_seconds, "timeout_seconds")
    max_seed_attempts = _optional_positive_int_option(options.max_seed_attempts, "max_seed_attempts")
    max_density_failures = _optional_positive_int_option(options.max_density_failures, "max_density_failures")
    max_total_objective_evaluations = _optional_positive_int_option(
        options.max_total_objective_evaluations, "max_total_objective_evaluations"
    )
    return EquilibriumOptions(
        max_iterations=max_iterations,
        tolerance=tolerance,
        min_composition=min_composition,
        include_phase_diagnostics=options.include_phase_diagnostics,
        stability_precheck=options.stability_precheck,
        ignored_legacy_options=ignored_legacy_options,
        density_diagnostics=density_diagnostics,  # type: ignore[arg-type]
        experimental_coupled_density_lle=options.experimental_coupled_density_lle,
        jacobian_backend=jacobian_backend,  # type: ignore[arg-type]
        solver_backend=solver_backend,  # type: ignore[arg-type]
        timeout_seconds=timeout_seconds,
        max_seed_attempts=max_seed_attempts,
        max_density_failures=max_density_failures,
        max_total_objective_evaluations=max_total_objective_evaluations,
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


def _optional_positive_int_option(value: Any, label: str) -> int | None:
    if value is None:
        return None
    if isinstance(value, bool) or not isinstance(value, Integral):
        raise InputError(f"options.{label} must be an integer greater than zero when provided.")
    out = int(value)
    if out <= 0:
        raise InputError(f"options.{label} must be an integer greater than zero when provided.")
    return out


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
        "salt_molality": dict(_normalize_salt_molality(salt_molality)),
        "solvent_feed": _json_like(solvent_feed),
    }


def _normalize_solvent_feed(species: list[str], charges: np.ndarray, solvent_feed: Any) -> np.ndarray:
    neutral_indices = [i for i, charge in enumerate(charges) if abs(float(charge)) <= 1.0e-12]
    if not neutral_indices:
        raise InputError("solvent_feed requires at least one neutral solvent species.")
    solvent_x = np.zeros(len(species), dtype=float)
    if isinstance(solvent_feed, dict):
        for label, value in solvent_feed.items():
            try:
                index = species.index(str(label))
            except ValueError as exc:
                raise InputError(f"solvent_feed species '{label}' is not present in the mixture.") from exc
            if index not in neutral_indices:
                raise InputError(f"solvent_feed species '{label}' is ionic; expected neutral solvents only.")
            solvent_x[index] = float(value)
    else:
        values = np.asarray(solvent_feed, dtype=float).flatten()
        if values.size != len(neutral_indices):
            raise InputError("solvent_feed must be a dict or a vector with one entry per neutral solvent species.")
        for index, value in zip(neutral_indices, values):
            solvent_x[index] = float(value)
    if not np.all(np.isfinite(solvent_x)) or np.any(solvent_x < 0.0):
        raise InputError("solvent_feed must contain finite non-negative values.")
    total = float(np.sum(solvent_x))
    if total <= 0.0:
        raise InputError("solvent_feed must have a positive sum.")
    return solvent_x / total


def _normalize_salt_molality(salt_molality: Any) -> dict[str, float]:
    if not isinstance(salt_molality, dict) or not salt_molality:
        raise InputError("salt_molality must be a non-empty dict like {'NaCl': 1.0}.")
    out: dict[str, float] = {}
    for label, value in salt_molality.items():
        molality = float(value)
        if not np.isfinite(molality) or molality < 0.0:
            raise InputError("salt_molality values must be finite and non-negative.")
        out[str(label)] = molality
    return out


def _species_pair_for_salt(species: list[str], charges: np.ndarray, salt_label: str) -> tuple[int, int]:
    normalized = _salt_label_token(salt_label)
    cation_indices = [i for i, charge in enumerate(charges) if float(charge) > 1.0e-12]
    anion_indices = [i for i, charge in enumerate(charges) if float(charge) < -1.0e-12]
    matches: list[tuple[int, int]] = []
    for cation_i in cation_indices:
        for anion_i in anion_indices:
            cation_stoich, anion_stoich = _salt_stoichiometry(charges[cation_i], charges[anion_i])
            cation_label = _ion_stem(species[cation_i], charges[cation_i])
            anion_label = _ion_stem(species[anion_i], charges[anion_i])
            pair_label = _salt_label_token(
                cation_label
                + ("" if cation_stoich == 1 else str(cation_stoich))
                + anion_label
                + ("" if anion_stoich == 1 else str(anion_stoich))
            )
            if pair_label == normalized:
                matches.append((cation_i, anion_i))
    if len(matches) != 1:
        raise InputError(f"Could not uniquely map salt_molality key '{salt_label}' onto mixture ions.")
    return matches[0]


def _salt_label_token(label: Any) -> str:
    return "".join(ch for ch in str(label) if ch.isalnum()).lower()


def _ion_stem(label: str, charge: float | None = None) -> str:
    text = str(label)
    stripped = text.replace("+", "").replace("-", "")
    if charge is not None and ("+" in text or "-" in text):
        charge_int = round(abs(float(charge)))
        suffix = str(charge_int)
        if charge_int > 1 and stripped.endswith(suffix):
            stripped = stripped[: -len(suffix)]
    return stripped


def _salt_stoichiometry(cation_charge: float, anion_charge: float) -> tuple[int, int]:
    cation_int = round(abs(float(cation_charge)))
    anion_int = round(abs(float(anion_charge)))
    if cation_int <= 0 or anion_int <= 0:
        raise InputError("electrolyte salt stoichiometry requires non-zero ion charges.")
    if (
        abs(float(cation_int) - abs(float(cation_charge))) > 1.0e-12
        or abs(float(anion_int) - abs(float(anion_charge))) > 1.0e-12
    ):
        raise InputError("electrolyte salt stoichiometry currently requires integer ion charges.")
    gcd = int(np.gcd(cation_int, anion_int))
    return anion_int // gcd, cation_int // gcd


def _electrolyte_formula_basis(mixture: Any, feed: np.ndarray, feed_diagnostics: dict[str, Any]) -> dict[str, Any]:
    charges = _mixture_charges(mixture)
    species = list(mixture.species)
    neutral_indices = [i for i, charge in enumerate(charges) if abs(float(charge)) <= 1.0e-12]
    cation_indices = [i for i, charge in enumerate(charges) if float(charge) > 1.0e-12]
    anion_indices = [i for i, charge in enumerate(charges) if float(charge) < -1.0e-12]
    if len(neutral_indices) < 2:
        raise InputError("electrolyte_lle requires at least two neutral solvent species.")
    if not cation_indices or not anion_indices:
        raise InputError("electrolyte_lle requires at least one cation and one anion.")
    salt_labels = tuple(feed_diagnostics.get("salt_molality", {}).keys())
    basis = build_electrolyte_basis(species, charges, feed, salt_labels=tuple(salt_labels))
    payload = basis.to_dict()
    return {
        "neutral_indices": tuple(payload["neutral_indices"]),
        "charged_indices": tuple(payload["charged_indices"]),
        "cation_indices": tuple(payload["cation_indices"]),
        "anion_indices": tuple(payload["anion_indices"]),
        "salt_pairs": tuple(payload["salt_pairs"]),
        "formula_feed": np.asarray(payload["formula_feed"], dtype=float),
        "e_matrix": np.asarray(payload["e_matrix"], dtype=float),
        "basis_rank": int(payload["basis_rank"]),
        "variable_model": str(payload["variable_model"]),
    }


def _formula_to_explicit_composition(
    formula_composition: np.ndarray, basis: dict[str, Any], ncomp: int
) -> tuple[np.ndarray, float]:
    formula = np.asarray(formula_composition, dtype=float)
    explicit = np.zeros(int(ncomp), dtype=float)
    neutral_indices = basis["neutral_indices"]
    salt_pairs = basis["salt_pairs"]
    for pos, index in enumerate(neutral_indices):
        explicit[int(index)] += float(formula[pos])
    offset = len(neutral_indices)
    for salt_pos, pair in enumerate(salt_pairs):
        amount = float(formula[offset + salt_pos])
        explicit[int(pair["cation"])] += amount * float(pair.get("cation_stoich", 1.0))
        explicit[int(pair["anion"])] += amount * float(pair.get("anion_stoich", 1.0))
    total = float(np.sum(explicit))
    if total <= 0.0:
        raise SolutionError("Formula-basis electrolyte phase expanded to a non-positive explicit composition.")
    return explicit / total, total


def _explicit_to_formula_composition(composition: np.ndarray, basis: dict[str, Any]) -> np.ndarray:
    comp = np.asarray(composition, dtype=float)
    values = [float(comp[index]) for index in basis["neutral_indices"]]
    values.extend(
        float(comp[int(pair["cation"])]) / float(pair.get("cation_stoich", 1.0)) for pair in basis["salt_pairs"]
    )
    out = np.asarray(values, dtype=float)
    total = float(np.sum(out))
    if total <= 0.0:
        raise InputError("Explicit electrolyte composition cannot be represented on the formula basis.")
    return out / total


def _electrolyte_initial_phase_seed(
    mixture: Any,
    feed: np.ndarray,
    basis: dict[str, Any],
    initial_phases: Any,
    options: EquilibriumOptions,
) -> dict[str, Any]:
    if not isinstance(initial_phases, dict):
        raise InputError("initial_phases for electrolyte_lle must be a dict with 'aq', 'org', and 'phase_fraction'.")
    required = {"aq", "org", "phase_fraction"}
    keys = set(initial_phases)
    if keys != required:
        raise InputError("initial_phases for electrolyte_lle must contain exactly 'aq', 'org', and 'phase_fraction'.")
    aq_comp = _normalize_initial_phase(initial_phases["aq"], mixture.ncomp, options.min_composition, "aq")
    org_comp = _normalize_initial_phase(initial_phases["org"], mixture.ncomp, options.min_composition, "org")
    beta_org = float(initial_phases["phase_fraction"])
    if not np.isfinite(beta_org) or beta_org <= 0.0 or beta_org >= 1.0:
        raise InputError("initial_phases phase_fraction must be > 0 and < 1.")
    charges = _mixture_charges(mixture)
    phase_charge_error = max(abs(float(np.dot(aq_comp, charges))), abs(float(np.dot(org_comp, charges))))
    if phase_charge_error > 1.0e-8:
        raise InputError("initial_phases aq and org must be charge neutral for electrolyte_lle.")
    material_error = float(np.max(np.abs((1.0 - beta_org) * aq_comp + beta_org * org_comp - feed)))
    if material_error > 1.0e-7:
        raise InputError("initial_phases aq/org/phase_fraction must reconstruct the electrolyte_lle feed.")
    aq_formula = _explicit_to_formula_composition(aq_comp, basis)
    org_formula = _explicit_to_formula_composition(org_comp, basis)
    beta_formula = _explicit_beta_to_formula_beta(beta_org, aq_formula, org_formula, basis, mixture.ncomp)
    return {
        "seed_name": "initial_phases",
        "beta_formula": beta_formula,
        "aq_formula": aq_formula,
        "org_formula": org_formula,
        "fixture": None,
        "diagnostics": {
            "initial_phase_material_balance_error": material_error,
            "initial_phase_charge_balance_error": phase_charge_error,
        },
    }


def _explicit_beta_to_formula_beta(
    beta_explicit: float, aq_formula: np.ndarray, org_formula: np.ndarray, basis: dict[str, Any], ncomp: int
) -> float:
    _aq_exp, aq_scale = _formula_to_explicit_composition(aq_formula, basis, ncomp)
    _org_exp, org_scale = _formula_to_explicit_composition(org_formula, basis, ncomp)
    numerator = float(beta_explicit) / org_scale
    denominator = numerator + (1.0 - float(beta_explicit)) / aq_scale
    if denominator <= 0.0:
        return float(beta_explicit)
    return float(np.clip(numerator / denominator, 1.0e-12, 1.0 - 1.0e-12))


def _phase_state(
    mixture: Any,
    T: float,
    P: float,
    composition: np.ndarray,
    label: str,
    options: EquilibriumOptions,
    context: str,
) -> dict[str, Any]:
    try:
        state = mixture.state(T=T, P=P, x=composition, phase=label)
    except SolutionError:
        raise
    except (InputError, ValueError, RuntimeError, ArithmeticError) as exc:
        raise SolutionError(f"Failed to construct {label} phase during {context}: {exc}") from exc
    diagnostics = state.state_diagnostics(species=mixture.species) if options.include_phase_diagnostics else None
    return {
        "state": state,
        "ln_phi": np.asarray(state.fugacity_coefficient(), dtype=float),
        "density": float(state.density()),
        "diagnostics": diagnostics,
    }


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


def _normalize_initial_phase(value: Any, ncomp: int, min_composition: float, label: str) -> np.ndarray:
    composition = np.asarray(value, dtype=float).flatten()
    if composition.size != int(ncomp):
        raise InputError(
            f"initial_phases {label} length ({composition.size}) must match mixture component count ({ncomp})."
        )
    if not np.all(np.isfinite(composition)):
        raise InputError(f"initial_phases {label} must contain only finite values.")
    if np.any(composition < 0.0):
        raise InputError(f"initial_phases {label} must be non-negative.")
    total = float(np.sum(composition))
    if total <= 0.0:
        raise InputError(f"initial_phases {label} must have a positive sum.")
    composition = composition / total
    if np.any(composition < min_composition):
        raise InputError(f"initial_phases {label} entries must be >= min_composition.")
    return composition


def _normalize_neutral_initial_phases(initial_phases: Any, ncomp: int, min_composition: float) -> dict[str, Any]:
    if not isinstance(initial_phases, dict):
        raise InputError("initial_phases must be a dict with 'liq1', 'liq2', and 'phase_fraction'.")
    missing = {"liq1", "liq2", "phase_fraction"} - set(initial_phases)
    if missing:
        raise InputError("initial_phases is missing required key(s): {}.".format(", ".join(sorted(missing))))
    comp1 = _normalize_initial_phase(initial_phases["liq1"], ncomp, min_composition, "liq1")
    comp2 = _normalize_initial_phase(initial_phases["liq2"], ncomp, min_composition, "liq2")
    beta = float(initial_phases["phase_fraction"])
    if not np.isfinite(beta) or not (0.0 < beta < 1.0):
        raise InputError("initial_phases phase_fraction must be > 0 and < 1.")
    return {"liq1": comp1.tolist(), "liq2": comp2.tolist(), "phase_fraction": beta}


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


def _options_to_native_dict(options: EquilibriumOptions) -> dict[str, Any]:
    return {
        "max_iterations": int(options.max_iterations),
        "tolerance": float(options.tolerance),
        "min_composition": float(options.min_composition),
        "include_phase_diagnostics": bool(options.include_phase_diagnostics),
        "stability_precheck": bool(options.stability_precheck),
        "density_diagnostics": str(options.density_diagnostics),
        "experimental_coupled_density_lle": bool(options.experimental_coupled_density_lle),
        "jacobian_backend": str(options.jacobian_backend),
        "solver_backend": str(options.solver_backend),
        "timeout_seconds": None if options.timeout_seconds is None else float(options.timeout_seconds),
        "max_seed_attempts": None if options.max_seed_attempts is None else int(options.max_seed_attempts),
        "max_density_failures": None if options.max_density_failures is None else int(options.max_density_failures),
        "max_total_objective_evaluations": (
            None if options.max_total_objective_evaluations is None else int(options.max_total_objective_evaluations)
        ),
    }


def _native_phase_from_payload(payload: dict[str, Any]) -> EquilibriumPhase:
    return EquilibriumPhase(
        label=payload["label"],
        composition=payload["composition"],
        density=payload["density"],
        temperature=payload["temperature"],
        pressure=payload["pressure"],
        phase_fraction=payload["phase_fraction"],
        ln_fugacity_coefficient=payload.get("ln_fugacity_coefficient"),
        diagnostics=payload.get("diagnostics") or None,
    )


def _native_trial_from_payload(payload: dict[str, Any]) -> StabilityTrial:
    return StabilityTrial(
        parent_phase=payload["parent_phase"],
        trial_phase=payload["trial_phase"],
        seed_name=payload["seed_name"],
        composition=payload["composition"],
        tpd=payload["tpd"],
        iterations=payload["iterations"],
        converged=payload["converged"],
        unstable=payload["unstable"],
        diagnostics=payload.get("diagnostics") or {},
    )


def _residual_norm_from_diagnostics(diagnostics: dict[str, Any]) -> float | None:
    for key in ("solver_residual_norm", "fugacity_residual_norm", "residual_norm", "best_fugacity_residual_norm"):
        if key in diagnostics:
            try:
                return float(diagnostics[key])
            except (TypeError, ValueError):
                return None
    return None


def _phase_charge_balance_from_phases(feed: np.ndarray, phases: Any, charges: np.ndarray) -> dict[str, float]:
    charge_balance = {"feed": float(np.dot(feed, charges))}
    max_abs = abs(charge_balance["feed"])
    for index, phase in enumerate(phases):
        label = str(getattr(phase, "label", "") or f"phase_{index}")
        value = float(np.dot(np.asarray(phase.composition, dtype=float), charges))
        charge_balance[label] = value
        max_abs = max(max_abs, abs(value))
    charge_balance["max_abs"] = max_abs
    return charge_balance


def _solved_internal_states(problem_kind: str) -> list[str]:
    by_kind = {
        "tp_flash": ["density_roots", "phase_compositions"],
        "bubble_p": ["density_roots", "bubble_pressure_root", "phase_compositions"],
        "bubble_t": ["density_roots", "bubble_temperature_root", "phase_compositions"],
        "dew_p": ["density_roots", "dew_pressure_root", "phase_compositions"],
        "dew_t": ["density_roots", "dew_temperature_root", "phase_compositions"],
        "lle_flash": ["density_roots", "phase_compositions"],
        "electrolyte_lle": ["density_roots", "charge_constrained_formula_basis", "phase_compositions"],
        "electrolyte_lle_flash": ["density_roots", "charge_constrained_formula_basis", "phase_compositions"],
        "stability": ["density_roots", "tpd_trial_compositions"],
        "electrolyte_stability": ["density_roots", "charge_constrained_formula_basis", "tpd_trial_compositions"],
    }
    return list(by_kind.get(problem_kind, ["density_roots"]))


def _derivative_backend_blocks(problem_kind: str, derivative_backend: str) -> dict[str, str]:
    if derivative_backend in {"analytic_implicit", "cppad_implicit"}:
        density_backend = "implicit_density_root"
    else:
        density_backend = "not_applicable"
    blocks: dict[str, str] = {
        "density_root": density_backend,
        "eos_state_properties": "analytic",
    }
    if "bubble" in problem_kind:
        blocks["bubble_pressure_root" if problem_kind == "bubble_p" else "bubble_or_dew_root"] = derivative_backend
    elif "dew" in problem_kind:
        blocks["dew_root"] = derivative_backend
    elif "lle" in problem_kind:
        blocks["lle_residual"] = derivative_backend
    elif "stability" in problem_kind:
        blocks["tpd_trial_residual"] = derivative_backend
    else:
        blocks["phase_equilibrium_residual"] = derivative_backend
    return blocks


def _route_diagnostics_for_problem_kind(problem_kind: str) -> dict[str, str]:
    if problem_kind in {"tp_flash", "bubble_p", "bubble_t", "dew_p", "dew_t"}:
        return {"equilibrium_route": "neutral_vle", "route_reason": "requested vapor-liquid path"}
    if problem_kind == "lle_flash":
        return {"equilibrium_route": "neutral_lle", "route_reason": "requested neutral liquid-liquid path"}
    if problem_kind == "stability":
        return {"equilibrium_route": "neutral_tpd", "route_reason": "requested neutral stability path"}
    if problem_kind in {"electrolyte_lle", "electrolyte_lle_flash"}:
        return {"equilibrium_route": "electrolyte_lle", "route_reason": "requested electrolyte liquid-liquid path"}
    if problem_kind == "electrolyte_stability":
        return {"equilibrium_route": "electrolyte_lle", "route_reason": "requested electrolyte stability path"}
    if problem_kind == "electrolyte_bubble_pressure":
        return {"equilibrium_route": "electrolyte_bubble", "route_reason": "requested electrolyte bubble-pressure path"}
    return {"equilibrium_route": "unsupported", "route_reason": "unsupported equilibrium kind"}


def _normalize_derivative_diagnostics(
    diagnostics: dict[str, Any],
    *,
    problem_kind: str,
    phase_count: int = 0,
) -> dict[str, Any]:
    diagnostics = dict(diagnostics)
    route_diagnostics = _route_diagnostics_for_problem_kind(problem_kind)
    diagnostics.setdefault("equilibrium_route", route_diagnostics["equilibrium_route"])
    diagnostics.setdefault("route_reason", route_diagnostics["route_reason"])
    derivative_backend = str(diagnostics.get("derivative_backend", "not_applicable"))
    diagnostics.setdefault("thermodynamic_backend", "epcsaft_state_fugacity_activity_property_api")
    diagnostics.setdefault(
        "solver_backend", diagnostics.get("nonlinear_solver", diagnostics.get("selected_solver_backend", "native"))
    )
    diagnostics.setdefault("requested_jacobian_backend", "auto")
    diagnostics.setdefault("derivative_backend", derivative_backend)
    diagnostics.setdefault("derivative_status", derivative_backend)
    if derivative_backend == "not_applicable":
        diagnostics.setdefault("derivative_available", False)
        diagnostics.setdefault("jacobian_available", False)
    residual_norm = _residual_norm_from_diagnostics(diagnostics)
    if residual_norm is not None:
        diagnostics.setdefault("residual_norm", residual_norm)
    diagnostics.setdefault("solved_internal_states", _solved_internal_states(problem_kind))
    diagnostics.setdefault("derivative_backend_by_block", _derivative_backend_blocks(problem_kind, derivative_backend))
    diagnostics.setdefault("implicit_sensitivity_blocks", [])
    residual_by_block: dict[str, float] = {}
    if "solver_residual_norm" in diagnostics:
        residual_by_block["solver"] = float(diagnostics["solver_residual_norm"])
    if "fugacity_residual_norm" in diagnostics:
        residual_by_block["fugacity"] = float(diagnostics["fugacity_residual_norm"])
    if "material_balance_error" in diagnostics:
        residual_by_block["material_balance"] = float(diagnostics["material_balance_error"])
    if "charge_balance_error" in diagnostics:
        residual_by_block["charge_balance"] = float(diagnostics["charge_balance_error"])
    diagnostics.setdefault("residual_norm_by_block", residual_by_block)
    best_available = bool(
        phase_count
        or diagnostics.get("best_noncollapsed_candidate") == "accepted"
        or diagnostics.get("best_noncollapsed_candidate") == "available"
        or "best_P" in diagnostics
    )
    diagnostics.setdefault("best_state_available", best_available)
    if best_available:
        diagnostics.setdefault("best_state", {"phase_count": int(phase_count), "source": "native_equilibrium_result"})
    diagnostics.setdefault(
        "row_failure_count", int(diagnostics.get("state_failure_count", diagnostics.get("density_failure_count", 0)))
    )
    diagnostics.setdefault("association_solver_status", "not_coupled")
    return diagnostics


def _native_result_from_payload(payload: dict[str, Any]) -> EquilibriumResult | StabilityResult:
    diagnostics = _normalize_derivative_diagnostics(
        payload.get("diagnostics") or {},
        problem_kind=str(payload.get("problem_kind", "")),
        phase_count=len(payload.get("phases", ()) or ()),
    )
    payload = dict(payload)
    payload["diagnostics"] = diagnostics
    if payload.get("result_type") == "stability":
        return StabilityResult(
            backend=payload["backend"],
            problem_kind=payload["problem_kind"],
            stable=payload["stable"],
            min_tpd=payload["min_tpd"],
            parent_phase=payload["parent_phase"],
            trial_phase=payload["trial_phase"],
            trial_composition=payload["trial_composition"],
            trials=tuple(_native_trial_from_payload(item) for item in payload.get("trials", ())),
            diagnostics=payload.get("diagnostics") or {},
        )
    return EquilibriumResult(
        backend=payload["backend"],
        problem_kind=payload["problem_kind"],
        phases=tuple(_native_phase_from_payload(item) for item in payload.get("phases", ())),
        stable=payload["stable"],
        split_detected=payload["split_detected"],
        diagnostics=payload.get("diagnostics") or {},
    )


def _call_native_equilibrium(
    mixture: Any,
    *,
    kind: str,
    T: float,
    P: float,
    z: Any,
    options: EquilibriumOptions,
    initial_phases: Any = None,
    parent_phase: Any = None,
    trial_phases: Any = None,
) -> EquilibriumResult | StabilityResult:
    if options.solver_backend == "ipopt":
        _raise_native_ipopt_not_routed(kind)
    from . import _core

    request: dict[str, Any] = {
        "kind": kind,
        "T": float(T),
        "P": float(P),
        "z": np.asarray(z, dtype=float).flatten().tolist(),
        "species": list(getattr(mixture, "species", [])),
        "options": _options_to_native_dict(options),
        "initial_phases": initial_phases,
    }
    if parent_phase is not None:
        request["parent_phases"] = list(_normalize_parent_phases(parent_phase))
    if trial_phases is not None:
        request["trial_phases"] = list(_normalize_trial_phases(trial_phases))
    try:
        payload = _core._solve_equilibrium_native(mixture._native, request)
    except _core.NativeValueError as exc:
        raise InputError(str(exc)) from exc
    except _core.NativeSolutionError as exc:
        message = str(exc.args[0]) if getattr(exc, "args", ()) else str(exc)
        diagnostics = exc.args[1] if len(getattr(exc, "args", ())) > 1 and isinstance(exc.args[1], dict) else None
        diagnostics = _diagnostics_with_options(options=options, diagnostics=diagnostics)
        raise SolutionError(message, diagnostics) from exc
    if kind in {"electrolyte_lle", "electrolyte_lle_flash"}:
        diagnostics = dict(payload.get("diagnostics") or {})
        _add_legacy_option_diagnostics(diagnostics, options)
        payload["diagnostics"] = diagnostics
    if payload.get("result_type") == "stability":
        diagnostics = dict(payload.get("diagnostics") or {})
        diagnostics.setdefault("parent_phases", request.get("parent_phases", ["liq", "vap"]))
        diagnostics.setdefault("trial_phases", request.get("trial_phases", ["liq", "vap"]))
        payload["diagnostics"] = diagnostics
    return _native_result_from_payload(payload)


def _add_legacy_option_diagnostics(diagnostics: dict[str, Any], options: EquilibriumOptions) -> None:
    diagnostics.setdefault("ignored_legacy_options", list(options.ignored_legacy_options))
    diagnostics.setdefault("density_diagnostics_mode", str(options.density_diagnostics))
    diagnostics.setdefault("experimental_coupled_density_lle", bool(options.experimental_coupled_density_lle))
    diagnostics.setdefault("requested_solver_backend", str(options.solver_backend))
    diagnostics.setdefault("selected_solver_backend", "native")
    diagnostics.setdefault(
        "solver_selection_reason", "default_native" if options.solver_backend == "auto" else "explicit_request"
    )
    diagnostics.setdefault("default_auto_uses_ipopt", False)
    diagnostics.setdefault("density_failure_count", 0)
    diagnostics.setdefault("density_failure_contexts", [])
    diagnostics.setdefault("density_scan_summary", {})
    diagnostics.setdefault("density_candidate_roots", [])
    diagnostics.setdefault("density_best_near_root", {})
    diagnostics.setdefault("density_best_candidate_refinement_used", False)
    diagnostics.setdefault("density_best_candidate_rejection_reason", "")
    diagnostics.setdefault("density_warm_start_source", "")
    diagnostics.setdefault("density_validity_gate", "not_evaluated")


def _diagnostics_with_options(
    *,
    options: EquilibriumOptions,
    diagnostics: dict[str, Any] | None,
) -> dict[str, Any] | None:
    if diagnostics is None:
        return None
    out = dict(diagnostics)
    _add_legacy_option_diagnostics(out, options)
    return out


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
    """Run the native coupled two-liquid reactive phase-equilibrium solve."""
    from . import _core
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
    reactive_options, solver_options = _reactive_phase_option_pair(
        options=options,
        phase_options=phase_options,
        normalize_reactive_options=_normalize_reactive_options,
    )
    temperature = _positive_scalar(T, "T", "reactive_phase_equilibrium")
    pressure = _positive_scalar(P, "P", "reactive_phase_equilibrium")
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
        basis_payload = _electrolyte_formula_basis(mixture, feed, feed_diagnostics)
        native_initial_phases = None
        if extra_phase_kwargs.get("initial_phases") is not None:
            initial_phases = extra_phase_kwargs["initial_phases"]
            _electrolyte_initial_phase_seed(mixture, feed, basis_payload, initial_phases, solver_options)
            native_initial_phases = {
                "aq": np.asarray(initial_phases["aq"], dtype=float).flatten().tolist(),
                "org": np.asarray(initial_phases["org"], dtype=float).flatten().tolist(),
                "phase_fraction": float(initial_phases["phase_fraction"]),
            }
    else:
        if extra_phase_kwargs.get("solvent_feed") is not None or extra_phase_kwargs.get("salt_molality") is not None:
            raise InputError("solvent_feed and salt_molality require reactive_electrolyte_lle.")
        _reject_ion_containing_mixture(mixture)
        charges = np.zeros(int(mixture.ncomp), dtype=float)
        feed_source = z if z is not None else initial_x
        feed = _normalize_feed(feed_source, int(mixture.ncomp), solver_options.min_composition, "reactive_lle")
        feed_diagnostics = {"composition_basis": "mole_fraction"}
        native_initial_phases = None
        if extra_phase_kwargs.get("initial_phases") is not None:
            native_initial_phases = _normalize_neutral_initial_phases(
                extra_phase_kwargs["initial_phases"],
                feed.size,
                solver_options.min_composition,
            )

    balance_matrix, total_vector, balance_names = _normalize_reactive_balances(species, balances, totals)
    reaction_defs = _normalize_reactions(species, reactions)
    reaction_matrix = np.asarray(
        [[float(reaction.stoichiometry.get(label, 0.0)) for label in species] for reaction in reaction_defs],
        dtype=float,
    )
    reaction_phase_matrix, reaction_phase_scope = _reaction_phase_stoichiometry_matrix(
        species,
        reaction_defs,
        route,
    )
    request = {
        "T": temperature,
        "P": pressure,
        "z": feed.tolist(),
        "initial_phases": native_initial_phases,
        "balance_matrix": np.asarray(balance_matrix, dtype=float).reshape(-1).tolist(),
        "balance_rows": int(balance_matrix.shape[0]),
        "total_vector": np.asarray(total_vector, dtype=float).tolist(),
        "reaction_stoichiometry": reaction_matrix.reshape(-1).tolist(),
        "reaction_rows": int(reaction_matrix.shape[0]),
        "log_equilibrium_constants": [float(reaction.log_equilibrium_constant) for reaction in reaction_defs],
        "reaction_standard_states": [reaction.convention.native_standard_state_code for reaction in reaction_defs],
        "reaction_phase_stoichiometry": None
        if reaction_phase_matrix is None
        else reaction_phase_matrix.reshape(-1).tolist(),
        "options": _options_to_native_dict(solver_options),
    }
    try:
        payload = _core._solve_reactive_phase_equilibrium_native(mixture._native, request)
    except _core.NativeValueError as exc:
        raise InputError(str(exc)) from exc
    except _core.NativeSolutionError as exc:
        message = str(exc.args[0]) if getattr(exc, "args", ()) else str(exc)
        diagnostics = exc.args[1] if len(getattr(exc, "args", ())) > 1 and isinstance(exc.args[1], dict) else None
        raise SolutionError(message, diagnostics) from exc
    result = _native_result_from_payload(payload)
    assert isinstance(result, EquilibriumResult)
    diagnostics = _reactive_phase_result_diagnostics(
        result,
        route=route,
        feed=feed,
        feed_diagnostics=feed_diagnostics,
        balance_names=balance_names,
        reactions=reaction_defs,
        reaction_matrix=reaction_matrix,
        reaction_phase_scope=reaction_phase_scope,
        reactive_options=reactive_options,
        charges=charges,
    )
    return EquilibriumResult(
        backend=result.backend,
        problem_kind=result.problem_kind,
        phases=result.phases,
        stable=result.stable,
        split_detected=result.split_detected,
        diagnostics=diagnostics,
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
        if kwargs.get("solvent_feed") is not None or kwargs.get("salt_molality") is not None:
            return "electrolyte_lle"
        charges = np.asarray(getattr(mixture, "parameters", {}).get("z", []), dtype=float).flatten()
        if charges.size == int(getattr(mixture, "ncomp", 0)) and np.any(np.abs(charges) > 1.0e-12):
            return "electrolyte_lle"
        return "lle_flash"
    if token not in {"lle_flash", "electrolyte_lle"}:
        raise InputError(
            "ReactivePhaseEquilibriumProblem production solves currently support phase_kind='lle_flash' "
            "or phase_kind='electrolyte_lle'. Use reactive_staged_equilibrium for explicit staged compatibility routes."
        )
    return token


def _reject_reactive_phase_kwargs(phase_kwargs: Mapping[str, Any], route: str) -> None:
    allowed = {"initial_phases", "solvent_feed", "salt_molality"}
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


def _reactive_phase_result_diagnostics(
    result: EquilibriumResult,
    *,
    route: str,
    feed: np.ndarray,
    feed_diagnostics: Mapping[str, Any],
    balance_names: list[str],
    reactions: list[Any],
    reaction_matrix: np.ndarray,
    reaction_phase_scope: str,
    reactive_options: Any,
    charges: np.ndarray | None,
) -> dict[str, Any]:
    diagnostics = dict(result.diagnostics)
    reaction_names = [reaction.name or f"reaction_{index}" for index, reaction in enumerate(reactions)]
    phase_compositions = {phase.label: phase.composition.tolist() for phase in result.phases}
    phase_amounts = {phase.label: float(phase.phase_fraction) for phase in result.phases}
    overall = np.zeros_like(feed, dtype=float)
    for phase in result.phases:
        overall += float(phase.phase_fraction) * np.asarray(phase.composition, dtype=float)
    element_residual = list(diagnostics.get("element_balance_residual", []))
    phase1_reaction = list(diagnostics.get("reaction_residual_phase1", []))
    phase2_reaction = list(diagnostics.get("reaction_residual_phase2", []))
    cross_phase_reaction = list(diagnostics.get("reaction_residual_cross_phase", []))
    named_reaction_residuals: dict[str, float] = {}
    for index, name in enumerate(reaction_names):
        if index < len(cross_phase_reaction):
            named_reaction_residuals[name] = abs(float(cross_phase_reaction[index]))
        else:
            values = []
            if index < len(phase1_reaction):
                values.append(float(phase1_reaction[index]))
            if index < len(phase2_reaction):
                values.append(float(phase2_reaction[index]))
            named_reaction_residuals[name] = max((abs(value) for value in values), default=0.0)
    ionic_residual = list(diagnostics.get("ionic_equilibrium_residual", []))
    diagnostics.update(
        {
            "equilibrium_route": "reactive_phase_equilibrium",
            "phase_kind": route,
            "phase_problem_kind": route,
            "reactive_phase_method": "native_coupled_reactive_phase_equilibrium",
            "reactive_workflow_class": "coupled_native",
            "coupling_level": "single_native_residual_state",
            "production_route": "native_coupled_reactive_phase_equilibrium",
            "staged_route_used": False,
            "reaction_and_phase_residuals_share_state": True,
            "reaction_phase_scope": str(diagnostics.get("reaction_phase_scope", reaction_phase_scope)),
            "phase_tagged_reaction_stoichiometry": reaction_phase_scope == "phase_tagged_cross_phase",
            "reaction_constant_policy": "fixed_literature_constants_first",
            "reaction_constant_sources": {
                name: str(reaction.metadata.get("constant_source", "unspecified"))
                for name, reaction in zip(reaction_names, reactions)
            },
            "reaction_constant_conventions": {
                name: reaction.convention.to_dict() for name, reaction in zip(reaction_names, reactions)
            },
            "reaction_coordinates": {
                "status": "implicit_in_coupled_phase_amounts",
                "reaction_count": len(reaction_names),
                "named_reactions": reaction_names,
            },
            "element_balance_residuals": {name: float(value) for name, value in zip(balance_names, element_residual)},
            "reaction_equilibrium_residuals": named_reaction_residuals,
            "phase_compositions": phase_compositions,
            "phase_amounts": phase_amounts,
            "phase_fraction_sum": float(sum(phase_amounts.values())),
            "overall_composition": overall.tolist(),
            "feed_composition": feed.tolist(),
            "reaction_extents": _named_reaction_extents(feed, overall, reaction_matrix, reaction_names),
            "reactive_options": {
                "max_iterations": int(reactive_options.max_iterations),
                "tolerance": float(reactive_options.tolerance),
                "min_mole_fraction": float(reactive_options.min_mole_fraction),
            },
        }
    )
    diagnostics.setdefault("material_balance_norm", float(diagnostics.get("element_balance_norm", 0.0)))
    diagnostics.setdefault("ionic_equilibrium_residual_norm", _max_abs(ionic_residual))
    if charges is not None:
        diagnostics.setdefault("phase_charge_balance", _phase_charge_balance_from_phases(feed, result.phases, charges))
    diagnostics.update(dict(feed_diagnostics))
    return _json_like(diagnostics)


def _named_reaction_extents(
    feed: np.ndarray,
    overall: np.ndarray,
    reaction_matrix: np.ndarray,
    reaction_names: list[str],
) -> dict[str, float]:
    if reaction_matrix.size == 0 or not reaction_names:
        return {}
    try:
        delta = overall - feed
        gram = reaction_matrix @ reaction_matrix.T
        rhs = reaction_matrix @ delta
        extents = np.linalg.solve(gram, rhs)
    except np.linalg.LinAlgError:
        return {name: float("nan") for name in reaction_names}
    return {name: float(value) for name, value in zip(reaction_names, extents)}


def _max_abs(values: Any) -> float:
    try:
        return max((abs(float(value)) for value in values), default=0.0)
    except TypeError:
        return 0.0


def initial_phases_from_result(result: EquilibriumResult) -> dict[str, object]:
    """Build electrolyte LLE ``initial_phases`` from an accepted aq/org result."""
    phases = {phase.label: phase for phase in result.phases}
    if "aq" not in phases or "org" not in phases:
        raise InputError("initial_phases_from_result requires an electrolyte LLE result with aq and org phases.")
    return {
        "aq": phases["aq"].composition.copy(),
        "org": phases["org"].composition.copy(),
        "phase_fraction": float(phases["org"].phase_fraction),
    }


def bubble_p(mixture: Any, *, T: float, x: Any, options: EquilibriumOptions | None = None) -> EquilibriumResult:
    """Solve a neutral bubble pressure at fixed liquid composition and temperature."""
    opts = _normalize_options(options)
    _normalize_feed(x, mixture.ncomp, opts.min_composition, "bubble_p")
    _reject_ion_containing_mixture(mixture)
    _positive_scalar(T, "T", "bubble_p")
    _raise_native_ipopt_equilibrium_required("bubble_p")


def dew_p(mixture: Any, *, T: float, y: Any, options: EquilibriumOptions | None = None) -> EquilibriumResult:
    """Solve a neutral dew pressure at fixed vapor composition and temperature."""
    opts = _normalize_options(options)
    _normalize_feed(y, mixture.ncomp, opts.min_composition, "dew_p")
    _reject_ion_containing_mixture(mixture)
    _positive_scalar(T, "T", "dew_p")
    _raise_native_ipopt_equilibrium_required("dew_p")


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
    """Solve a neutral TP flash through the native C++ equilibrium backend."""
    opts = _normalize_options(options)
    feed = _normalize_feed(z, mixture.ncomp, opts.min_composition, "tp_flash")
    _reject_ion_containing_mixture(mixture)
    temperature = _positive_scalar(T, "T", "tp_flash")
    pressure = _positive_scalar(P, "P", "tp_flash")
    result = _call_native_equilibrium(mixture, kind="tp_flash", T=temperature, P=pressure, z=feed, options=opts)
    assert isinstance(result, EquilibriumResult)
    return result


def lle_flash(
    mixture: Any,
    *,
    T: float,
    P: float,
    z: Any,
    options: EquilibriumOptions | None = None,
    initial_phases: Any = None,
) -> EquilibriumResult:
    """Solve a neutral LLE flash through the native C++ equilibrium backend."""
    opts = _normalize_options(options)
    feed = _normalize_feed(z, mixture.ncomp, opts.min_composition, "lle_flash")
    _reject_ion_containing_mixture(mixture)
    temperature = _positive_scalar(T, "T", "lle_flash")
    pressure = _positive_scalar(P, "P", "lle_flash")
    native_initial_phases = None
    if initial_phases is not None:
        native_initial_phases = _normalize_neutral_initial_phases(initial_phases, feed.size, opts.min_composition)
    result = _call_native_equilibrium(
        mixture,
        kind="lle_flash",
        T=temperature,
        P=pressure,
        z=feed,
        options=opts,
        initial_phases=native_initial_phases,
    )
    assert isinstance(result, EquilibriumResult)
    diagnostics = dict(result.diagnostics)
    if (
        result.split_detected is False
        and diagnostics.get("solution_accepted") is False
        and diagnostics.get("stability_stable") is False
        and "initial liquid phases are compositionally identical" not in str(diagnostics.get("message", ""))
    ):
        seed = diagnostics.get("seed_name", "unknown")
        reason = diagnostics.get("point_solver_message") or diagnostics.get("message") or "not accepted"
        raise SolutionError(f"neutral LLE flash did not converge; best_seed={seed}; {reason}", _json_like(diagnostics))
    return result


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
    """Run native C++ neutral tangent-plane-distance stability analysis."""
    opts = _normalize_options(options)
    feed = _normalize_feed(z, mixture.ncomp, opts.min_composition, "stability")
    _reject_ion_containing_mixture(mixture)
    temperature = _positive_scalar(T, "T", "stability")
    pressure = _positive_scalar(P, "P", "stability")
    result = _call_native_equilibrium(
        mixture,
        kind="stability",
        T=temperature,
        P=pressure,
        z=feed,
        options=opts,
        parent_phase=parent_phase,
        trial_phases=trial_phases,
    )
    assert isinstance(result, StabilityResult)
    return result


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
    """Run native C++ electrolyte transformed-basis stability analysis."""
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
    basis_payload = _electrolyte_formula_basis(mixture, feed, feed_diagnostics)
    result = _call_native_equilibrium(
        mixture,
        kind="electrolyte_stability",
        T=temperature,
        P=pressure,
        z=feed,
        options=opts,
    )
    assert isinstance(result, StabilityResult)
    diagnostics = dict(result.diagnostics)
    diagnostics.setdefault("algorithm_reference", dict(_ASCANI_2022_REFERENCE))
    diagnostics.setdefault("basis_rank", int(basis_payload["basis_rank"]))
    diagnostics.setdefault("e_matrix", np.asarray(basis_payload["e_matrix"], dtype=float).tolist())
    diagnostics.setdefault("salt_pairs", [dict(pair) for pair in basis_payload["salt_pairs"]])
    diagnostics.update(feed_diagnostics)
    return StabilityResult(
        backend=result.backend,
        problem_kind=result.problem_kind,
        stable=result.stable,
        min_tpd=result.min_tpd,
        parent_phase=result.parent_phase,
        trial_phase=result.trial_phase,
        trial_composition=result.trial_composition,
        trials=result.trials,
        diagnostics=diagnostics,
    )


def electrolyte_lle_flash_native(
    mixture: Any,
    *,
    T: float,
    P: float,
    z: Any = None,
    solvent_feed: Any = None,
    salt_molality: Any = None,
    initial_phases: Any = None,
    options: EquilibriumOptions | None = None,
) -> EquilibriumResult:
    """Run native C++ electrolyte LLE for native-ready request forms."""
    opts = _normalize_options(options)
    _require_ion_containing_mixture(mixture, "electrolyte_lle")
    temperature = _positive_scalar(T, "T", "electrolyte_lle")
    pressure = _positive_scalar(P, "P", "electrolyte_lle")
    feed, feed_diagnostics = _normalize_electrolyte_feed(
        mixture,
        z=z,
        solvent_feed=solvent_feed,
        salt_molality=salt_molality,
        options=opts,
    )
    charges = _mixture_charges(mixture)
    _require_charge_neutral(feed, charges, "electrolyte_lle feed")
    basis_payload = _electrolyte_formula_basis(mixture, feed, feed_diagnostics)
    native_initial_phases = None
    if initial_phases is not None:
        seed = _electrolyte_initial_phase_seed(mixture, feed, basis_payload, initial_phases, opts)
        native_initial_phases = {
            "aq": np.asarray(initial_phases["aq"], dtype=float).flatten().tolist(),
            "org": np.asarray(initial_phases["org"], dtype=float).flatten().tolist(),
            "phase_fraction": float(initial_phases["phase_fraction"]),
        }
        _ = seed
    if opts.solver_backend == "ipopt":
        _raise_native_ipopt_not_routed("electrolyte_lle")
    try:
        result = _call_native_equilibrium(
            mixture,
            kind="electrolyte_lle",
            T=temperature,
            P=pressure,
            z=feed,
            options=opts,
            initial_phases=native_initial_phases,
        )
    except SolutionError as exc:
        if isinstance(getattr(exc, "diagnostics", None), dict):
            diagnostics = dict(exc.diagnostics)
            diagnostics.setdefault("algorithm_reference", dict(_ASCANI_2022_REFERENCE))
            diagnostics.setdefault("basis_rank", int(basis_payload["basis_rank"]))
            diagnostics.setdefault("e_matrix", np.asarray(basis_payload["e_matrix"], dtype=float).tolist())
            diagnostics.setdefault("salt_pairs", [dict(pair) for pair in basis_payload["salt_pairs"]])
            diagnostics.update(feed_diagnostics)
            raise SolutionError(exc.message, _json_like(diagnostics)) from exc
        diagnostics = {
            "phase_equilibrium_model": "electrolyte_lle_v5_native_charge_constrained_solve",
            "equilibrium_route": "electrolyte_lle",
            "route_reason": "ion-containing mixture",
            "variable_model": "ascani_transformed_salt_pairs",
            "basis_rank": int(basis_payload["basis_rank"]),
            "e_matrix": np.asarray(basis_payload["e_matrix"], dtype=float).tolist(),
            "salt_pairs": [dict(pair) for pair in basis_payload["salt_pairs"]],
            "algorithm_reference": dict(_ASCANI_2022_REFERENCE),
            "feed_composition": feed.tolist(),
            "stability_analysis": "electrolyte_tpd",
            "stability_checked": True,
            "solver_backend": "ceres",
            "selected_solver_backend": "ceres",
            "solver_method": "ceres_trust_region_residual_solve",
            "solver_language": "c++",
            "native_entrypoint": "_solve_equilibrium_native",
            "jacobian_backend": "cppad_implicit",
            "derivative_backend": "cppad_implicit",
            "jacobian_available": True,
            "derivative_available": True,
            "tpd_method": "native_tpd_global_search",
            "gibbs_seed_method": "native_golden_section",
            "acceptance_gate": "predictive_solve_failed",
            "solver_residual_norm": 1.0,
            "best_failure_reason": str(exc),
        }
        diagnostics.update(feed_diagnostics)
        diagnostics = _diagnostics_with_options(options=opts, diagnostics=diagnostics)
        if diagnostics is not None:
            diagnostics = _normalize_derivative_diagnostics(
                diagnostics,
                problem_kind="electrolyte_lle",
                phase_count=0,
            )
        raise SolutionError("electrolyte LLE flash did not converge", _json_like(diagnostics)) from exc
    assert isinstance(result, EquilibriumResult)
    diagnostics = dict(result.diagnostics)
    diagnostics.setdefault("algorithm_reference", dict(_ASCANI_2022_REFERENCE))
    diagnostics.setdefault("basis_rank", int(basis_payload["basis_rank"]))
    diagnostics.setdefault("e_matrix", np.asarray(basis_payload["e_matrix"], dtype=float).tolist())
    diagnostics.setdefault("salt_pairs", [dict(pair) for pair in basis_payload["salt_pairs"]])
    diagnostics.setdefault("phase_charge_balance", _phase_charge_balance_from_phases(feed, result.phases, charges))
    diagnostics.update(feed_diagnostics)
    diagnostics = _normalize_derivative_diagnostics(
        diagnostics,
        problem_kind=result.problem_kind,
        phase_count=len(result.phases),
    )
    return EquilibriumResult(
        backend=result.backend,
        problem_kind=result.problem_kind,
        phases=result.phases,
        stable=result.stable,
        split_detected=result.split_detected,
        diagnostics=diagnostics,
    )
