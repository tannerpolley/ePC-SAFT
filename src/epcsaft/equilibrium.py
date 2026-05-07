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


@dataclass(frozen=True, slots=True)
class EquilibriumOptions:
    """Numerical controls for equilibrium solvers."""

    max_iterations: int = 180
    tolerance: float = 1.0e-6
    damping: float = 0.5
    min_composition: float = 1.0e-12
    include_phase_diagnostics: bool = False
    stability_precheck: bool = True
    legacy_candidate_mode: str = "auto"
    legacy_candidate_residual_tolerance: float = 2.0e-1
    legacy_candidate_split_tolerance: float = 1.0e-4
    legacy_candidate_max_iterations: int = 80
    ignored_legacy_options: tuple[str, ...] = ()
    density_diagnostics: Literal["auto", "off", "full"] = "auto"
    experimental_coupled_density_lle: bool = False
    jacobian_backend: Literal["auto", "autodiff", "finite_difference"] = "auto"
    solver_backend: Literal["auto", "newton", "ipopt"] = "auto"
    hessian_strategy: Literal["gauss_newton", "lbfgs"] = "gauss_newton"


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
        z_cat = int(round(abs(float(charges[cation_i]))))
        z_an = int(round(abs(float(charges[anion_i]))))
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
            "split_tol": "legacy_candidate_split_tolerance",
            "solver_accept_norm": "legacy_candidate_residual_tolerance",
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
            "damping",
            "min_composition",
            "include_phase_diagnostics",
            "stability_precheck",
            "legacy_candidate_mode",
            "legacy_candidate_residual_tolerance",
            "legacy_candidate_split_tolerance",
            "legacy_candidate_max_iterations",
            "density_diagnostics",
            "experimental_coupled_density_lle",
            "jacobian_backend",
            "solver_backend",
            "hessian_strategy",
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
    damping = _finite_float_option(options.damping, "damping")
    if not (0.0 < damping <= 1.0):
        raise InputError("options.damping must be > 0 and <= 1.")
    min_composition = _finite_float_option(options.min_composition, "min_composition")
    if min_composition <= 0.0:
        raise InputError("options.min_composition must be positive.")
    if not isinstance(options.include_phase_diagnostics, bool):
        raise InputError("options.include_phase_diagnostics must be boolean.")
    if not isinstance(options.stability_precheck, bool):
        raise InputError("options.stability_precheck must be boolean.")
    legacy_candidate_mode = str(options.legacy_candidate_mode).strip().lower()
    if legacy_candidate_mode not in {"auto", "off"}:
        raise InputError("options.legacy_candidate_mode must be 'auto' or 'off'.")
    legacy_candidate_residual_tolerance = _finite_float_option(
        options.legacy_candidate_residual_tolerance, "legacy_candidate_residual_tolerance"
    )
    if legacy_candidate_residual_tolerance <= 0.0:
        raise InputError("options.legacy_candidate_residual_tolerance must be positive.")
    legacy_candidate_split_tolerance = _finite_float_option(
        options.legacy_candidate_split_tolerance, "legacy_candidate_split_tolerance"
    )
    if legacy_candidate_split_tolerance <= 0.0:
        raise InputError("options.legacy_candidate_split_tolerance must be positive.")
    if isinstance(options.legacy_candidate_max_iterations, bool) or not isinstance(
        options.legacy_candidate_max_iterations, Integral
    ):
        raise InputError("options.legacy_candidate_max_iterations must be an integer greater than zero.")
    legacy_candidate_max_iterations = int(options.legacy_candidate_max_iterations)
    if legacy_candidate_max_iterations <= 0:
        raise InputError("options.legacy_candidate_max_iterations must be an integer greater than zero.")
    ignored_legacy_options = tuple(str(item) for item in options.ignored_legacy_options)
    density_diagnostics = str(options.density_diagnostics).strip().lower()
    if density_diagnostics not in {"auto", "off", "full"}:
        raise InputError("options.density_diagnostics must be 'auto', 'off', or 'full'.")
    if not isinstance(options.experimental_coupled_density_lle, bool):
        raise InputError("options.experimental_coupled_density_lle must be boolean.")
    jacobian_backend = str(options.jacobian_backend).strip().lower()
    aliases = {"numerical": "finite_difference", "fd": "finite_difference"}
    jacobian_backend = aliases.get(jacobian_backend, jacobian_backend)
    if jacobian_backend not in {"auto", "autodiff", "finite_difference"}:
        raise InputError("options.jacobian_backend must be 'auto', 'autodiff', or 'finite_difference'.")
    solver_backend = str(options.solver_backend).strip().lower()
    if solver_backend not in {"auto", "newton", "ipopt"}:
        raise InputError("options.solver_backend must be 'auto', 'newton', or 'ipopt'.")
    hessian_strategy = str(options.hessian_strategy).strip().lower()
    hessian_aliases = {"gn": "gauss_newton", "gauss-newton": "gauss_newton", "bfgs": "lbfgs"}
    hessian_strategy = hessian_aliases.get(hessian_strategy, hessian_strategy)
    if hessian_strategy not in {"gauss_newton", "lbfgs"}:
        raise InputError("options.hessian_strategy must be 'gauss_newton' or 'lbfgs'.")
    return EquilibriumOptions(
        max_iterations=max_iterations,
        tolerance=tolerance,
        damping=damping,
        min_composition=min_composition,
        include_phase_diagnostics=options.include_phase_diagnostics,
        stability_precheck=options.stability_precheck,
        legacy_candidate_mode=legacy_candidate_mode,
        legacy_candidate_residual_tolerance=legacy_candidate_residual_tolerance,
        legacy_candidate_split_tolerance=legacy_candidate_split_tolerance,
        legacy_candidate_max_iterations=legacy_candidate_max_iterations,
        ignored_legacy_options=ignored_legacy_options,
        density_diagnostics=density_diagnostics,  # type: ignore[arg-type]
        experimental_coupled_density_lle=options.experimental_coupled_density_lle,
        jacobian_backend=jacobian_backend,  # type: ignore[arg-type]
        solver_backend=solver_backend,  # type: ignore[arg-type]
        hessian_strategy=hessian_strategy,  # type: ignore[arg-type]
    )


def _finite_float_option(value: Any, label: str) -> float:
    if isinstance(value, bool) or not isinstance(value, Real):
        raise InputError("options.{} must be a finite real number.".format(label))
    out = float(value)
    if not np.isfinite(out):
        raise InputError("options.{} must be finite.".format(label))
    return out


def _normalize_feed(z: Any, ncomp: int, min_composition: float, kind: str) -> np.ndarray:
    if z is None:
        raise InputError("z is required for kind='{}'.".format(kind))
    feed = np.asarray(z, dtype=float).flatten()
    if feed.size != int(ncomp):
        raise InputError(
            "Feed composition length ({}) must match mixture component count ({}).".format(feed.size, ncomp)
        )
    if not np.all(np.isfinite(feed)):
        raise InputError("Feed composition z must contain only finite values.")
    if np.any(feed < 0.0):
        raise InputError("Feed composition z must be non-negative.")
    total = float(np.sum(feed))
    if total <= 0.0:
        raise InputError("Feed composition z must have a positive sum.")
    feed = feed / total
    if np.any(feed < min_composition):
        raise InputError("{} requires each feed composition entry to be >= min_composition.".format(kind))
    return feed


def _positive_scalar(value: Any, label: str, kind: str) -> float:
    if value is None:
        raise InputError("{} is required for kind='{}'.".format(label, kind))
    out = float(value)
    if not np.isfinite(out) or out <= 0.0:
        raise InputError("{} must be a positive finite scalar.".format(label))
    return out


def _reject_ion_containing_mixture(mixture: Any) -> None:
    z = np.asarray(mixture.parameters.get("z", []), dtype=float).flatten()
    if z.size and np.any(np.abs(z) > 1.0e-12):
        raise InputError("Neutral equilibrium does not support ion-containing mixtures.")


def _require_ion_containing_mixture(mixture: Any, kind: str) -> None:
    charges = _mixture_charges(mixture)
    if not np.any(np.abs(charges) > 1.0e-12):
        raise InputError("{} requires an ion-containing mixture.".format(kind))


def _mixture_charges(mixture: Any) -> np.ndarray:
    charges = np.asarray(mixture.parameters.get("z", []), dtype=float).flatten()
    if charges.size != int(mixture.ncomp):
        raise InputError("mixture parameters must include one charge value per species in params['z'].")
    return charges


def _require_charge_neutral(composition: np.ndarray, charges: np.ndarray, label: str) -> None:
    charge = float(np.dot(np.asarray(composition, dtype=float), np.asarray(charges, dtype=float)))
    if abs(charge) > 1.0e-10:
        raise InputError("{} must be charge neutral; charge balance is {}.".format(label, charge))


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
                raise InputError("solvent_feed species '{}' is not present in the mixture.".format(label)) from exc
            if index not in neutral_indices:
                raise InputError("solvent_feed species '{}' is ionic; expected neutral solvents only.".format(label))
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
        raise InputError("Could not uniquely map salt_molality key '{}' onto mixture ions.".format(salt_label))
    return matches[0]


def _salt_label_token(label: Any) -> str:
    return "".join(ch for ch in str(label) if ch.isalnum()).lower()


def _ion_stem(label: str, charge: float | None = None) -> str:
    text = str(label)
    stripped = text.replace("+", "").replace("-", "")
    if charge is not None and ("+" in text or "-" in text):
        charge_int = int(round(abs(float(charge))))
        suffix = str(charge_int)
        if charge_int > 1 and stripped.endswith(suffix):
            stripped = stripped[: -len(suffix)]
    return stripped


def _salt_stoichiometry(cation_charge: float, anion_charge: float) -> tuple[int, int]:
    cation_int = int(round(abs(float(cation_charge))))
    anion_int = int(round(abs(float(anion_charge))))
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
    except Exception as exc:
        raise SolutionError("Failed to construct {} phase during {}: {}".format(label, context, exc)) from exc
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
        raise InputError("{} must be None, 'liq', or 'vap'.".format(label))
    return token


def _normalize_initial_phase(value: Any, ncomp: int, min_composition: float, label: str) -> np.ndarray:
    composition = np.asarray(value, dtype=float).flatten()
    if composition.size != int(ncomp):
        raise InputError(
            "initial_phases {} length ({}) must match mixture component count ({}).".format(
                label, composition.size, ncomp
            )
        )
    if not np.all(np.isfinite(composition)):
        raise InputError("initial_phases {} must contain only finite values.".format(label))
    if np.any(composition < 0.0):
        raise InputError("initial_phases {} must be non-negative.".format(label))
    total = float(np.sum(composition))
    if total <= 0.0:
        raise InputError("initial_phases {} must have a positive sum.".format(label))
    composition = composition / total
    if np.any(composition < min_composition):
        raise InputError("initial_phases {} entries must be >= min_composition.".format(label))
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
        "damping": float(options.damping),
        "min_composition": float(options.min_composition),
        "include_phase_diagnostics": bool(options.include_phase_diagnostics),
        "stability_precheck": bool(options.stability_precheck),
        "density_diagnostics": str(options.density_diagnostics),
        "experimental_coupled_density_lle": bool(options.experimental_coupled_density_lle),
        "jacobian_backend": str(options.jacobian_backend),
        "solver_backend": str(options.solver_backend),
        "hessian_strategy": str(options.hessian_strategy),
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


def _native_result_from_payload(payload: dict[str, Any]) -> EquilibriumResult | StabilityResult:
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
        from .ipopt_backend import unsupported_ipopt_route

        unsupported_ipopt_route(kind)
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
        diagnostics = _diagnostics_with_legacy_candidate(
            mixture,
            kind=kind,
            T=float(T),
            P=float(P),
            feed=np.asarray(z, dtype=float).flatten(),
            options=options,
            diagnostics=diagnostics,
        )
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
    diagnostics.setdefault("legacy_candidate_mode", str(options.legacy_candidate_mode))
    diagnostics.setdefault("legacy_candidate_residual_tolerance", float(options.legacy_candidate_residual_tolerance))
    diagnostics.setdefault("legacy_candidate_split_tolerance", float(options.legacy_candidate_split_tolerance))
    diagnostics.setdefault("ignored_legacy_options", list(options.ignored_legacy_options))
    diagnostics.setdefault("density_diagnostics_mode", str(options.density_diagnostics))
    diagnostics.setdefault("experimental_coupled_density_lle", bool(options.experimental_coupled_density_lle))
    diagnostics.setdefault("requested_solver_backend", str(options.solver_backend))
    diagnostics.setdefault("requested_hessian_strategy", str(options.hessian_strategy))
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
    diagnostics.setdefault("density_fallback_used", False)
    diagnostics.setdefault("density_fallback_rejected_reason", "")
    diagnostics.setdefault("density_warm_start_source", "")
    diagnostics.setdefault("density_validity_gate", "not_evaluated")


def _diagnostics_with_legacy_candidate(
    mixture: Any,
    *,
    kind: str,
    T: float,
    P: float,
    feed: np.ndarray,
    options: EquilibriumOptions,
    diagnostics: dict[str, Any] | None,
) -> dict[str, Any] | None:
    if diagnostics is None:
        return None
    out = dict(diagnostics)
    _add_legacy_option_diagnostics(out, options)
    if kind not in {"electrolyte_lle", "electrolyte_lle_flash"}:
        return out
    if str(options.legacy_candidate_mode) != "auto":
        out.setdefault("legacy_candidate_found", False)
        out.setdefault(
            "legacy_candidate_message", "legacy candidate fallback disabled by options.legacy_candidate_mode"
        )
        return out
    min_tpd = out.get("stability_min_tpd", out.get("min_tpd"))
    try:
        unstable = float(min_tpd) < -max(float(options.tolerance), 1.0e-8)
    except (TypeError, ValueError):
        unstable = out.get("acceptance_gate") == "predictive_solve_failed"
    collapsed_or_failed = (
        out.get("best_failure_reason") == "candidate collapsed to one phase"
        or out.get("acceptance_gate") == "predictive_solve_failed"
    )
    if not unstable or not collapsed_or_failed:
        out.setdefault("legacy_candidate_found", False)
        out.setdefault("legacy_candidate_message", "legacy candidate fallback was not triggered")
        return out
    out.setdefault("legacy_candidate_found", False)
    out.setdefault(
        "legacy_candidate_message",
        "legacy candidate fallback disabled; Python equilibrium candidates are not exposed.",
    )
    return out


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
    liquid = _normalize_feed(x, mixture.ncomp, opts.min_composition, "bubble_p")
    _reject_ion_containing_mixture(mixture)
    temperature = _positive_scalar(T, "T", "bubble_p")
    return _neutral_bubble_dew_outer(
        mixture,
        problem_kind="bubble_p",
        fixed_name="P",
        fixed_value=temperature,
        source_composition=liquid,
        source_phase="liq",
        incipient_phase="vap",
        options=opts,
    )


def dew_p(mixture: Any, *, T: float, y: Any, options: EquilibriumOptions | None = None) -> EquilibriumResult:
    """Solve a neutral dew pressure at fixed vapor composition and temperature."""
    opts = _normalize_options(options)
    vapor = _normalize_feed(y, mixture.ncomp, opts.min_composition, "dew_p")
    _reject_ion_containing_mixture(mixture)
    temperature = _positive_scalar(T, "T", "dew_p")
    return _neutral_bubble_dew_outer(
        mixture,
        problem_kind="dew_p",
        fixed_name="P",
        fixed_value=temperature,
        source_composition=vapor,
        source_phase="vap",
        incipient_phase="liq",
        options=opts,
    )


def bubble_t(mixture: Any, *, P: float, x: Any, options: EquilibriumOptions | None = None) -> EquilibriumResult:
    """Solve a neutral bubble temperature at fixed liquid composition and pressure."""
    opts = _normalize_options(options)
    liquid = _normalize_feed(x, mixture.ncomp, opts.min_composition, "bubble_t")
    _reject_ion_containing_mixture(mixture)
    pressure = _positive_scalar(P, "P", "bubble_t")
    return _neutral_bubble_dew_outer(
        mixture,
        problem_kind="bubble_t",
        fixed_name="T",
        fixed_value=pressure,
        source_composition=liquid,
        source_phase="liq",
        incipient_phase="vap",
        options=opts,
    )


def dew_t(mixture: Any, *, P: float, y: Any, options: EquilibriumOptions | None = None) -> EquilibriumResult:
    """Solve a neutral dew temperature at fixed vapor composition and pressure."""
    opts = _normalize_options(options)
    vapor = _normalize_feed(y, mixture.ncomp, opts.min_composition, "dew_t")
    _reject_ion_containing_mixture(mixture)
    pressure = _positive_scalar(P, "P", "dew_t")
    return _neutral_bubble_dew_outer(
        mixture,
        problem_kind="dew_t",
        fixed_name="T",
        fixed_value=pressure,
        source_composition=vapor,
        source_phase="vap",
        incipient_phase="liq",
        options=opts,
    )


def _neutral_bubble_dew_outer(
    mixture: Any,
    *,
    problem_kind: str,
    fixed_name: str,
    fixed_value: float,
    source_composition: np.ndarray,
    source_phase: str,
    incipient_phase: str,
    options: EquilibriumOptions,
) -> EquilibriumResult:
    solve_pressure = fixed_name == "P"
    grid = np.geomspace(1.0, 1.0e8, 81) if solve_pressure else np.linspace(120.0, 700.0, 117)
    evaluations: list[dict[str, Any]] = []
    failures: list[dict[str, str]] = []
    for candidate in grid:
        try:
            evaluation = _neutral_bubble_dew_evaluate(
                mixture,
                problem_kind=problem_kind,
                variable=float(candidate),
                fixed_value=fixed_value,
                source_composition=source_composition,
                source_phase=source_phase,
                incipient_phase=incipient_phase,
                options=options,
            )
        except Exception as exc:
            failures.append({"variable": str(float(candidate)), "message": str(exc)})
            continue
        evaluations.append(evaluation)
        if abs(float(evaluation["residual"])) <= options.tolerance:
            return _neutral_bubble_dew_result(mixture, problem_kind, evaluation, options, evaluations, failures)
    bracket: tuple[dict[str, Any], dict[str, Any]] | None = None
    for left, right in zip(evaluations, evaluations[1:]):
        if float(left["residual"]) * float(right["residual"]) <= 0.0:
            bracket = (left, right)
            break
    if bracket is None:
        diagnostics = {
            "message": "failed to bracket neutral bubble/dew scalar residual",
            "problem_kind": problem_kind,
            "residual_samples": [{"variable": item["variable"], "residual": item["residual"]} for item in evaluations],
            "state_failures": failures[:10],
        }
        raise SolutionError("neutral {} did not bracket a scalar root".format(problem_kind), diagnostics)
    left, right = bracket
    history = [left, right]
    best = min((left, right), key=lambda item: abs(float(item["residual"])))
    for _ in range(options.max_iterations):
        if solve_pressure:
            midpoint = float(np.exp(0.5 * (np.log(float(left["variable"])) + np.log(float(right["variable"])))))
        else:
            midpoint = 0.5 * (float(left["variable"]) + float(right["variable"]))
        current = _neutral_bubble_dew_evaluate(
            mixture,
            problem_kind=problem_kind,
            variable=midpoint,
            fixed_value=fixed_value,
            source_composition=source_composition,
            source_phase=source_phase,
            incipient_phase=incipient_phase,
            options=options,
        )
        history.append(current)
        if abs(float(current["residual"])) < abs(float(best["residual"])):
            best = current
        if abs(float(current["residual"])) <= options.tolerance:
            return _neutral_bubble_dew_result(mixture, problem_kind, current, options, history, failures)
        if float(left["residual"]) * float(current["residual"]) <= 0.0:
            right = current
        else:
            left = current
    if abs(float(best["residual"])) <= 10.0 * options.tolerance:
        return _neutral_bubble_dew_result(mixture, problem_kind, best, options, history, failures)
    diagnostics = {
        "message": "neutral bubble/dew scalar solve reached max_iterations",
        "problem_kind": problem_kind,
        "best_variable": best["variable"],
        "best_residual": best["residual"],
        "residual_history": [{"variable": item["variable"], "residual": item["residual"]} for item in history],
        "state_failures": failures[:10],
    }
    raise SolutionError("neutral {} did not converge".format(problem_kind), diagnostics)


def _neutral_bubble_dew_evaluate(
    mixture: Any,
    *,
    problem_kind: str,
    variable: float,
    fixed_value: float,
    source_composition: np.ndarray,
    source_phase: str,
    incipient_phase: str,
    options: EquilibriumOptions,
) -> dict[str, Any]:
    if problem_kind.endswith("_p"):
        temperature = fixed_value
        pressure = variable
    else:
        temperature = variable
        pressure = fixed_value
    incipient = np.array(source_composition, dtype=float, copy=True)
    last_residual = float("inf")
    for inner_iteration in range(max(1, options.max_iterations)):
        source_state = _phase_state(
            mixture, temperature, pressure, source_composition, source_phase, options, problem_kind
        )
        incipient_state = _phase_state(
            mixture, temperature, pressure, incipient, incipient_phase, options, problem_kind
        )
        ln_k = np.asarray(source_state["ln_phi"], dtype=float) - np.asarray(incipient_state["ln_phi"], dtype=float)
        k_values = np.exp(np.clip(ln_k, -100.0, 100.0))
        raw = source_composition * k_values
        residual = float(np.sum(raw) - 1.0)
        total = float(np.sum(raw))
        if not np.isfinite(total) or total <= 0.0:
            raise SolutionError("neutral {} produced a non-positive incipient composition sum".format(problem_kind))
        updated = np.maximum(raw / total, options.min_composition)
        updated = updated / float(np.sum(updated))
        delta = float(np.max(np.abs(updated - incipient)))
        incipient = updated
        last_residual = residual
        if delta <= max(options.tolerance, 1.0e-12):
            break
    source_state = _phase_state(mixture, temperature, pressure, source_composition, source_phase, options, problem_kind)
    incipient_state = _phase_state(mixture, temperature, pressure, incipient, incipient_phase, options, problem_kind)
    fugacity_residual = (
        np.log(np.maximum(incipient, options.min_composition))
        + np.asarray(incipient_state["ln_phi"], dtype=float)
        - np.log(np.maximum(source_composition, options.min_composition))
        - np.asarray(source_state["ln_phi"], dtype=float)
    )
    return {
        "variable": float(variable),
        "T": float(temperature),
        "P": float(pressure),
        "residual": float(last_residual),
        "source_composition": source_composition,
        "incipient_composition": incipient,
        "source_phase": source_phase,
        "incipient_phase": incipient_phase,
        "source_state": source_state,
        "incipient_state": incipient_state,
        "fugacity_residual": fugacity_residual,
        "fugacity_residual_norm": float(np.max(np.abs(fugacity_residual))),
        "inner_iterations": inner_iteration + 1,
    }


def _neutral_bubble_dew_result(
    mixture: Any,
    problem_kind: str,
    evaluation: dict[str, Any],
    options: EquilibriumOptions,
    history: list[dict[str, Any]],
    failures: list[dict[str, str]],
) -> EquilibriumResult:
    source_phase = str(evaluation["source_phase"])
    incipient_phase = str(evaluation["incipient_phase"])
    source = EquilibriumPhase(
        label=source_phase,
        composition=evaluation["source_composition"],
        density=evaluation["source_state"]["density"],
        temperature=evaluation["T"],
        pressure=evaluation["P"],
        phase_fraction=1.0,
        ln_fugacity_coefficient=evaluation["source_state"]["ln_phi"],
        diagnostics=evaluation["source_state"]["diagnostics"],
    )
    incipient = EquilibriumPhase(
        label=incipient_phase,
        composition=evaluation["incipient_composition"],
        density=evaluation["incipient_state"]["density"],
        temperature=evaluation["T"],
        pressure=evaluation["P"],
        phase_fraction=0.0,
        ln_fugacity_coefficient=evaluation["incipient_state"]["ln_phi"],
        diagnostics=evaluation["incipient_state"]["diagnostics"],
    )
    phases = (source, incipient) if source_phase == "liq" else (incipient, source)
    diagnostics = {
        "solver_method": "scalar_outer_composition_inner_update",
        "problem_kind": problem_kind,
        "T": float(evaluation["T"]),
        "P": float(evaluation["P"]),
        "scalar_residual": float(evaluation["residual"]),
        "fugacity_residual_norm": float(evaluation["fugacity_residual_norm"]),
        "fugacity_residual": np.asarray(evaluation["fugacity_residual"], dtype=float).tolist(),
        "inner_iterations": int(evaluation["inner_iterations"]),
        "outer_iterations": len(history),
        "residual_history": [{"variable": item["variable"], "residual": item["residual"]} for item in history],
        "state_failures": failures[:10],
        "min_composition": float(options.min_composition),
        "species": list(getattr(mixture, "species", [])),
    }
    return EquilibriumResult(
        backend="neutral_bubble_dew",
        problem_kind=problem_kind,
        phases=phases,
        stable=False,
        split_detected=False,
        diagnostics=diagnostics,
    )


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
        from .ipopt_backend import solve_electrolyte_lle_ipopt

        return solve_electrolyte_lle_ipopt(
            mixture=mixture,
            T=temperature,
            P=pressure,
            feed=feed,
            feed_diagnostics=feed_diagnostics,
            basis_payload=basis_payload,
            initial_phases=native_initial_phases,
            options=opts,
        )
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
            "solver_method": "native_transformed_newton",
            "solver_language": "c++",
            "native_entrypoint": "_solve_equilibrium_native",
            "tpd_method": "native_tpd_global_search",
            "gibbs_seed_method": "native_golden_section",
            "acceptance_gate": "predictive_solve_failed",
            "solver_residual_norm": 1.0,
            "best_failure_reason": str(exc),
        }
        diagnostics.update(feed_diagnostics)
        diagnostics = _diagnostics_with_legacy_candidate(
            mixture,
            kind="electrolyte_lle",
            T=temperature,
            P=pressure,
            feed=feed,
            options=opts,
            diagnostics=diagnostics,
        )
        raise SolutionError("electrolyte LLE flash did not converge", _json_like(diagnostics)) from exc
    assert isinstance(result, EquilibriumResult)
    diagnostics = dict(result.diagnostics)
    diagnostics.setdefault("algorithm_reference", dict(_ASCANI_2022_REFERENCE))
    diagnostics.setdefault("basis_rank", int(basis_payload["basis_rank"]))
    diagnostics.setdefault("e_matrix", np.asarray(basis_payload["e_matrix"], dtype=float).tolist())
    diagnostics.setdefault("salt_pairs", [dict(pair) for pair in basis_payload["salt_pairs"]])
    diagnostics.update(feed_diagnostics)
    return EquilibriumResult(
        backend=result.backend,
        problem_kind=result.problem_kind,
        phases=result.phases,
        stable=result.stable,
        split_detected=result.split_detected,
        diagnostics=diagnostics,
    )
