"""Python-first phase-equilibrium helpers built on native ePC-SAFT states."""

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


def tp_flash(
    mixture: Any, *, T: float, P: float, z: Any, options: EquilibriumOptions | None = None
) -> EquilibriumResult:
    """Solve a V1 neutral TP flash with Rachford-Rice and fugacity updates."""
    opts = _normalize_options(options)
    feed = _normalize_feed(z, mixture.ncomp, opts.min_composition, "tp_flash")
    _reject_ion_containing_mixture(mixture)
    temperature = _positive_scalar(T, "T", "tp_flash")
    pressure = _positive_scalar(P, "P", "tp_flash")

    liquid_seed = _phase_state(mixture, temperature, pressure, feed, "liq", opts, "TP flash")
    vapor_seed = _phase_state(mixture, temperature, pressure, feed, "vap", opts, "TP flash")
    stability_precheck = _equilibrium_stability_precheck(
        mixture,
        T=temperature,
        P=pressure,
        feed=feed,
        options=opts,
        parent_phase=None,
        trial_phases=("liq", "vap"),
    )
    stability_diagnostics = stability_precheck["diagnostics"]
    ln_k = liquid_seed["ln_phi"] - vapor_seed["ln_phi"]

    split, beta, no_split_message = _rachford_rice_beta(feed, np.exp(ln_k))
    if not split:
        phase_label = "liq" if beta <= 0.0 else "vap"
        state_payload = liquid_seed if phase_label == "liq" else vapor_seed
        phase = _phase_from_state(
            phase_label,
            feed,
            1.0,
            state_payload,
            opts,
        )
        diagnostics = {
            "iterations": 0,
            "fugacity_residual_norm": 0.0,
            "material_balance_error": 0.0,
            "vapor_fraction": float(beta),
            "message": no_split_message,
            "point_solver_split_detected": False,
            "point_solver_message": no_split_message,
        }
        diagnostics.update(stability_diagnostics)
        return EquilibriumResult(
            backend="neutral_vle",
            problem_kind="tp_flash",
            phases=(phase,),
            stable=_no_split_stable_from_precheck(stability_diagnostics),
            split_detected=False,
            diagnostics=diagnostics,
        )

    best: dict[str, Any] | None = None
    for iteration in range(1, opts.max_iterations + 1):
        k_values = np.exp(ln_k)
        split, beta, no_split_message = _rachford_rice_beta(feed, k_values)
        if not split:
            raise SolutionError(
                "TP flash lost its two-phase Rachford-Rice bracket after {} iterations: {}".format(
                    iteration,
                    no_split_message,
                )
            )
        x_liq, y_vap = _phase_compositions(feed, k_values, beta, opts.min_composition)
        liquid = _phase_state(mixture, temperature, pressure, x_liq, "liq", opts, "TP flash")
        vapor = _phase_state(mixture, temperature, pressure, y_vap, "vap", opts, "TP flash")
        fugacity_residual = np.log(y_vap) + vapor["ln_phi"] - np.log(x_liq) - liquid["ln_phi"]
        residual_norm = float(np.max(np.abs(fugacity_residual)))
        material_error = float(np.max(np.abs((1.0 - beta) * x_liq + beta * y_vap - feed)))
        best = {
            "iteration": iteration,
            "beta": float(beta),
            "x_liq": x_liq,
            "y_vap": y_vap,
            "liquid": liquid,
            "vapor": vapor,
            "fugacity_residual": fugacity_residual,
            "residual_norm": residual_norm,
            "material_error": material_error,
        }
        if residual_norm <= opts.tolerance and material_error <= max(opts.tolerance, 1.0e-10):
            return _two_phase_result(best, opts, "converged", stability_diagnostics)
        ln_k_target = liquid["ln_phi"] - vapor["ln_phi"]
        ln_k = (1.0 - opts.damping) * ln_k + opts.damping * ln_k_target

    assert best is not None
    raise SolutionError(
        "neutral TP flash did not converge after {} iterations; residual_norm={}, material_balance_error={}".format(
            opts.max_iterations,
            best["residual_norm"],
            best["material_error"],
        )
    )


def lle_flash(
    mixture: Any,
    *,
    T: float,
    P: float,
    z: Any,
    options: EquilibriumOptions | None = None,
    initial_phases: Any = None,
) -> EquilibriumResult:
    """Solve a V2 neutral liquid-liquid TP flash with damped finite-difference Newton updates."""
    opts = _normalize_options(options)
    feed = _normalize_feed(z, mixture.ncomp, opts.min_composition, "lle_flash")
    _reject_ion_containing_mixture(mixture)
    temperature = _positive_scalar(T, "T", "lle_flash")
    pressure = _positive_scalar(P, "P", "lle_flash")
    stability_precheck = _equilibrium_stability_precheck(
        mixture,
        T=temperature,
        P=pressure,
        feed=feed,
        options=opts,
        parent_phase="liq",
        trial_phases=("liq",),
    )
    stability_diagnostics = stability_precheck["diagnostics"]
    seeds = _initial_lle_guesses(initial_phases, feed, opts, stability_precheck)
    state_feed = _phase_state(mixture, temperature, pressure, feed, "liq", opts, "LLE flash")

    degenerate_attempts: list[dict[str, Any]] = []
    failed_attempts: list[dict[str, Any]] = []
    for attempt_count, seed in enumerate(seeds, start=1):
        attempt = _solve_lle_attempt(
            mixture,
            temperature,
            pressure,
            feed,
            opts,
            seed,
            attempt_count,
        )
        if attempt["status"] == "converged":
            return _lle_two_phase_result(
                attempt["candidate"],
                opts,
                "converged",
                seed["seed_name"],
                attempt_count,
                stability_diagnostics,
            )
        if attempt["status"] == "degenerate":
            degenerate_attempts.append(attempt)
        else:
            failed_attempts.append(attempt)

    if failed_attempts:
        best_failed = min(failed_attempts, key=lambda item: item["candidate"]["objective"])
        best = best_failed["candidate"]
        raise SolutionError(
            "neutral LLE flash did not converge after {} attempt(s); best_seed={}, reason={}, residual_norm={}, material_balance_error={}, phase_distance={}".format(
                len(seeds),
                best_failed["seed_name"],
                best_failed["message"],
                best["fugacity_residual_norm"],
                best["material_error"],
                _phase_distance(best["comp1"], best["comp2"]),
            )
        )

    best_degenerate = min(degenerate_attempts, key=lambda item: item["candidate"]["objective"])
    best = best_degenerate["candidate"]
    return _lle_no_split_result(
        feed,
        state_feed,
        opts,
        best_degenerate["message"],
        int(best["iteration"]),
        best["fugacity_residual_norm"],
        best["material_error"],
        best["beta"],
        best["comp1"],
        best["comp2"],
        best_degenerate["seed_name"],
        best_degenerate["attempt_count"],
        stability_diagnostics,
    )


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
    """Run V3 neutral tangent-plane-distance stability analysis."""
    opts = _normalize_options(options)
    feed = _normalize_feed(z, mixture.ncomp, opts.min_composition, "stability")
    _reject_ion_containing_mixture(mixture)
    temperature = _positive_scalar(T, "T", "stability")
    pressure = _positive_scalar(P, "P", "stability")
    parent_tokens = _normalize_parent_phases(parent_phase)
    trial_tokens = _normalize_trial_phases(trial_phases)
    threshold = _tpd_instability_threshold(opts)

    trials: list[StabilityTrial] = []
    seeds = _tpd_trial_seeds(feed, opts)
    for parent_token in parent_tokens:
        parent = _phase_state(mixture, temperature, pressure, feed, parent_token, opts, "TPD stability parent")
        for trial_token in trial_tokens:
            for seed in seeds:
                trials.append(
                    _solve_tpd_trial(
                        mixture,
                        temperature,
                        pressure,
                        feed,
                        parent["ln_phi"],
                        parent_token,
                        trial_token,
                        seed,
                        opts,
                        threshold,
                    )
                )

    best = min(trials, key=lambda trial: trial.tpd)
    stable = not any(trial.unstable for trial in trials)
    return StabilityResult(
        backend="neutral_tpd",
        problem_kind="stability",
        stable=stable,
        min_tpd=best.tpd,
        parent_phase=best.parent_phase,
        trial_phase=best.trial_phase,
        trial_composition=best.composition,
        trials=tuple(trials),
        diagnostics={
            "stability_analysis": "neutral_tpd",
            "parent_phases": list(parent_tokens),
            "trial_phases": list(trial_tokens),
            "trial_count": len(trials),
            "seed_count": len(seeds),
            "tpd_threshold": threshold,
            "min_seed_name": best.seed_name,
            "message": "unstable trial phase detected" if not stable else "no negative TPD trial found",
        },
    )


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
    )


def _stability_precheck_options(options: EquilibriumOptions) -> EquilibriumOptions:
    return EquilibriumOptions(
        max_iterations=min(int(options.max_iterations), 40),
        tolerance=max(float(options.tolerance), 1.0e-8),
        damping=options.damping,
        min_composition=options.min_composition,
        include_phase_diagnostics=False,
        stability_precheck=True,
        density_diagnostics=options.density_diagnostics,
        experimental_coupled_density_lle=options.experimental_coupled_density_lle,
        jacobian_backend=options.jacobian_backend,
    )


def _equilibrium_stability_precheck(
    mixture: Any,
    *,
    T: float,
    P: float,
    feed: np.ndarray,
    options: EquilibriumOptions,
    parent_phase: Any,
    trial_phases: Any,
) -> dict[str, Any]:
    if not options.stability_precheck:
        return {
            "diagnostics": {
                "stability_analysis": "not_run",
                "stability_checked": False,
                "stability_stable": None,
                "stability_message": "stability precheck skipped by options.stability_precheck=False",
            },
            "unstable_trial_composition": None,
        }
    precheck_options = _stability_precheck_options(options)
    stability_mixture = type(mixture).from_params(mixture.parameters, species=mixture.species)
    result = neutral_stability(
        stability_mixture,
        T=T,
        P=P,
        z=feed,
        options=precheck_options,
        parent_phase=parent_phase,
        trial_phases=trial_phases,
    )
    unstable_trial_count = sum(1 for trial in result.trials if trial.unstable)
    unstable_trials = [trial for trial in result.trials if trial.unstable]
    best_unstable_trial = min(unstable_trials, key=lambda trial: trial.tpd) if unstable_trials else None
    return {
        "diagnostics": {
            "stability_analysis": "neutral_tpd",
            "stability_checked": True,
            "stability_stable": bool(result.stable),
            "min_tpd": float(result.min_tpd),
            "parent_phase": result.parent_phase,
            "trial_phase": result.trial_phase,
            "trial_composition": result.trial_composition.tolist(),
            "unstable_trial_count": int(unstable_trial_count),
            "trial_count": len(result.trials),
            "stability_max_iterations": int(precheck_options.max_iterations),
            "stability_tolerance": float(precheck_options.tolerance),
        },
        "unstable_trial_composition": None if best_unstable_trial is None else best_unstable_trial.composition.copy(),
    }


def _no_split_stable_from_precheck(stability_diagnostics: dict[str, Any]) -> bool:
    if stability_diagnostics.get("stability_checked") is False:
        return False
    return bool(stability_diagnostics.get("stability_stable", True))


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


def _tpd_instability_threshold(options: EquilibriumOptions) -> float:
    return max(float(options.tolerance), 1.0e-8)


def _tpd_trial_seeds(feed: np.ndarray, options: EquilibriumOptions) -> list[dict[str, Any]]:
    seeds = [{"seed_name": "feed", "composition": np.asarray(feed, dtype=float)}]
    for index in range(int(feed.size)):
        seeds.append(
            {
                "seed_name": "component_{}_rich".format(index),
                "composition": _component_rich_lle_composition(feed.size, index, options.min_composition),
            }
        )
    return seeds


def _solve_tpd_trial(
    mixture: Any,
    T: float,
    P: float,
    feed: np.ndarray,
    parent_ln_phi: np.ndarray,
    parent_phase: str,
    trial_phase: str,
    seed: dict[str, Any],
    options: EquilibriumOptions,
    threshold: float,
) -> StabilityTrial:
    composition = _clip_normalize_composition(seed["composition"], options.min_composition)
    best_tpd = np.inf
    best_composition = composition.copy()
    best_iteration = 0
    converged = False
    max_delta = np.inf
    for iteration in range(1, options.max_iterations + 1):
        trial = _phase_state(mixture, T, P, composition, trial_phase, options, "TPD stability trial")
        tpd_value = _tpd_value(composition, trial["ln_phi"], feed, parent_ln_phi)
        if tpd_value < best_tpd:
            best_tpd = tpd_value
            best_composition = composition.copy()
            best_iteration = iteration

        target = _composition_from_log_weights(np.log(feed) + parent_ln_phi - trial["ln_phi"], options.min_composition)
        next_composition = _clip_normalize_composition(
            (1.0 - options.damping) * composition + options.damping * target,
            options.min_composition,
        )
        max_delta = float(np.max(np.abs(next_composition - composition)))
        if max_delta <= max(options.tolerance, 1.0e-10):
            converged = True
            break
        composition = next_composition

    return StabilityTrial(
        parent_phase=parent_phase,
        trial_phase=trial_phase,
        seed_name=seed["seed_name"],
        composition=best_composition,
        tpd=best_tpd,
        iterations=best_iteration,
        converged=converged,
        unstable=best_tpd < -threshold,
        diagnostics={
            "tpd_threshold": threshold,
            "final_max_composition_delta": max_delta,
            "message": "negative TPD trial" if best_tpd < -threshold else "non-negative TPD trial",
        },
    )


def _tpd_value(composition: np.ndarray, trial_ln_phi: np.ndarray, feed: np.ndarray, parent_ln_phi: np.ndarray) -> float:
    return float(np.sum(composition * (np.log(composition) + trial_ln_phi - np.log(feed) - parent_ln_phi)))


def _composition_from_log_weights(log_weights: np.ndarray, min_composition: float) -> np.ndarray:
    shifted = np.asarray(log_weights, dtype=float)
    shifted = shifted - float(np.max(shifted))
    weights = np.exp(np.clip(shifted, -700.0, 700.0))
    return _clip_normalize_composition(weights, min_composition)


def _clip_normalize_composition(composition: np.ndarray, min_composition: float) -> np.ndarray:
    clipped = np.maximum(np.asarray(composition, dtype=float), float(min_composition))
    return clipped / float(np.sum(clipped))


def _rachford_rice_beta(feed: np.ndarray, k_values: np.ndarray) -> tuple[bool, float, str]:
    def rr(beta: float) -> float:
        return float(np.sum(feed * (k_values - 1.0) / (1.0 + beta * (k_values - 1.0))))

    f0 = rr(0.0)
    f1 = rr(1.0)
    if f0 <= 0.0:
        return False, 0.0, "no two-phase Rachford-Rice bracket; liquid-like single phase"
    if f1 >= 0.0:
        return False, 1.0, "no two-phase Rachford-Rice bracket; vapor-like single phase"
    lo = 0.0
    hi = 1.0
    for _ in range(100):
        mid = 0.5 * (lo + hi)
        if rr(mid) > 0.0:
            lo = mid
        else:
            hi = mid
    return True, 0.5 * (lo + hi), ""


def _phase_compositions(
    feed: np.ndarray, k_values: np.ndarray, beta: float, min_composition: float
) -> tuple[np.ndarray, np.ndarray]:
    x_liq = feed / (1.0 + beta * (k_values - 1.0))
    y_vap = k_values * x_liq
    x_liq = np.maximum(x_liq, min_composition)
    y_vap = np.maximum(y_vap, min_composition)
    x_liq = x_liq / np.sum(x_liq)
    y_vap = y_vap / np.sum(y_vap)
    return x_liq, y_vap


def _phase_from_state(
    label: str,
    composition: np.ndarray,
    phase_fraction: float,
    state_payload: dict[str, Any],
    options: EquilibriumOptions,
) -> EquilibriumPhase:
    return EquilibriumPhase(
        label=label,
        composition=composition,
        density=state_payload["density"],
        temperature=state_payload["state"].T,
        pressure=state_payload["state"].pressure(),
        phase_fraction=phase_fraction,
        ln_fugacity_coefficient=state_payload["ln_phi"],
        diagnostics=state_payload["diagnostics"] if options.include_phase_diagnostics else None,
    )


def _two_phase_result(
    best: dict[str, Any],
    options: EquilibriumOptions,
    message: str,
    stability_diagnostics: dict[str, Any],
) -> EquilibriumResult:
    beta = best["beta"]
    liquid_phase = _phase_from_state("liq", best["x_liq"], 1.0 - beta, best["liquid"], options)
    vapor_phase = _phase_from_state("vap", best["y_vap"], beta, best["vapor"], options)
    diagnostics = {
        "iterations": int(best["iteration"]),
        "fugacity_residual_norm": float(best["residual_norm"]),
        "fugacity_residual": best["fugacity_residual"].tolist(),
        "material_balance_error": float(best["material_error"]),
        "vapor_fraction": float(beta),
        "message": message,
        "point_solver_split_detected": True,
        "point_solver_message": message,
    }
    diagnostics.update(stability_diagnostics)
    return EquilibriumResult(
        backend="neutral_vle",
        problem_kind="tp_flash",
        phases=(liquid_phase, vapor_phase),
        stable=False,
        split_detected=True,
        diagnostics=diagnostics,
    )


def _initial_lle_guesses(
    initial_phases: Any,
    feed: np.ndarray,
    options: EquilibriumOptions,
    stability_precheck: dict[str, Any],
) -> list[dict[str, Any]]:
    if not isinstance(initial_phases, dict):
        if initial_phases is None:
            guesses = _default_lle_guesses(feed, options)
            tpd_seed = _tpd_lle_guess(feed, stability_precheck.get("unstable_trial_composition"), options)
            if tpd_seed is not None:
                return [tpd_seed, *guesses]
            return guesses
        raise InputError("initial_phases must be a dict with 'liq1', 'liq2', and 'phase_fraction'.")
    missing = {"liq1", "liq2", "phase_fraction"} - set(initial_phases)
    if missing:
        raise InputError("initial_phases is missing required key(s): {}.".format(", ".join(sorted(missing))))
    comp1 = _normalize_initial_phase(initial_phases["liq1"], feed.size, options.min_composition, "liq1")
    comp2 = _normalize_initial_phase(initial_phases["liq2"], feed.size, options.min_composition, "liq2")
    beta = float(initial_phases["phase_fraction"])
    if not np.isfinite(beta) or not (0.0 < beta < 1.0):
        raise InputError("initial_phases phase_fraction must be > 0 and < 1.")
    return [{"seed_name": "user", "beta": beta, "comp1": comp1, "comp2": comp2}]


def _tpd_lle_guess(
    feed: np.ndarray,
    trial_composition: Any,
    options: EquilibriumOptions,
) -> dict[str, Any] | None:
    if trial_composition is None:
        return None
    comp2 = _clip_normalize_composition(np.asarray(trial_composition, dtype=float), options.min_composition)
    if comp2.size != feed.size or _phase_distance(feed, comp2) <= _split_distance_tolerance(options):
        return None
    for beta in (0.5, 0.25, 0.75, 0.1, 0.9):
        comp1_raw = (feed - beta * comp2) / (1.0 - beta)
        if np.all(np.isfinite(comp1_raw)) and np.all(comp1_raw >= options.min_composition):
            comp1 = _clip_normalize_composition(comp1_raw, options.min_composition)
            if _phase_distance(comp1, comp2) > _split_distance_tolerance(options):
                return {"seed_name": "tpd_liq_trial", "beta": float(beta), "comp1": comp1, "comp2": comp2}
    lean_index = int(np.argmin(comp2 - feed))
    comp1 = _blend_lle_seed(
        feed,
        _component_rich_lle_composition(feed.size, lean_index, options.min_composition),
        0.9,
        options.min_composition,
    )
    if _phase_distance(comp1, comp2) <= _split_distance_tolerance(options):
        return None
    return {"seed_name": "tpd_liq_trial", "beta": 0.5, "comp1": comp1, "comp2": comp2}


def _default_lle_guesses(feed: np.ndarray, options: EquilibriumOptions) -> list[dict[str, Any]]:
    guesses: list[dict[str, Any]] = []
    ncomp = int(feed.size)
    if ncomp > 1:
        for strength in (0.9, 0.7):
            for component1 in range(ncomp):
                for component2 in range(ncomp):
                    if component1 == component2:
                        continue
                    comp1 = _blend_lle_seed(
                        feed,
                        _component_rich_lle_composition(ncomp, component2, options.min_composition),
                        strength,
                        options.min_composition,
                    )
                    comp2 = _blend_lle_seed(
                        feed,
                        _component_rich_lle_composition(ncomp, component1, options.min_composition),
                        strength,
                        options.min_composition,
                    )
                    guesses.append(
                        {
                            "seed_name": "auto_pair_{}_{}_s{}".format(component2, component1, int(100 * strength)),
                            "beta": 0.5,
                            "comp1": comp1,
                            "comp2": comp2,
                        }
                    )

    comp1, comp2 = _feed_perturb_lle_guess(feed, options)
    guesses.append({"seed_name": "auto_feed_perturb", "beta": 0.5, "comp1": comp1, "comp2": comp2})
    return guesses


def _component_rich_lle_composition(ncomp: int, rich_index: int, min_composition: float) -> np.ndarray:
    composition = np.full(int(ncomp), float(min_composition), dtype=float)
    composition[int(rich_index)] = max(float(min_composition), 1.0 - float(min_composition) * (int(ncomp) - 1))
    return composition / np.sum(composition)


def _blend_lle_seed(feed: np.ndarray, target: np.ndarray, strength: float, min_composition: float) -> np.ndarray:
    composition = (1.0 - float(strength)) * feed + float(strength) * target
    composition = np.maximum(composition, min_composition)
    return composition / np.sum(composition)


def _feed_perturb_lle_guess(feed: np.ndarray, options: EquilibriumOptions) -> tuple[np.ndarray, np.ndarray]:
    delta = min(0.2, 0.5 * float(np.min(feed)))
    comp1 = feed.copy()
    comp2 = feed.copy()
    if feed.size > 1 and delta > options.min_composition:
        comp1[0] = max(options.min_composition, comp1[0] - delta)
        comp1[1] = comp1[1] + delta
        comp2[0] = comp2[0] + delta
        comp2[1] = max(options.min_composition, comp2[1] - delta)
        comp1 = comp1 / np.sum(comp1)
        comp2 = comp2 / np.sum(comp2)
    return comp1, comp2


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


def _pack_lle_variables(beta: float, comp1: np.ndarray, comp2: np.ndarray) -> np.ndarray:
    beta = float(np.clip(beta, 1.0e-12, 1.0 - 1.0e-12))
    return np.concatenate(
        [
            np.asarray([np.log(beta / (1.0 - beta))]),
            _composition_to_logits(comp1),
            _composition_to_logits(comp2),
        ]
    )


def _unpack_lle_variables(variables: np.ndarray, ncomp: int) -> tuple[float, np.ndarray, np.ndarray]:
    raw = np.asarray(variables, dtype=float).flatten()
    if raw.size != 1 + 2 * (int(ncomp) - 1):
        raise SolutionError("Unexpected LLE variable vector size.")
    beta = float(1.0 / (1.0 + np.exp(-np.clip(raw[0], -700.0, 700.0))))
    offset = 1
    comp1 = _logits_to_composition(raw[offset : offset + ncomp - 1])
    offset += ncomp - 1
    comp2 = _logits_to_composition(raw[offset : offset + ncomp - 1])
    return beta, comp1, comp2


def _composition_to_logits(composition: np.ndarray) -> np.ndarray:
    comp = np.asarray(composition, dtype=float)
    return np.log(comp[:-1] / comp[-1])


def _logits_to_composition(logits: np.ndarray) -> np.ndarray:
    shifted = np.clip(np.asarray(logits, dtype=float), -700.0, 700.0)
    weights = np.concatenate([np.exp(shifted), np.asarray([1.0])])
    return weights / np.sum(weights)


def _evaluate_lle_variables(
    mixture: Any,
    T: float,
    P: float,
    feed: np.ndarray,
    variables: np.ndarray,
    options: EquilibriumOptions,
) -> dict[str, Any]:
    beta, comp1, comp2 = _unpack_lle_variables(variables, feed.size)
    state1 = _phase_state(mixture, T, P, comp1, "liq", options, "LLE flash")
    state2 = _phase_state(mixture, T, P, comp2, "liq", options, "LLE flash")
    fugacity_residual = np.log(comp2) + state2["ln_phi"] - np.log(comp1) - state1["ln_phi"]
    material_residual = (1.0 - beta) * comp1 + beta * comp2 - feed
    residual = np.concatenate([fugacity_residual, material_residual])
    return {
        "beta": beta,
        "comp1": comp1,
        "comp2": comp2,
        "state1": state1,
        "state2": state2,
        "fugacity_residual": fugacity_residual,
        "material_residual": material_residual,
        "residual": residual,
        "fugacity_residual_norm": float(np.max(np.abs(fugacity_residual))),
        "material_error": float(np.max(np.abs(material_residual))),
    }


def _solve_lle_attempt(
    mixture: Any,
    T: float,
    P: float,
    feed: np.ndarray,
    options: EquilibriumOptions,
    seed: dict[str, Any],
    attempt_count: int,
) -> dict[str, Any]:
    beta = float(seed["beta"])
    comp1 = seed["comp1"]
    comp2 = seed["comp2"]
    seed_name = str(seed["seed_name"])

    if _phase_distance(comp1, comp2) <= _split_distance_tolerance(options):
        return {
            "status": "degenerate",
            "seed_name": seed_name,
            "attempt_count": int(attempt_count),
            "message": "no V2 LLE split found; initial liquid phases are compositionally identical",
            "candidate": {
                "iteration": 0,
                "objective": 0.0,
                "fugacity_residual_norm": 0.0,
                "material_error": 0.0,
                "beta": beta,
                "comp1": comp1,
                "comp2": comp2,
            },
        }

    variables = _pack_lle_variables(beta, comp1, comp2)
    best: dict[str, Any] | None = None
    best_objective = np.inf
    for iteration in range(1, options.max_iterations + 1):
        current = _evaluate_lle_variables(mixture, T, P, feed, variables, options)
        objective = _lle_objective(current["residual"])
        if objective < best_objective:
            best = current | {"iteration": iteration, "objective": objective}
            best_objective = objective

        if _lle_converged(current, options):
            return {
                "status": "converged",
                "seed_name": seed_name,
                "attempt_count": int(attempt_count),
                "candidate": current | {"iteration": iteration, "objective": objective},
            }
        if _lle_degenerate(current, options):
            return {
                "status": "degenerate",
                "seed_name": seed_name,
                "attempt_count": int(attempt_count),
                "message": "no V2 LLE split found; phase split collapsed during iteration",
                "candidate": current | {"iteration": iteration, "objective": objective},
            }

        step = _lle_newton_step(
            lambda candidate: _evaluate_lle_variables(mixture, T, P, feed, candidate, options)["residual"],
            variables,
            current["residual"],
        )
        accepted: np.ndarray | None = None
        for scale in _damping_schedule(options.damping):
            candidate = variables + scale * step
            candidate_eval = _evaluate_lle_variables(mixture, T, P, feed, candidate, options)
            if _lle_objective(candidate_eval["residual"]) < objective:
                accepted = candidate
                break
        if accepted is None:
            assert best is not None
            return {
                "status": "failed",
                "seed_name": seed_name,
                "attempt_count": int(attempt_count),
                "message": "residual improvement stalled",
                "candidate": best,
            }
        variables = accepted

    assert best is not None
    if _lle_degenerate(best, options):
        return {
            "status": "degenerate",
            "seed_name": seed_name,
            "attempt_count": int(attempt_count),
            "message": "no V2 LLE split found; best candidate collapsed to one liquid phase",
            "candidate": best,
        }
    return {
        "status": "failed",
        "seed_name": seed_name,
        "attempt_count": int(attempt_count),
        "message": "maximum iterations reached",
        "candidate": best,
    }


def _lle_newton_step(residual_fn: Any, variables: np.ndarray, residual: np.ndarray) -> np.ndarray:
    jacobian = np.empty((residual.size, variables.size), dtype=float)
    for column in range(variables.size):
        step = 1.0e-5 * max(1.0, abs(float(variables[column])))
        forward = variables.copy()
        backward = variables.copy()
        forward[column] += step
        backward[column] -= step
        jacobian[:, column] = (residual_fn(forward) - residual_fn(backward)) / (2.0 * step)
    delta, *_ = np.linalg.lstsq(jacobian, -residual, rcond=None)
    return delta


def _damping_schedule(damping: float) -> tuple[float, ...]:
    start = float(np.clip(damping, 1.0e-6, 1.0))
    return tuple(start * factor for factor in (1.0, 0.5, 0.25, 0.1, 0.05, 0.01, 0.001))


def _lle_objective(residual: np.ndarray) -> float:
    return float(np.linalg.norm(np.asarray(residual, dtype=float), ord=2))


def _phase_distance(comp1: np.ndarray, comp2: np.ndarray) -> float:
    return float(np.max(np.abs(np.asarray(comp1, dtype=float) - np.asarray(comp2, dtype=float))))


def _split_distance_tolerance(options: EquilibriumOptions) -> float:
    return max(1.0e-8, 100.0 * options.min_composition)


def _lle_converged(candidate: dict[str, Any], options: EquilibriumOptions) -> bool:
    return (
        candidate["fugacity_residual_norm"] <= options.tolerance
        and candidate["material_error"] <= max(options.tolerance, 1.0e-10)
        and not _lle_degenerate(candidate, options)
    )


def _lle_degenerate(candidate: dict[str, Any], options: EquilibriumOptions) -> bool:
    beta = float(candidate["beta"])
    return (
        beta <= options.min_composition
        or beta >= 1.0 - options.min_composition
        or _phase_distance(candidate["comp1"], candidate["comp2"]) <= _split_distance_tolerance(options)
    )


def _lle_two_phase_result(
    best: dict[str, Any],
    options: EquilibriumOptions,
    message: str,
    seed_name: str,
    attempt_count: int,
    stability_diagnostics: dict[str, Any],
) -> EquilibriumResult:
    beta = best["beta"]
    phase1 = _phase_from_state("liq1", best["comp1"], 1.0 - beta, best["state1"], options)
    phase2 = _phase_from_state("liq2", best["comp2"], beta, best["state2"], options)
    diagnostics = {
        "iterations": int(best["iteration"]),
        "fugacity_residual_norm": float(best["fugacity_residual_norm"]),
        "fugacity_residual": best["fugacity_residual"].tolist(),
        "material_balance_error": float(best["material_error"]),
        "liquid2_phase_fraction": float(beta),
        "phase_distance": _phase_distance(best["comp1"], best["comp2"]),
        "seed_name": str(seed_name),
        "attempt_count": int(attempt_count),
        "message": message,
        "point_solver_split_detected": True,
        "point_solver_message": message,
    }
    diagnostics.update(stability_diagnostics)
    return EquilibriumResult(
        backend="neutral_lle",
        problem_kind="lle_flash",
        phases=(phase1, phase2),
        stable=False,
        split_detected=True,
        diagnostics=diagnostics,
    )


def _lle_no_split_result(
    feed: np.ndarray,
    state_feed: dict[str, Any],
    options: EquilibriumOptions,
    message: str,
    iterations: int,
    residual_norm: float,
    material_error: float,
    beta: float,
    comp1: np.ndarray,
    comp2: np.ndarray,
    seed_name: str,
    attempt_count: int,
    stability_diagnostics: dict[str, Any],
) -> EquilibriumResult:
    phase = _phase_from_state("liq", feed, 1.0, state_feed, options)
    diagnostics = {
        "iterations": int(iterations),
        "fugacity_residual_norm": float(residual_norm),
        "material_balance_error": float(material_error),
        "liquid2_phase_fraction": float(beta),
        "phase_distance": _phase_distance(comp1, comp2),
        "seed_name": str(seed_name),
        "attempt_count": int(attempt_count),
        "message": message,
        "point_solver_split_detected": False,
        "point_solver_message": message,
    }
    diagnostics.update(stability_diagnostics)
    return EquilibriumResult(
        backend="neutral_lle",
        problem_kind="lle_flash",
        phases=(phase,),
        stable=_no_split_stable_from_precheck(stability_diagnostics),
        split_detected=False,
        diagnostics=diagnostics,
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
    try:
        from .equilibrium_core.legacy_electrolyte_candidate import legacy_electrolyte_lle_candidate

        legacy = legacy_electrolyte_lle_candidate(mixture, T=T, P=P, feed=feed, options=options)
    except Exception as exc:
        legacy = {
            "legacy_candidate_found": False,
            "legacy_candidate_message": "legacy candidate fallback failed: {}".format(exc),
        }
    out.update(_json_like(legacy))
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
    native_initial_phases = initial_phases
    if initial_phases is not None:
        seed = _initial_lle_guesses(initial_phases, feed, opts, {"unstable_trial_composition": None})[0]
        native_initial_phases = {
            "liq1": seed["comp1"].tolist(),
            "liq2": seed["comp2"].tolist(),
            "phase_fraction": float(seed["beta"]),
        }
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
