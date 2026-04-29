"""Python-first phase-equilibrium helpers built on native ePC-SAFT states."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np

from ._types import InputError, SolutionError


@dataclass(frozen=True, slots=True)
class EquilibriumOptions:
    """Numerical controls for equilibrium solvers."""

    max_iterations: int = 80
    tolerance: float = 1.0e-6
    damping: float = 0.5
    min_composition: float = 1.0e-12
    include_phase_diagnostics: bool = False


@dataclass(frozen=True, slots=True)
class EquilibriumPhase:
    """One phase returned by an equilibrium calculation."""

    label: str
    composition: np.ndarray
    density: float
    temperature: float
    pressure: float
    phase_fraction: float
    fugacity_coefficient: np.ndarray | None = None
    diagnostics: dict[str, Any] | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "label", str(self.label))
        object.__setattr__(self, "composition", np.asarray(self.composition, dtype=float))
        object.__setattr__(self, "density", float(self.density))
        object.__setattr__(self, "temperature", float(self.temperature))
        object.__setattr__(self, "pressure", float(self.pressure))
        object.__setattr__(self, "phase_fraction", float(self.phase_fraction))
        if self.fugacity_coefficient is not None:
            object.__setattr__(self, "fugacity_coefficient", np.asarray(self.fugacity_coefficient, dtype=float))
        if self.diagnostics is not None:
            object.__setattr__(self, "diagnostics", dict(self.diagnostics))

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-like phase payload."""
        return {
            "label": self.label,
            "composition": self.composition.tolist(),
            "density": self.density,
            "temperature": self.temperature,
            "pressure": self.pressure,
            "phase_fraction": self.phase_fraction,
            "fugacity_coefficient": None
            if self.fugacity_coefficient is None
            else self.fugacity_coefficient.tolist(),
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


def tp_flash(mixture: Any, *, T: float, P: float, z: Any, options: EquilibriumOptions | None = None) -> EquilibriumResult:
    """Solve a V1 neutral TP flash with Rachford-Rice and fugacity updates."""
    opts = _normalize_options(options)
    feed = _normalize_feed(z, mixture.ncomp, opts.min_composition, "tp_flash")
    _reject_ion_containing_mixture(mixture)
    temperature = _positive_scalar(T, "T", "tp_flash")
    pressure = _positive_scalar(P, "P", "tp_flash")

    liquid_seed = _phase_state(mixture, temperature, pressure, feed, "liq", opts, "TP flash")
    vapor_seed = _phase_state(mixture, temperature, pressure, feed, "vap", opts, "TP flash")
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
        return EquilibriumResult(
            backend="neutral_vle",
            problem_kind="tp_flash",
            phases=(phase,),
            stable=True,
            split_detected=False,
            diagnostics={
                "iterations": 0,
                "fugacity_residual_norm": 0.0,
                "material_balance_error": 0.0,
                "vapor_fraction": float(beta),
                "message": no_split_message,
            },
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
            return _two_phase_result(best, opts, "converged")
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
    beta, comp1, comp2 = _initial_lle_guess(initial_phases, feed, opts)
    state_feed = _phase_state(mixture, temperature, pressure, feed, "liq", opts, "LLE flash")

    if _phase_distance(comp1, comp2) <= _split_distance_tolerance(opts):
        return _lle_no_split_result(
            feed,
            state_feed,
            opts,
            "no V2 LLE split found; initial liquid phases are compositionally identical",
            0,
            0.0,
            0.0,
            beta,
            comp1,
            comp2,
        )

    variables = _pack_lle_variables(beta, comp1, comp2)
    best: dict[str, Any] | None = None
    best_objective = np.inf
    stalled = False
    for iteration in range(1, opts.max_iterations + 1):
        current = _evaluate_lle_variables(mixture, temperature, pressure, feed, variables, opts)
        objective = _lle_objective(current["residual"])
        if objective < best_objective:
            best = current | {"iteration": iteration, "objective": objective}
            best_objective = objective

        if _lle_converged(current, opts):
            return _lle_two_phase_result(current | {"iteration": iteration}, opts, "converged")
        if _lle_degenerate(current, opts):
            return _lle_no_split_result(
                feed,
                state_feed,
                opts,
                "no V2 LLE split found; phase split collapsed during iteration",
                iteration,
                current["fugacity_residual_norm"],
                current["material_error"],
                current["beta"],
                current["comp1"],
                current["comp2"],
            )

        step = _lle_newton_step(
            lambda candidate: _evaluate_lle_variables(mixture, temperature, pressure, feed, candidate, opts)["residual"],
            variables,
            current["residual"],
        )
        accepted: np.ndarray | None = None
        for scale in _damping_schedule(opts.damping):
            candidate = variables + scale * step
            candidate_eval = _evaluate_lle_variables(mixture, temperature, pressure, feed, candidate, opts)
            if _lle_objective(candidate_eval["residual"]) < objective:
                accepted = candidate
                break
        if accepted is None:
            stalled = True
            break
        variables = accepted

    assert best is not None
    if stalled:
        return _lle_no_split_result(
            feed,
            state_feed,
            opts,
            "no V2 LLE split found; residual improvement stalled",
            int(best["iteration"]),
            best["fugacity_residual_norm"],
            best["material_error"],
            best["beta"],
            best["comp1"],
            best["comp2"],
        )
    if _lle_degenerate(best, opts):
        return _lle_no_split_result(
            feed,
            state_feed,
            opts,
            "no V2 LLE split found; best candidate collapsed to one liquid phase",
            int(best["iteration"]),
            best["fugacity_residual_norm"],
            best["material_error"],
            best["beta"],
            best["comp1"],
            best["comp2"],
        )
    raise SolutionError(
        "neutral LLE flash did not converge after {} iterations; residual_norm={}, material_balance_error={}, phase_distance={}".format(
            opts.max_iterations,
            best["fugacity_residual_norm"],
            best["material_error"],
            _phase_distance(best["comp1"], best["comp2"]),
        )
    )


def _normalize_options(options: EquilibriumOptions | None) -> EquilibriumOptions:
    if options is None:
        return EquilibriumOptions()
    if not isinstance(options, EquilibriumOptions):
        raise InputError("options must be an EquilibriumOptions instance.")
    if options.max_iterations <= 0:
        raise InputError("options.max_iterations must be positive.")
    if options.tolerance <= 0.0:
        raise InputError("options.tolerance must be positive.")
    if not (0.0 < options.damping <= 1.0):
        raise InputError("options.damping must be > 0 and <= 1.")
    if options.min_composition <= 0.0:
        raise InputError("options.min_composition must be positive.")
    return options


def _normalize_feed(z: Any, ncomp: int, min_composition: float, kind: str) -> np.ndarray:
    if z is None:
        raise InputError("z is required for kind='{}'.".format(kind))
    feed = np.asarray(z, dtype=float).flatten()
    if feed.size != int(ncomp):
        raise InputError("Feed composition length ({}) must match mixture component count ({}).".format(feed.size, ncomp))
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


def _phase_compositions(feed: np.ndarray, k_values: np.ndarray, beta: float, min_composition: float) -> tuple[np.ndarray, np.ndarray]:
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
        fugacity_coefficient=state_payload["ln_phi"],
        diagnostics=state_payload["diagnostics"] if options.include_phase_diagnostics else None,
    )


def _two_phase_result(best: dict[str, Any], options: EquilibriumOptions, message: str) -> EquilibriumResult:
    beta = best["beta"]
    liquid_phase = _phase_from_state("liq", best["x_liq"], 1.0 - beta, best["liquid"], options)
    vapor_phase = _phase_from_state("vap", best["y_vap"], beta, best["vapor"], options)
    return EquilibriumResult(
        backend="neutral_vle",
        problem_kind="tp_flash",
        phases=(liquid_phase, vapor_phase),
        stable=False,
        split_detected=True,
        diagnostics={
            "iterations": int(best["iteration"]),
            "fugacity_residual_norm": float(best["residual_norm"]),
            "fugacity_residual": best["fugacity_residual"].tolist(),
            "material_balance_error": float(best["material_error"]),
            "vapor_fraction": float(beta),
            "message": message,
        },
    )


def _initial_lle_guess(initial_phases: Any, feed: np.ndarray, options: EquilibriumOptions) -> tuple[float, np.ndarray, np.ndarray]:
    if initial_phases is None:
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
        return 0.5, comp1, comp2

    if not isinstance(initial_phases, dict):
        raise InputError("initial_phases must be a dict with 'liq1', 'liq2', and 'phase_fraction'.")
    missing = {"liq1", "liq2", "phase_fraction"} - set(initial_phases)
    if missing:
        raise InputError("initial_phases is missing required key(s): {}.".format(", ".join(sorted(missing))))
    comp1 = _normalize_initial_phase(initial_phases["liq1"], feed.size, options.min_composition, "liq1")
    comp2 = _normalize_initial_phase(initial_phases["liq2"], feed.size, options.min_composition, "liq2")
    beta = float(initial_phases["phase_fraction"])
    if not np.isfinite(beta) or not (0.0 < beta < 1.0):
        raise InputError("initial_phases phase_fraction must be > 0 and < 1.")
    return beta, comp1, comp2


def _normalize_initial_phase(value: Any, ncomp: int, min_composition: float, label: str) -> np.ndarray:
    composition = np.asarray(value, dtype=float).flatten()
    if composition.size != int(ncomp):
        raise InputError("initial_phases {} length ({}) must match mixture component count ({}).".format(label, composition.size, ncomp))
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


def _lle_two_phase_result(best: dict[str, Any], options: EquilibriumOptions, message: str) -> EquilibriumResult:
    beta = best["beta"]
    phase1 = _phase_from_state("liq1", best["comp1"], 1.0 - beta, best["state1"], options)
    phase2 = _phase_from_state("liq2", best["comp2"], beta, best["state2"], options)
    return EquilibriumResult(
        backend="neutral_lle",
        problem_kind="lle_flash",
        phases=(phase1, phase2),
        stable=False,
        split_detected=True,
        diagnostics={
            "iterations": int(best["iteration"]),
            "fugacity_residual_norm": float(best["fugacity_residual_norm"]),
            "fugacity_residual": best["fugacity_residual"].tolist(),
            "material_balance_error": float(best["material_error"]),
            "liquid2_phase_fraction": float(beta),
            "phase_distance": _phase_distance(best["comp1"], best["comp2"]),
            "message": message,
        },
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
) -> EquilibriumResult:
    phase = _phase_from_state("liq", feed, 1.0, state_feed, options)
    return EquilibriumResult(
        backend="neutral_lle",
        problem_kind="lle_flash",
        phases=(phase,),
        stable=True,
        split_detected=False,
        diagnostics={
            "iterations": int(iterations),
            "fugacity_residual_norm": float(residual_norm),
            "material_balance_error": float(material_error),
            "liquid2_phase_fraction": float(beta),
            "phase_distance": _phase_distance(comp1, comp2),
            "message": message,
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
