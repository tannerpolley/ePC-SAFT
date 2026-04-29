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
    feed = _normalize_feed(z, mixture.ncomp, opts.min_composition)
    _reject_ion_containing_mixture(mixture)
    temperature = _positive_scalar(T, "T")
    pressure = _positive_scalar(P, "P")

    liquid_seed = _phase_state(mixture, temperature, pressure, feed, "liq", opts)
    vapor_seed = _phase_state(mixture, temperature, pressure, feed, "vap", opts)
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
        liquid = _phase_state(mixture, temperature, pressure, x_liq, "liq", opts)
        vapor = _phase_state(mixture, temperature, pressure, y_vap, "vap", opts)
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


def _normalize_feed(z: Any, ncomp: int, min_composition: float) -> np.ndarray:
    if z is None:
        raise InputError("z is required for kind='tp_flash'.")
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
        raise InputError("V1 TP flash requires each feed composition entry to be >= min_composition.")
    return feed


def _positive_scalar(value: Any, label: str) -> float:
    if value is None:
        raise InputError("{} is required for kind='tp_flash'.".format(label))
    out = float(value)
    if not np.isfinite(out) or out <= 0.0:
        raise InputError("{} must be a positive finite scalar.".format(label))
    return out


def _reject_ion_containing_mixture(mixture: Any) -> None:
    z = np.asarray(mixture.parameters.get("z", []), dtype=float).flatten()
    if z.size and np.any(np.abs(z) > 1.0e-12):
        raise InputError("V1 equilibrium does not support ion-containing mixtures.")


def _phase_state(mixture: Any, T: float, P: float, composition: np.ndarray, label: str, options: EquilibriumOptions) -> dict[str, Any]:
    try:
        state = mixture.state(T=T, P=P, x=composition, phase=label)
    except SolutionError:
        raise
    except Exception as exc:
        raise SolutionError("Failed to construct {} phase during TP flash: {}".format(label, exc)) from exc
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
