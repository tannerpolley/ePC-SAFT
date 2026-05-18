"""Internal conversion helpers for native equilibrium payloads."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any

import numpy as np

from .._types import SolutionError


_ROUTE_STRING_DIAGNOSTIC_KEYS = (
    "solver_status",
    "application_status",
    "last_callback_exception",
    "problem_name",
    "adapter_kind",
    "gradient_approximation",
    "jacobian_approximation",
    "hessian_approximation",
    "hessian_backend",
    "scaling_method",
)

_ROUTE_BOOL_DIAGNOSTIC_KEYS = (
    "exact_gradient_required",
    "exact_jacobian_required",
    "exact_hessian_available",
    "warm_start_requested",
    "warm_start_used",
)

_ROUTE_INT_DIAGNOSTIC_KEYS = (
    "iteration_count",
    "iteration_history_limit",
    "iteration_history_size",
    "variable_scaling_count",
    "constraint_scaling_count",
    "eval_h_calls",
)


def _diagnostics(payload: Mapping[str, Any]) -> dict[str, Any]:
    diagnostics = payload.get("diagnostics", {})
    if isinstance(diagnostics, Mapping):
        return dict(diagnostics)
    return {}


def native_route_diagnostics(
    route: Mapping[str, Any],
    *,
    route_status_key: str = "route_status",
    defaults: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Return diagnostics for a native route acceptance gate."""
    postsolve = route.get("postsolve", {})
    diagnostics = dict(postsolve) if isinstance(postsolve, Mapping) else {}
    default_values = dict(defaults or {})
    diagnostics[route_status_key] = str(route.get("status", default_values.get("status", "")))
    diagnostics["solver_backend"] = str(route.get("backend", default_values.get("solver_backend", "")))
    for key in _ROUTE_STRING_DIAGNOSTIC_KEYS:
        if key in route or key in default_values:
            diagnostics[key] = str(route.get(key, default_values.get(key, "")))
    for key in _ROUTE_BOOL_DIAGNOSTIC_KEYS:
        if key in route or key in default_values:
            diagnostics[key] = bool(route.get(key, default_values.get(key, False)))
    for key in _ROUTE_INT_DIAGNOSTIC_KEYS:
        if key in route or key in default_values:
            diagnostics[key] = int(route.get(key, default_values.get(key, 0)))
    return diagnostics


def raise_native_route_rejected(
    route: Mapping[str, Any],
    message: str,
    *,
    route_status_key: str = "route_status",
    defaults: Mapping[str, Any] | None = None,
) -> None:
    """Raise a SolutionError with the shared native route diagnostics shape."""
    raise SolutionError(
        message,
        native_route_diagnostics(route, route_status_key=route_status_key, defaults=defaults),
    )


def _phase_payload_to_public(phase: Mapping[str, Any], *, label: str | None = None):
    from ..equilibrium import EquilibriumPhase

    ln_fugacity = phase.get("ln_fugacity_coefficient")
    fugacity = phase.get("fugacity_coefficient")
    return EquilibriumPhase(
        str(phase["label"] if label is None else label),
        composition=np.asarray(phase["composition"], dtype=float),
        density=float(phase["density"]),
        temperature=float(phase["temperature"]),
        pressure=float(phase["pressure"]),
        phase_fraction=float(phase["phase_fraction"]),
        ln_fugacity_coefficient=None if ln_fugacity is None else np.asarray(ln_fugacity, dtype=float),
        fugacity_coefficient=None if fugacity is None else np.asarray(fugacity, dtype=float),
        diagnostics=_diagnostics(phase),
    )


def neutral_two_phase_payload_to_result(
    payload: Mapping[str, Any],
    *,
    problem_kind: str | None = None,
    phase_labels: Sequence[str] | None = None,
):
    """Convert an accepted native neutral two-phase payload into public dataclasses."""
    from ..equilibrium import EquilibriumResult

    diagnostics = _diagnostics(payload)
    if not bool(payload.get("accepted", False)):
        reason = str(payload.get("rejection_reason", diagnostics.get("rejection_reason", "native_rejected")))
        if "rejection_reason" not in diagnostics:
            diagnostics["rejection_reason"] = reason
        raise SolutionError(f"Native neutral two-phase EOS result was rejected: {reason}", diagnostics)

    phases_raw = payload.get("phases", ())
    if not isinstance(phases_raw, Sequence) or isinstance(phases_raw, (str, bytes)):
        raise SolutionError("Native neutral two-phase EOS result did not contain a phase sequence.", diagnostics)
    if phase_labels is not None and len(phase_labels) != len(phases_raw):
        raise SolutionError("Native neutral two-phase EOS result label count did not match phase payloads.", diagnostics)
    labels = [None] * len(phases_raw) if phase_labels is None else list(phase_labels)
    phases = tuple(
        _phase_payload_to_public(phase, label=labels[index]) for index, phase in enumerate(phases_raw)
    )
    if not phases:
        raise SolutionError("Native neutral two-phase EOS result accepted without phase payloads.", diagnostics)

    resolved_problem_kind = payload.get("problem_kind", "neutral_two_phase_eos") if problem_kind is None else problem_kind
    return EquilibriumResult(
        backend=str(payload.get("backend", "native_equilibrium_nlp")),
        problem_kind=str(resolved_problem_kind),
        phases=phases,
        stable=bool(payload.get("stable", False)),
        split_detected=bool(payload.get("split_detected", True)),
        diagnostics=diagnostics,
    )


def native_route_summed_phase_amounts(route: Mapping[str, Any], ncomp: int, route_label: str) -> np.ndarray:
    """Return the positive feed implied by a native two-phase route result."""
    try:
        phase_amounts = np.asarray(route["phase_amounts"], dtype=float)
    except (KeyError, TypeError, ValueError) as exc:
        raise SolutionError(f"Native neutral {route_label} route did not return phase amounts.") from exc
    if phase_amounts.ndim != 2 or phase_amounts.shape[1] != int(ncomp):
        raise SolutionError(f"Native neutral {route_label} route phase amounts had an invalid shape.")
    feed = np.sum(phase_amounts, axis=0)
    if not np.all(np.isfinite(feed)) or np.any(feed <= 0.0):
        raise SolutionError(f"Native neutral {route_label} route phase amounts did not define a positive feed.")
    return feed


def native_route_solved_pressure(route: Mapping[str, Any], route_label: str) -> float:
    """Return the final pressure variable from a native fixed-temperature route."""
    try:
        variables = np.asarray(route["variables"], dtype=float).flatten()
    except (KeyError, TypeError, ValueError) as exc:
        raise SolutionError(f"Native neutral {route_label} route did not return solver variables.") from exc
    if variables.size == 0:
        raise SolutionError(f"Native neutral {route_label} route returned no solver variables.")
    pressure = float(variables[-1])
    if not np.isfinite(pressure) or pressure <= 0.0:
        raise SolutionError(f"Native neutral {route_label} route must be a finite positive P value.")
    return pressure


def native_route_solved_temperature(route: Mapping[str, Any], route_label: str) -> float:
    """Return the final temperature variable from a native fixed-pressure route."""
    try:
        variables = np.asarray(route["variables"], dtype=float).flatten()
    except (KeyError, TypeError, ValueError) as exc:
        raise SolutionError(f"Native neutral {route_label} route did not return solver variables.") from exc
    if variables.size == 0:
        raise SolutionError(f"Native neutral {route_label} route returned no solver variables.")
    temperature = float(variables[-1])
    if not np.isfinite(temperature) or temperature <= 0.0:
        raise SolutionError(f"Native neutral {route_label} route must be a finite positive T value.")
    return temperature
