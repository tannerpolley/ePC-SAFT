"""Generic implicit sensitivity helpers for solved internal states."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field
from typing import Any

import numpy as np

from ._types import InputError

_IMPLICIT_BACKENDS = {"analytic_implicit", "cppad_implicit", "not_available"}
_EXPLICIT_TO_IMPLICIT_BACKEND = {
    "analytic": "analytic_implicit",
    "cppad": "cppad_implicit",
    "analytic_implicit": "analytic_implicit",
    "cppad_implicit": "cppad_implicit",
}


@dataclass(frozen=True, slots=True)
class ImplicitSolveResult:
    """Structured derivative payload for a solved state block."""

    state: Any
    residual: Any
    jacobians: Mapping[str, Any] = field(default_factory=dict)
    sensitivity: Any = None
    backend: str = "not_available"
    status: str = "not_available"
    diagnostics: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        backend = str(self.backend).strip().lower()
        if backend not in _IMPLICIT_BACKENDS:
            supported = "', '".join(sorted(_IMPLICIT_BACKENDS))
            raise InputError(f"ImplicitSolveResult.backend must be one of '{supported}'.")
        status = str(self.status).strip().lower()
        if "numerical_derivative" in backend or "numerical_derivative" in status:
            raise InputError("finite difference is not an allowed implicit sensitivity backend.")

        state = _array_payload(self.state, name="state")
        residual = _array_payload(self.residual, name="residual")
        sensitivity = None if self.sensitivity is None else _array_payload(self.sensitivity, name="sensitivity")
        jacobians = {str(key): _json_like(value) for key, value in dict(self.jacobians).items()}
        diagnostics = {str(key): _json_like(value) for key, value in dict(self.diagnostics).items()}

        object.__setattr__(self, "state", state)
        object.__setattr__(self, "residual", residual)
        object.__setattr__(self, "jacobians", jacobians)
        object.__setattr__(self, "sensitivity", sensitivity)
        object.__setattr__(self, "backend", backend)
        object.__setattr__(self, "status", status)
        object.__setattr__(self, "diagnostics", diagnostics)

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-like implicit solve payload."""
        return {
            "state": _json_like(self.state),
            "residual": _json_like(self.residual),
            "jacobians": _json_like(dict(self.jacobians)),
            "sensitivity": None if self.sensitivity is None else _json_like(self.sensitivity),
            "backend": self.backend,
            "status": self.status,
            "diagnostics": _json_like(dict(self.diagnostics)),
        }


def implicit_backend_for_residual_backend(backend: str) -> str:
    """Map an explicit residual-Jacobian backend onto its implicit solved-state backend."""

    return _EXPLICIT_TO_IMPLICIT_BACKEND.get(str(backend).strip().lower(), "not_available")


def implicit_sensitivity_from_jacobians(
    *,
    state: Any,
    residual: Any,
    residual_state_jacobian: Any,
    residual_parameter_jacobian: Any,
    backend: str = "analytic_implicit",
    diagnostics: Mapping[str, Any] | None = None,
) -> ImplicitSolveResult:
    """Compute solved-state sensitivities from residual Jacobians.

    For a converged residual system ``R(u, theta) = 0``, the local sensitivity is
    ``du/dtheta = -R_u^{-1} R_theta``.
    """

    normalized_backend = str(backend).strip().lower()
    if normalized_backend not in {"analytic_implicit", "cppad_implicit"}:
        raise InputError("implicit_sensitivity_from_jacobians requires analytic_implicit or cppad_implicit backend.")
    r_u = _matrix_payload(residual_state_jacobian, name="residual_state_jacobian")
    r_theta = _matrix_payload(residual_parameter_jacobian, name="residual_parameter_jacobian")
    if r_u.shape[0] != r_u.shape[1]:
        raise InputError("residual_state_jacobian must be square.")
    if r_theta.shape[0] != r_u.shape[0]:
        raise InputError("residual_parameter_jacobian row count must match residual_state_jacobian.")
    sensitivity = -np.linalg.solve(r_u, r_theta)
    return ImplicitSolveResult(
        state=state,
        residual=residual,
        jacobians={
            "residual_state": r_u,
            "residual_parameter": r_theta,
        },
        sensitivity=sensitivity,
        backend=normalized_backend,
        status="ok",
        diagnostics=dict(diagnostics or {}),
    )


def not_available_implicit_result(
    *,
    state: Any = (),
    residual: Any = (),
    reason: str,
    diagnostics: Mapping[str, Any] | None = None,
) -> ImplicitSolveResult:
    """Create a structured payload for an explicitly unavailable implicit derivative."""

    merged = dict(diagnostics or {})
    merged.setdefault("reason", str(reason))
    return ImplicitSolveResult(
        state=state,
        residual=residual,
        jacobians={},
        sensitivity=None,
        backend="not_available",
        status="not_available",
        diagnostics=merged,
    )


def _array_payload(value: Any, *, name: str) -> np.ndarray:
    array = np.asarray(value, dtype=float)
    if not np.all(np.isfinite(array)):
        raise InputError(f"{name} must contain only finite values.")
    return array


def _matrix_payload(value: Any, *, name: str) -> np.ndarray:
    matrix = _array_payload(value, name=name)
    if matrix.ndim != 2:
        raise InputError(f"{name} must be a two-dimensional matrix.")
    return matrix


def _json_like(value: Any) -> Any:
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, Mapping):
        return {str(key): _json_like(item) for key, item in value.items()}
    if isinstance(value, tuple):
        return [_json_like(item) for item in value]
    if isinstance(value, list):
        return [_json_like(item) for item in value]
    return value


__all__ = [
    "ImplicitSolveResult",
    "not_available_implicit_result",
    "implicit_backend_for_residual_backend",
    "implicit_sensitivity_from_jacobians",
]
