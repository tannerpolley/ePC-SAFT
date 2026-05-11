"""Native regression contract helpers.

This module is intentionally contract-only until the Ceres backend tranche is
implemented. It lets Python callers and tests inspect the native status and
fixed-shape schema that production native regression results must follow.
"""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any

from . import _core


def native_regression_contract_schema() -> dict[str, Any]:
    """Return the native regression status and fixed-shape residual contract."""

    return dict(_core._native_regression_contract_schema())


def evaluate_native_regression_residual_records(
    records: list[dict[str, Any]],
    *,
    penalty_residual: float = 1.0e6,
) -> dict[str, Any]:
    """Evaluate fixed-shape native residual records.

    This low-level helper validates native residual packing and row diagnostic
    semantics. It is not a production optimizer entrypoint.
    """

    return dict(_core._evaluate_native_regression_residual_records(records, float(penalty_residual)))


def solve_native_regression_residual_records(
    records: list[dict[str, Any]],
    parameters: list[dict[str, Any]],
    *,
    options: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Run the native fixed-shape residual contract solve surface.

    The current supported slice validates production derivative policy and
    status handling for already-evaluated residual records. It rejects
    finite-difference production derivatives instead of silently falling back.
    """

    return dict(_core._solve_native_regression_residual_records(records, parameters, dict(options or {})))


def evaluate_native_thermo_regression_rows(mixture: Any, request: Mapping[str, Any]) -> dict[str, Any]:
    """Evaluate supported thermodynamic regression rows in native C++.

    This is the issue #53 hot-loop boundary: Python supplies a validated native
    mixture and serialized row payload, while C++ runs the thermodynamic row
    solves and packs fixed-shape residuals.
    """

    native_mixture = getattr(mixture, "_native", mixture)
    return dict(_core._evaluate_native_thermo_regression_rows(native_mixture, dict(request)))


def fit_native_thermo_regression(mixture: Any, request: Mapping[str, Any]) -> dict[str, Any]:
    """Fit supported thermodynamic regression rows in the native hot loop.

    Python serializes the problem once. Native C++ applies parameters, evaluates
    thermodynamic rows, and owns the optimizer/derivative loop. Unsupported
    row/parameter derivative combinations return ``backend_unavailable`` rather
    than falling back to finite differences.
    """

    native_mixture = getattr(mixture, "_native", mixture)
    return dict(_core._fit_native_thermo_regression(native_mixture, dict(request)))


CANONICAL_NATIVE_REGRESSION_STATUSES: tuple[str, ...] = tuple(native_regression_contract_schema()["statuses"])


__all__ = [
    "CANONICAL_NATIVE_REGRESSION_STATUSES",
    "evaluate_native_regression_residual_records",
    "evaluate_native_thermo_regression_rows",
    "fit_native_thermo_regression",
    "native_regression_contract_schema",
    "solve_native_regression_residual_records",
]
