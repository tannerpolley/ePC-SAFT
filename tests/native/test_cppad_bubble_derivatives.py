from __future__ import annotations

import math

from epcsaft import _core


def test_cppad_pressure_residual_derivative_check_is_finite_when_available() -> None:
    payload = _core._native_autodiff_derivative_checks()

    assert payload["unsupported_derivative_used"] is False
    assert "pressure_log_residual" in payload["checked_residuals"]
    if payload["status"] == "unsupported_derivative":
        assert payload["derivative_backend"] == "cppad_unavailable"
        return

    derivative = payload["derivative_by_residual"]["pressure_log_residual"]
    assert payload["derivative_backend"] == "cppad"
    assert math.isfinite(derivative)
    assert derivative == 0.5


def test_native_implicit_sensitivity_reports_singular_bubble_linearization() -> None:
    payload = _core._solve_native_implicit_sensitivity(
        [1.0, 2.0, 2.0, 4.0],
        2,
        2,
        [1.0, 1.0],
        1,
    )

    assert payload["success"] is False
    assert payload["status"] == "singular_jacobian"
    assert payload["unsupported_derivative_used"] is False



