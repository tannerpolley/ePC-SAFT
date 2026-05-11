from __future__ import annotations

import math

from epcsaft import _core


def test_cppad_eos_derivative_check_is_not_finite_difference() -> None:
    payload = _core._native_autodiff_derivative_checks()

    assert payload["finite_difference_used"] is False
    assert payload["derivative_backend"] in {"cppad", "cppad_unavailable"}
    assert "scaled_residual" in payload["checked_residuals"]
    if payload["status"] == "backend_unavailable":
        assert payload["cppad_compiled"] is False
        return

    assert payload["status"] == "ok"
    assert payload["cppad_compiled"] is True
    assert payload["cppad_used"] is True
    assert math.isfinite(payload["max_abs_error"])
    assert payload["max_abs_error"] < 1.0e-12
