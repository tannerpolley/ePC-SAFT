from __future__ import annotations

import math

import pytest

from epcsaft import _core


def test_cppad_reaction_residual_derivative_check_is_finite_when_available() -> None:
    payload = _core._native_autodiff_derivative_checks()

    assert payload["finite_difference_used"] is False
    assert "reaction_log_residual" in payload["checked_residuals"]
    if payload["status"] == "backend_unavailable":
        assert payload["derivative_backend"] == "cppad_unavailable"
        return

    derivative = payload["derivative_by_residual"]["reaction_log_residual"]
    assert payload["derivative_backend"] == "cppad"
    assert math.isfinite(derivative)
    assert derivative == 0.5


def test_native_implicit_sensitivity_solves_nested_speciation_linearization() -> None:
    payload = _core._solve_native_implicit_sensitivity(
        [2.0, 0.0, 0.0, 4.0],
        2,
        2,
        [-6.0, 8.0],
        1,
    )

    assert payload["success"] is True
    assert payload["status"] == "converged"
    assert payload["finite_difference_used"] is False
    assert payload["sensitivity_backend"] == "analytic_implicit"
    assert payload["shape"] == (2, 1)
    assert payload["sensitivities_row_major"] == [3.0, -2.0]


def test_native_implicit_sensitivity_accepts_overdetermined_full_rank_system() -> None:
    payload = _core._solve_native_implicit_sensitivity(
        [
            1.0,
            0.0,
            0.0,
            1.0,
            1.0,
            1.0,
        ],
        3,
        2,
        [-2.0, 4.0, 2.0],
        1,
    )

    assert payload["success"] is True
    assert payload["status"] == "converged"
    assert payload["finite_difference_used"] is False
    assert payload["sensitivity_backend"] == "analytic_implicit"
    assert payload["shape"] == (2, 1)
    assert payload["sensitivities_row_major"] == pytest.approx([2.0, -4.0])
