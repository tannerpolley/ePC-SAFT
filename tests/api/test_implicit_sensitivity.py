from __future__ import annotations

import json

import numpy as np
import pytest

import epcsaft


def test_implicit_sensitivity_from_jacobians_solves_linearized_state_response() -> None:
    result = epcsaft.implicit_sensitivity_from_jacobians(
        state=[0.25, 0.75],
        residual=[0.0, 0.0],
        residual_state_jacobian=[[2.0, 0.0], [0.0, 4.0]],
        residual_parameter_jacobian=[[6.0], [8.0]],
        backend="analytic_implicit",
        diagnostics={"block": "unit_test"},
    )

    assert isinstance(result, epcsaft.ImplicitSolveResult)
    assert result.backend == "analytic_implicit"
    assert result.status == "ok"
    np.testing.assert_allclose(result.sensitivity, [[-3.0], [-2.0]])
    payload = result.to_dict()
    assert payload["jacobians"]["residual_state"] == [[2.0, 0.0], [0.0, 4.0]]
    assert "finite_difference" not in json.dumps(payload).lower()


def test_implicit_sensitivity_rejects_nonimplicit_backend() -> None:
    with pytest.raises(epcsaft.InputError, match="analytic_implicit or cppad_implicit"):
        epcsaft.implicit_sensitivity_from_jacobians(
            state=[1.0],
            residual=[0.0],
            residual_state_jacobian=[[1.0]],
            residual_parameter_jacobian=[[1.0]],
            backend="analytic",
        )


def test_backend_unavailable_implicit_result_is_structured() -> None:
    result = epcsaft.backend_unavailable_implicit_result(
        state=[1.0],
        residual=[0.0],
        reason="missing residual parameter jacobian",
    )

    assert result.backend == "backend_unavailable"
    assert result.status == "backend_unavailable"
    assert result.sensitivity is None
    assert result.to_dict()["diagnostics"]["reason"] == "missing residual parameter jacobian"
