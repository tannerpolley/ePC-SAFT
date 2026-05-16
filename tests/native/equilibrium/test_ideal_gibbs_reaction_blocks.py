from __future__ import annotations

import math

from epcsaft import _core


def test_ideal_reaction_smoke_satisfies_q_equals_k_and_stationarity() -> None:
    result = _core._native_ideal_reaction_smoke()

    assert result["model"] == "homogeneous_ideal_reaction"
    assert result["convex_kernel_scope"] == "homogeneous_ideal_reaction_validation"
    assert result["initial_amounts"] == [1.0, 0.0]
    assert result["extents"] == [0.75]
    assert result["amounts"] == [0.25, 0.75]
    assert result["mole_fractions"] == [0.25, 0.75]
    assert abs(result["log_q"][0] - math.log(3.0)) < 1.0e-14
    assert abs(result["residuals"][0]) < 1.0e-14
    assert abs(result["reaction_stationarity"]) < 1.0e-14
    assert abs(result["phase_validation_residuals"]["ideal_liquid"]) < 1.0e-14
    assert abs(result["phase_validation_residuals"]["ideal_vapor"]) < 1.0e-14


def test_ideal_gibbs_smoke_reports_exact_curvature_and_reaction_jacobian() -> None:
    result = _core._native_ideal_reaction_smoke()

    hessian = result["hessian_row_major"]
    expected_hessian = [3.0, -1.0, -1.0, 1.0 / 3.0]
    for actual, expected in zip(hessian, expected_hessian, strict=True):
        assert abs(actual - expected) < 1.0e-14
    jacobian = result["reaction_jacobian_row_major"]
    expected_jacobian = [-4.0, 4.0 / 3.0]
    for actual, expected in zip(jacobian, expected_jacobian, strict=True):
        assert abs(actual - expected) < 1.0e-14
