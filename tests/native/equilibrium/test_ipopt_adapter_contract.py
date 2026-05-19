from __future__ import annotations

import epcsaft
import pytest
from epcsaft import _core


def test_native_ipopt_smoke_reports_generic_adapter_contract() -> None:
    smoke = _core._native_ipopt_smoke()

    assert smoke["backend"] == "ipopt"
    assert smoke["adapter_source_available"] is True
    assert smoke["adapter_kind"] == "native_tnlp_adapter"
    assert smoke["requires_exact_gradient"] is True
    assert smoke["requires_exact_jacobian"] is True
    assert smoke["available"] is smoke["compiled"]
    assert smoke["adapter_available"] is smoke["compiled"]


def test_runtime_capabilities_report_public_ipopt_routes() -> None:
    capabilities = epcsaft.capabilities()
    ipopt = capabilities["optimizers"]["ipopt"]

    assert ipopt["adapter_source_available"] is True
    assert ipopt["adapter_kind"] == "native_tnlp_adapter"
    assert ipopt["public_routes"] == [
        "reactive_speciation:ideal_mole_fraction",
        "reactive_speciation:mole_fraction_activity",
        "reactive_speciation:concentration",
        "neutral_tp_flash",
        "neutral_stability",
        "electrolyte_stability",
        "neutral_lle_flash",
        "neutral_bubble_p",
        "neutral_bubble_t",
        "neutral_dew_p",
        "neutral_dew_t",
        "electrolyte_lle",
        "electrolyte_bubble_pressure",
    ]


def test_native_ipopt_quadratic_smoke_is_gated_by_compiled_dependency() -> None:
    smoke = _core._native_ipopt_quadratic_smoke()

    assert smoke["backend"] == "ipopt"
    assert smoke["adapter_kind"] == "native_tnlp_adapter"
    assert smoke["problem"] == "quadratic_linear_constraint_smoke"
    if not smoke["compiled"]:
        assert smoke["ran"] is False
        assert smoke["accepted"] is False
        assert smoke["status"] == "ipopt_dependency_required"
    else:
        assert smoke["ran"] is True
        assert smoke["accepted"] is True
        assert smoke["exact_gradient_required"] is True
        assert smoke["exact_jacobian_required"] is True
        assert abs(smoke["variables"][0] - 1.0) < 1.0e-6
        assert abs(smoke["variables"][1] - 2.0) < 1.0e-6
        assert abs(smoke["constraints"][0] - 3.0) < 1.0e-6


def test_native_ipopt_quadratic_exact_hessian_reports_callback_and_bounded_history() -> None:
    smoke = _core._native_ipopt_quadratic_smoke(hessian_mode="exact", iteration_history_limit=2)

    assert smoke["backend"] == "ipopt"
    assert smoke["adapter_kind"] == "native_tnlp_adapter"
    if not smoke["compiled"]:
        assert smoke["status"] == "ipopt_dependency_required"
        return

    assert smoke["accepted"] is True
    assert smoke["hessian_approximation"] == "exact"
    assert smoke["hessian_backend"] == "analytic"
    assert smoke["exact_hessian_available"] is True
    assert smoke["eval_h_calls"] > 0
    assert smoke["iteration_count"] >= len(smoke["iteration_history"])
    assert len(smoke["iteration_history"]) <= 2
    assert smoke["iteration_history_limit"] == 2
    assert smoke["scaling_method"] == "user-scaling"
    assert smoke["variable_scaling_count"] == 2
    assert smoke["constraint_scaling_count"] == 1


def test_native_ipopt_quadratic_limited_memory_reports_explicit_mode() -> None:
    smoke = _core._native_ipopt_quadratic_smoke(hessian_mode="limited-memory")

    assert smoke["backend"] == "ipopt"
    if not smoke["compiled"]:
        assert smoke["status"] == "ipopt_dependency_required"
        return

    assert smoke["accepted"] is True
    assert smoke["hessian_approximation"] == "limited-memory"
    assert smoke["hessian_backend"] == "limited-memory"
    assert smoke["eval_h_calls"] == 0


def test_native_ipopt_quadratic_reports_linear_solver_and_tolerance_controls() -> None:
    smoke = _core._native_ipopt_quadratic_smoke(
        linear_solver="mumps",
        acceptable_tolerance=9.0e-7,
        constraint_violation_tolerance=8.0e-8,
        dual_infeasibility_tolerance=7.0e-8,
        complementarity_tolerance=6.0e-8,
    )

    assert smoke["backend"] == "ipopt"
    if not smoke["compiled"]:
        assert smoke["status"] == "ipopt_dependency_required"
        return

    assert smoke["accepted"] is True
    assert smoke["linear_solver_requested"] == "mumps"
    assert smoke["acceptable_tolerance"] == pytest.approx(9.0e-7)
    assert smoke["constraint_violation_tolerance"] == pytest.approx(8.0e-8)
    assert smoke["dual_infeasibility_tolerance"] == pytest.approx(7.0e-8)
    assert smoke["complementarity_tolerance"] == pytest.approx(6.0e-8)


def test_native_ipopt_quadratic_warm_start_round_trips_primal_and_multipliers() -> None:
    first = _core._native_ipopt_quadratic_smoke(hessian_mode="exact")
    if not first["compiled"]:
        assert first["status"] == "ipopt_dependency_required"
        return

    continuation_state = first["continuation_state"]
    second = _core._native_ipopt_quadratic_smoke(
        hessian_mode="exact",
        continuation_state=continuation_state,
    )

    assert second["accepted"] is True
    assert second["warm_start_requested"] is True
    assert second["warm_start_used"] is True
    assert second["continuation_state"]["variables"] == second["variables"]
    assert len(second["continuation_state"]["bound_lower_multipliers"]) == 2
    assert len(second["continuation_state"]["bound_upper_multipliers"]) == 2
    assert len(second["continuation_state"]["constraint_multipliers"]) == 1


def test_native_ipopt_quadratic_rejects_incompatible_continuation_state() -> None:
    if not _core._native_ipopt_smoke()["compiled"]:
        return

    bad_state = {
        "variables": [1.0],
        "bound_lower_multipliers": [0.0],
        "bound_upper_multipliers": [0.0],
        "constraint_multipliers": [0.0],
    }

    try:
        _core._native_ipopt_quadratic_smoke(continuation_state=bad_state)
    except Exception as exc:  # pybind exposes native ValueError as a module exception type.
        assert "continuation_state.variables" in str(exc)
    else:
        raise AssertionError("incompatible continuation_state was accepted")
