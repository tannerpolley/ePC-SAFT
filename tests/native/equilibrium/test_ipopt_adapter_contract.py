from __future__ import annotations

import epcsaft
from epcsaft import _core


def test_native_ipopt_smoke_reports_generic_adapter_contract() -> None:
    smoke = _core._native_ipopt_smoke()

    assert smoke["backend"] == "ipopt"
    assert smoke["adapter_source_available"] is True
    assert smoke["adapter_kind"] == "native_tnlp_adapter"
    removed_solver_detail = "hessian" + "_strategy"
    assert removed_solver_detail not in smoke
    assert smoke["requires_exact_gradient"] is True
    assert smoke["requires_exact_jacobian"] is True
    removed_hessian_requirement = "requires" + "_exact" + "_hessian"
    assert removed_hessian_requirement not in smoke
    assert smoke["available"] is smoke["compiled"]
    assert smoke["adapter_available"] is smoke["compiled"]


def test_runtime_capabilities_report_public_ipopt_routes() -> None:
    capabilities = epcsaft.capabilities()
    ipopt = capabilities["optimizers"]["ipopt"]

    assert ipopt["adapter_source_available"] is True
    assert ipopt["adapter_kind"] == "native_tnlp_adapter"
    assert ipopt["public_routes"] == [
        "reactive_speciation:ideal_mole_fraction",
        "neutral_tp_flash",
        "neutral_lle_flash",
        "neutral_bubble_p",
        "neutral_dew_p",
        "electrolyte_lle",
        "electrolyte_bubble_pressure",
    ]
    assert ipopt["full_constrained_nlp_available"] is ipopt["available"]


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
