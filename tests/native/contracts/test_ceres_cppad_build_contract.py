from __future__ import annotations

import json

import epcsaft


def test_runtime_reports_ceres_build_contract() -> None:
    info = epcsaft.runtime_build_info()
    ceres = info["optional_dependencies"]["ceres"]

    assert ceres["backend"] == "ceres"
    assert ceres["status"] in {"disabled", "enabled_available", "not_configured"}
    assert ceres["compiled"] is (ceres["status"] == "enabled_available")
    assert ceres["available"] is ceres["compiled"]

    capabilities = epcsaft.capabilities()
    assert capabilities["optimizers"]["ceres"]["status"] == ceres["status"]
    assert capabilities["optimizers"]["ceres"]["compiled"] is ceres["compiled"]
    assert capabilities["optimizers"]["ceres"]["production"] is ceres["available"]
    assert capabilities["optimizers"]["ceres"]["native_hot_loop"] is ceres["available"]


def test_ceres_cppad_capability_claims_are_dependency_gated() -> None:
    capabilities = epcsaft.capabilities()
    ceres = capabilities["optimizers"]["ceres"]
    cppad = capabilities["derivatives"]["cppad"]
    coverage = capabilities["derivatives"]["coverage_matrix"]

    if not ceres["compiled"]:
        assert ceres["available"] is False
        assert ceres["reason"] == "dependency_not_compiled"
    if not cppad["compiled"]:
        assert cppad["available"] is False
        assert cppad["reason"] == "dependency_not_compiled"

    explicit = coverage["regression_ceres_explicit_cppad_jacobians"]
    assert explicit["available"] is (ceres["available"] and cppad["available"])
    assert explicit["production"] is False
    assert explicit["reason"] in {
        "dependency_not_compiled",
        "not_validated_for_production",
    }
    khb = next(row for row in coverage["rows"] if row["quantity"] == "binary_khb_ij")
    assert khb["classification"] == "blocker"
    assert khb["blocker_detail"] == "blocker_requires_implicit_association_sensitivity"
    assert khb["supported"] is False
    assert khb["future_owner"] == "Task C"
    assert khb["parameters"] == ["k_hb_ij"]
    assert "numerical_derivative" not in json.dumps({"ceres": ceres, "cppad": cppad, "coverage": coverage}).lower()
