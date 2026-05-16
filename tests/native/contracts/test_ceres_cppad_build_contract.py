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

    jacobians = coverage["regression_ceres_jacobians"]
    assert jacobians["available"] is (ceres["available"] and cppad["available"])
    assert jacobians["production"] is (ceres["available"] and cppad["available"])
    assert jacobians["routes"] == ["pure_neutral_parameters", "binary_kij"]
    assert {row["quantity"] for row in coverage["rows"]}.issuperset({"pure_neutral_parameters", "binary_kij"})
    assert "numerical" + "_derivative" not in json.dumps({"ceres": ceres, "cppad": cppad, "coverage": coverage}).lower()
