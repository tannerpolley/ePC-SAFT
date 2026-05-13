from __future__ import annotations

import json

import epcsaft


def test_disabled_default_build_does_not_claim_cppad_or_ceres_production() -> None:
    info = epcsaft.runtime_build_info()
    capabilities = epcsaft.capabilities()

    cppad = info["optional_dependencies"]["cppad"]
    ceres = info["optional_dependencies"]["ceres"]

    assert cppad["available"] is cppad["compiled"]
    assert ceres["available"] is ceres["compiled"]
    if not cppad["compiled"]:
        assert capabilities["derivatives"]["cppad"]["available"] is False
        assert capabilities["derivatives"]["cppad"]["production"] is False
        assert capabilities["derivatives"]["cppad"]["reason"] == "dependency_not_compiled"
    if not ceres["compiled"]:
        assert capabilities["optimizers"]["ceres"]["available"] is False
        assert capabilities["optimizers"]["ceres"]["production"] is False
        assert capabilities["optimizers"]["ceres"]["reason"] == "dependency_not_compiled"
    assert capabilities["optimizers"]["ceres"]["native_hot_loop"] is False


def test_capabilities_keep_eigen_forward_separate_from_cppad() -> None:
    derivatives = epcsaft.capabilities()["derivatives"]

    assert derivatives["finite_difference"] == {
        "available": False,
        "production": False,
        "reason": "finite_difference_derivatives_forbidden",
    }
    assert derivatives["eigen_forward"]["available"] is True
    assert derivatives["eigen_forward"]["scope"] == "legacy/local forward-mode AD"
    assert derivatives["cppad"]["scope"] == "package-wide AD substrate"
    assert "autodiff" not in derivatives


def test_derivative_coverage_matrix_has_required_contract_and_no_finite_difference() -> None:
    coverage = epcsaft.capabilities()["derivatives"]["coverage_matrix"]

    assert coverage["derivative_coverage_matrix_available"] is True
    assert coverage["minimum_columns"] == [
        "row_family",
        "subsystem",
        "quantity",
        "derivative",
        "backend",
        "supported",
        "not_applicable",
        "classification",
        "reason",
        "tests",
    ]
    rows = coverage["rows"]
    assert rows
    for row in rows:
        assert set(coverage["minimum_columns"]).issubset(row)
    assert "finite_difference" not in json.dumps(coverage).lower()


def test_derivative_coverage_matrix_classifies_runtime_rows_by_hard_gate_status() -> None:
    coverage = epcsaft.capabilities()["derivatives"]["coverage_matrix"]
    rows = coverage["rows"]

    classifications = {row["classification"] for row in rows}
    row_families = {row["row_family"] for row in rows}

    assert classifications == {"production_supported", "blocker", "out_of_scope"}
    assert {"regression", "solved_state", "electrolyte_property"}.issubset(row_families)
    for row in rows:
        if row["classification"] == "production_supported":
            assert row["supported"] is True
            assert row["not_applicable"] is False
            assert row["backend"] != "backend_unavailable"
        elif row["classification"] == "out_of_scope":
            assert row["supported"] is False
            assert row["not_applicable"] is True
        else:
            assert row["classification"] == "blocker"
            assert row["supported"] is False
            assert row["not_applicable"] is False


def test_issue_68_required_coverage_gate_fields_are_reported_honestly() -> None:
    coverage = epcsaft.capabilities()["derivatives"]["coverage_matrix"]
    required = {
        "association_direct_cppad_recording",
        "association_implicit_sensitivities",
        "density_root_implicit_sensitivities",
        "speciation_implicit_sensitivities",
        "bubble_pressure_implicit_sensitivities",
        "born_ssmds_liquid_derivatives",
        "regression_ceres_explicit_cppad_jacobians",
        "regression_ceres_implicit_jacobians",
    }
    assert required.issubset(coverage)
    assert coverage["association_direct_cppad_recording"] == {
        "available": False,
        "production": False,
        "reason": "active association uses solved site fractions; production derivative is implicit",
    }
    assert coverage["born_ssmds_liquid_derivatives"]["phase_scope"] == "liquid_electrolyte_only"
    assert coverage["born_ssmds_liquid_derivatives"]["vapor_support"] is False
    assert coverage["regression_ceres_implicit_jacobians"] == {
        "available": False,
        "production": False,
        "reason": "backend_unavailable",
    }


def test_reactive_batch_context_never_claims_ceres_native_hot_loop_in_default_contract() -> None:
    batch = epcsaft.capabilities()["regression"]["reactive_electrolyte_batch_context"]
    mixed = batch["bounded_mixed_pressure_speciation_regression"]

    assert mixed["native_hot_loop"] is False
    assert mixed["optimizer"] == "bounded_gauss_newton_least_squares"
    assert mixed["ceres"]["production"] is False
    assert "finite_difference" not in json.dumps(batch).lower()
