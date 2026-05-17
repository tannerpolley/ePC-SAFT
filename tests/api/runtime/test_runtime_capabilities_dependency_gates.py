from __future__ import annotations

import json

import epcsaft


def test_capabilities_report_cppad_without_legacy_forward_backend() -> None:
    derivatives = epcsaft.capabilities()["derivatives"]

    assert "numerical" + "_derivative" not in derivatives
    assert derivatives["cppad"]["scope"] == (
        "package-wide AD substrate; production derivative routes are listed in coverage_matrix"
    )
    assert "autodiff" not in derivatives
    assert "eigen_forward" not in derivatives


def test_derivative_coverage_matrix_has_required_contract_and_no_nonexact_derivative() -> None:
    coverage = epcsaft.capabilities()["derivatives"]["coverage_matrix"]

    assert coverage["derivative_coverage_matrix_available"] is True
    assert coverage["implemented_routes_only"] is True
    assert coverage["minimum_columns"] == [
        "row_family",
        "subsystem",
        "quantity",
        "derivative",
        "backend",
        "supported",
        "classification",
        "reason",
        "tests",
    ]
    rows = coverage["rows"]
    assert rows
    for row in rows:
        assert set(coverage["minimum_columns"]).issubset(row)
    assert "numerical" + "_derivative" not in json.dumps(coverage).lower()
    removed_label = "not" + "_" + "available"
    assert removed_label not in json.dumps(coverage).lower()


def test_derivative_coverage_matrix_reports_only_production_supported_routes() -> None:
    coverage = epcsaft.capabilities()["derivatives"]["coverage_matrix"]
    rows = coverage["rows"]

    classifications = {row["classification"] for row in rows}
    row_families = {row["row_family"] for row in rows}

    assert classifications == {"production_supported"}
    assert {"regression", "electrolyte_property"}.issubset(row_families)
    for row in rows:
        assert row["supported"] is True
        removed_column = "not" + "_applicable"
        assert removed_column not in row


def test_issue_68_required_coverage_gate_fields_are_reported_honestly() -> None:
    coverage = epcsaft.capabilities()["derivatives"]["coverage_matrix"]
    required = {
        "association_implicit_sensitivities",
        "density_root_implicit_sensitivities",
        "speciation_implicit_sensitivities",
        "born_ssmds_liquid_derivatives",
        "regression_ceres_jacobians",
    }
    assert required.issubset(coverage)
    assert coverage["association_implicit_sensitivities"]["production"] is True
    assert coverage["density_root_implicit_sensitivities"]["production"] is True
    assert coverage["born_ssmds_liquid_derivatives"]["phase_scope"] == "liquid_electrolyte_only"
    assert coverage["born_ssmds_liquid_derivatives"]["vapor_support"] is False
    assert coverage["regression_ceres_jacobians"]["routes"] == ["pure_neutral_parameters", "binary_kij"]


def test_reactive_batch_context_never_claims_ceres_native_hot_loop_in_default_build_contract() -> None:
    batch = epcsaft.capabilities()["regression"]["reactive_electrolyte_batch_context"]
    mixed = batch["mixed_pressure_speciation_residual_context"]

    assert mixed["production_optimizer"] is False
    assert mixed["optimizer"] is None
    assert mixed["native_hot_loop"] is False
    assert mixed["ceres"]["production"] is False
    assert mixed["python_role"] == "row orchestration and diagnostics"
    assert "numerical" + "_derivative" not in json.dumps(batch).lower()
