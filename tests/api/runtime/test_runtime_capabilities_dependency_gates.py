from __future__ import annotations

import json

import epcsaft
import epcsaft.capability_evidence as capability_evidence


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


def test_derivative_coverage_matrix_reports_only_production_supported_routes() -> None:
    coverage = epcsaft.capabilities()["derivatives"]["coverage_matrix"]
    rows = coverage["rows"]

    classifications = {row["classification"] for row in rows}
    row_families = {row["row_family"] for row in rows}

    assert classifications == {"production_supported"}
    assert {"regression", "electrolyte_property"}.issubset(row_families)
    for row in rows:
        assert row["supported"] is True


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


def test_capability_contract_is_derived_from_registered_evidence() -> None:
    capabilities = epcsaft.capabilities()
    evidence = capabilities["capability_evidence"]

    assert evidence["source"] == "registered_capability_evidence"
    assert evidence["ipopt_public_routes"] == capabilities["optimizers"]["ipopt"]["public_routes"]
    assert set(evidence["equilibrium_keys"]).issubset(capabilities["equilibrium"])
    assert evidence["equilibrium_route_derivative_row_count"] == len(
        capabilities["derivatives"]["equilibrium_route_evidence"]["rows"]
    )
    assert evidence["problem_object_classes"] == capabilities["equilibrium"]["problem_objects"]["classes"]
    assert set(evidence["regression_keys"]).issubset(capabilities["regression"])
    assert evidence["derivative_row_count"] == len(capabilities["derivatives"]["coverage_matrix"]["rows"])
    assert evidence["validation_lanes"] == list(capability_evidence.VALIDATION_LANES)
    assert evidence["pytest_slices"] == list(capability_evidence.TEST_SLICES)
    assert evidence["cheap_validation_lanes"] == [
        name for name, lane in capability_evidence.VALIDATION_LANES.items() if lane["cheap_by_default"]
    ]


def test_registered_evidence_links_capability_rows_to_executable_checks() -> None:
    for row in (*capability_evidence.DERIVATIVE_COVERAGE_ROWS, *capability_evidence.EQUILIBRIUM_ROUTE_DERIVATIVE_EVIDENCE):
        assert row["tests"]
        assert all(str(target).startswith("tests/") for target in row["tests"])

    for lane in capability_evidence.VALIDATION_LANES.values():
        commands = lane["commands"]
        assert commands
        assert all(command for command in commands)

    for slice_payload in capability_evidence.TEST_SLICES.values():
        targets = slice_payload["targets"]
        assert targets
        assert all(str(target).startswith("tests") for target in targets)


def test_route_level_equilibrium_derivative_rows_are_registered_without_capability_overclaims() -> None:
    route_evidence = epcsaft.capabilities()["derivatives"]["equilibrium_route_evidence"]
    by_quantity = {row["quantity"]: row for row in route_evidence["rows"]}

    assert by_quantity["neutral_two_phase_routes"]["classification"] == "production_supported"
    assert by_quantity["electrolyte_lle_and_stability"]["classification"] == "production_supported"
    assert (
        by_quantity["reactive_lle_and_reactive_electrolyte_lle"]["classification"]
        == "route_builder_supported_capability_pending"
    )
    assert "reactive_lle" not in epcsaft.capabilities()["equilibrium"]
    assert "reactive_electrolyte_lle" not in epcsaft.capabilities()["equilibrium"]
