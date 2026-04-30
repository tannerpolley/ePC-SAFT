from __future__ import annotations

import json

import numpy as np
import pytest

from epcsaft.equilibrium_core.thermo_diagnostics import (
    compare_khudaida_aad_tables,
    evaluate_khudaida_solver_gate,
    evaluate_khudaida_tieline,
    load_khudaida_tieline_case,
    summarize_khudaida_matrix,
)


def test_khudaida_fixture_loads_charge_neutral_explicit_ions() -> None:
    case = load_khudaida_tieline_case(figure=2, tie_line=1)

    assert case["dataset"] == "2026_Khudaida"
    assert case["temperature_K"] == 293.15
    assert case["salt_wtfrac"] == 0.05
    assert case["species"] == ["H2O", "Ethanol", "Butanol", "Na+", "Cl-"]
    assert case["formula_species"] == ["H2O", "Ethanol", "Butanol", "NaCl"]

    for source in ("experimental", "package"):
        payload = case[source]
        assert payload["source"]
        assert abs(payload["feed_charge_balance"]) <= 1.0e-12
        assert abs(payload["organic_charge_balance"]) <= 1.0e-12
        assert abs(payload["aqueous_charge_balance"]) <= 1.0e-12
        np.testing.assert_allclose(np.sum(payload["feed_composition"]), 1.0)
        np.testing.assert_allclose(np.sum(payload["organic_composition"]), 1.0)
        np.testing.assert_allclose(np.sum(payload["aqueous_composition"]), 1.0)

    json.dumps(case, allow_nan=False)


def test_khudaida_package_tieline_fixed_phase_residual_is_internally_consistent() -> None:
    diagnostics = evaluate_khudaida_tieline(figure=2, tie_line=1, source="package")

    assert diagnostics["source"] == "epcsaft_package"
    assert diagnostics["material_balance_error"] <= 1.0e-10
    assert diagnostics["charge_balance_error"] <= 1.0e-8
    assert diagnostics["cached_model_residual_norm"] <= 1.0e-6
    assert diagnostics["decision"] == "package_fixed_tieline_internally_consistent"
    assert diagnostics["fixed_phase_residual_norm"] <= 1.0e-6
    assert diagnostics["gibbs_delta"] < 0.0
    assert "NaCl" in diagnostics["mean_ionic_residuals"]
    assert {"H2O", "Ethanol", "Butanol"} <= set(diagnostics["neutral_residuals"])
    assert {"feed", "organic", "aqueous"} <= set(diagnostics["fugacity_contribution_terms"])
    json.dumps(diagnostics, allow_nan=False)


def test_khudaida_experimental_tieline_reports_thermodynamic_surface_diagnostics() -> None:
    diagnostics = evaluate_khudaida_tieline(figure=2, tie_line=1, source="experimental")

    assert diagnostics["source"] == "paper_table"
    assert diagnostics["decision"] in {
        "experimental_tieline_matches_current_surface",
        "thermodynamic_surface_differs_from_reference_tieline",
        "no_fixed_tieline_satisfies_current_surface",
    }
    assert diagnostics["material_balance_error"] <= 1.0e-10
    assert diagnostics["charge_balance_error"] <= 1.0e-8
    assert np.isfinite(diagnostics["fixed_phase_residual_norm"])
    assert "density_branch" in diagnostics
    assert "mean_ionic_residuals" in diagnostics
    assert "gibbs_delta" in diagnostics
    json.dumps(diagnostics, allow_nan=False)


def test_khudaida_aad_matrix_records_package_surface_gap_to_paper() -> None:
    diagnostics = compare_khudaida_aad_tables()

    assert diagnostics["dataset"] == "2026_Khudaida"
    assert diagnostics["decision"] == "package_aad_exceeds_paper_epcsaft_reference"
    assert diagnostics["tables"] == ["table_9", "table_10"]
    assert diagnostics["max_package_grand_aad"] > diagnostics["max_paper_epcsaft_grand_aad"]
    assert diagnostics["rows_compared"] == 6
    json.dumps(diagnostics, allow_nan=False)


def test_khudaida_full_matrix_fixture_is_complete_and_charge_neutral() -> None:
    diagnostics = summarize_khudaida_matrix()

    assert diagnostics["dataset"] == "2026_Khudaida"
    assert diagnostics["figures"] == [2, 3, 4, 5, 6, 7]
    assert diagnostics["temperatures_K"] == [293.15, 303.15, 313.15]
    assert diagnostics["salt_wtfracs"] == [0.05, 0.1]
    assert diagnostics["case_count"] > 0
    assert diagnostics["max_charge_balance_error"] <= 1.0e-12
    assert diagnostics["missing_data"] == []
    json.dumps(diagnostics, allow_nan=False)


def test_khudaida_solver_gate_reports_algorithm_or_thermo_failure() -> None:
    diagnostics = evaluate_khudaida_solver_gate(figure=2, tie_line=1, source="package")

    assert diagnostics["source"] == "epcsaft_package"
    assert np.isfinite(diagnostics["fixed_phase_residual_norm"])
    assert diagnostics["solver_outcome"] in {"accepted", "rejected"}
    assert diagnostics["acceptance_gate"] in {
        "predictive_nonlinear_solve",
        "predictive_solve_failed",
    }
    assert diagnostics["decision"] in {
        "fixed_tieline_consistent_solver_suspect",
        "solver_accepts_package_fixed_tieline_feed",
        "current_surface_inconsistent_before_solver",
    }
    assert diagnostics["fixed_phase_residual_norm"] <= 1.0e-6
    assert diagnostics["decision"] in {
        "fixed_tieline_consistent_solver_suspect",
        "solver_accepts_package_fixed_tieline_feed",
    }
    assert "ascani_case2_fixture_regression" not in json.dumps(diagnostics)
    assert "v4_partition_seed_api_compatibility" not in json.dumps(diagnostics)
    json.dumps(diagnostics, allow_nan=False)


@pytest.mark.xfail(reason="Unseeded Khudaida feed recovery remains a tracked v4 robustness follow-up.")
def test_khudaida_unseeded_solver_gate_is_not_yet_an_acceptance_test() -> None:
    diagnostics = evaluate_khudaida_solver_gate(figure=2, tie_line=1, source="package")

    assert diagnostics["solver_outcome"] == "accepted"
    assert diagnostics["acceptance_gate"] == "predictive_nonlinear_solve"
    assert diagnostics["decision"] == "solver_accepts_package_fixed_tieline_feed"


def test_khudaida_package_tieline_seeded_solver_gate_accepts_known_split() -> None:
    diagnostics = evaluate_khudaida_solver_gate(figure=2, tie_line=1, source="package", seeded=True)

    assert diagnostics["source"] == "epcsaft_package"
    assert diagnostics["fixed_phase_residual_norm"] <= 1.0e-6
    assert diagnostics["gibbs_delta"] < 0.0
    assert diagnostics["solver_outcome"] == "accepted"
    assert diagnostics["acceptance_gate"] == "predictive_nonlinear_solve"
    assert diagnostics["decision"] == "solver_accepts_package_fixed_tieline_feed"
    assert diagnostics["solver_diagnostics"]["solver_seed_name"] == "initial_phases"
    assert diagnostics["solver_diagnostics"]["solver_residual_norm"] <= 1.0e-6
    assert diagnostics["solver_diagnostics"]["gibbs_delta"] < 0.0
    assert diagnostics["solver_diagnostics"]["phase_labels_swapped"] is False
    json.dumps(diagnostics, allow_nan=False)
