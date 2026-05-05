from __future__ import annotations

import json

import numpy as np
import pytest

from epcsaft.equilibrium_core.thermo_diagnostics import (
    compare_khudaida_digitized_paper_to_package,
    compare_khudaida_aad_tables,
    evaluate_khudaida_solver_gate,
    evaluate_khudaida_tieline,
    load_khudaida_digitized_paper_epcsaft,
    load_khudaida_tieline_case,
    summarize_khudaida_digitized_paper_matrix,
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

    assert diagnostics["source"] == "epcsaft_native_v5"
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
    assert "package_missing_count" in diagnostics
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
    assert "package_invalid_model_count" in diagnostics
    assert diagnostics["digitized_feed_figures"] == [2, 3, 4, 5, 6, 7]
    assert diagnostics["digitized_paper_epcsaft_figures"] == [2, 3, 4, 5, 6, 7]
    json.dumps(diagnostics, allow_nan=False)


def test_khudaida_model_tielines_cached_package_rows_are_finite() -> None:
    diagnostics = summarize_khudaida_matrix()

    assert diagnostics["case_count"] == 39
    assert diagnostics["package_invalid_model_count"] == 0
    assert diagnostics["package_cached_converged_count"] == 39
    assert np.isfinite(diagnostics["package_cached_residual_norm_max"])
    assert diagnostics["package_cached_strict_residual_pass_count"] >= 1
    assert diagnostics["package_cached_diagnostic_residual_envelope"] == 5.0e-2
    assert diagnostics["package_cached_diagnostic_residual_over_envelope_count"] == 0
    assert diagnostics["package_cached_residual_norm_max"] <= diagnostics["package_cached_diagnostic_residual_envelope"]
    assert diagnostics["package_cached_residual_norm_max_case"] == {
        "figure": 7,
        "tie_line": 6,
        "residual_norm": pytest.approx(0.04169193672341311),
    }
    json.dumps(diagnostics, allow_nan=False)


def test_khudaida_digitized_paper_epcsaft_series_loads_from_figure_folders() -> None:
    rows = load_khudaida_digitized_paper_epcsaft(figure=2)

    assert len(rows) == 8
    assert rows[0]["figure"] == 2
    assert rows[0]["point_id"] == 1
    assert rows[0]["tie_line"] == 1
    assert rows[0]["source"] == "digitized_user_supplied"
    np.testing.assert_allclose(np.sum(rows[0]["salt_free_composition"]), 1.0)
    assert set(rows[0]["salt_free_components"]) == {"H2O", "Ethanol", "Butanol"}
    json.dumps(rows, allow_nan=False)


def test_khudaida_digitized_paper_epcsaft_compares_to_package_model_rows() -> None:
    diagnostics = compare_khudaida_digitized_paper_to_package(figure=2)

    assert diagnostics["dataset"] == "2026_Khudaida"
    assert diagnostics["figure"] == 2
    assert diagnostics["rows_compared"] == 8
    assert diagnostics["finite_rows_compared"] == 8
    assert diagnostics["source"] == "digitized_user_supplied"
    assert diagnostics["package_source"] == "epcsaft_native_v5"
    assert np.isfinite(diagnostics["organic_salt_free_grand_aad"])
    assert np.isfinite(diagnostics["organic_salt_free_max_abs_error"])
    assert diagnostics["organic_salt_free_max_abs_error"] >= diagnostics["organic_salt_free_grand_aad"]
    assert diagnostics["rows"][0]["paper_epcsaft_tie_line"] == 1
    assert {"water_abs_error", "ethanol_abs_error", "isobutanol_abs_error"} <= set(diagnostics["rows"][0])
    json.dumps(diagnostics, allow_nan=False)


def test_khudaida_digitized_paper_epcsaft_matrix_summary_covers_all_lle_figures() -> None:
    diagnostics = summarize_khudaida_digitized_paper_matrix()

    assert diagnostics["dataset"] == "2026_Khudaida"
    assert diagnostics["figures"] == [2, 3, 4, 5, 6, 7]
    assert diagnostics["missing_data"] == []
    assert diagnostics["rows_compared"] == 39
    assert 0 < diagnostics["finite_rows_compared"] <= diagnostics["rows_compared"]
    assert (
        diagnostics["package_missing_or_invalid_rows"]
        == diagnostics["rows_compared"] - diagnostics["finite_rows_compared"]
    )
    assert np.isfinite(diagnostics["max_organic_salt_free_grand_aad"])
    assert np.isfinite(diagnostics["max_organic_salt_free_max_abs_error"])
    assert diagnostics["decision"] in {
        "package_matches_digitized_paper_epcsaft_on_salt_free_basis",
        "package_differs_from_digitized_paper_epcsaft_on_salt_free_basis",
    }
    json.dumps(diagnostics, allow_nan=False)


def test_khudaida_solver_gate_reports_algorithm_or_thermo_failure() -> None:
    diagnostics = evaluate_khudaida_solver_gate(figure=2, tie_line=1, source="package")

    assert diagnostics["source"] == "epcsaft_native_v5"
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


def test_khudaida_unseeded_solver_gate_accepts_known_split() -> None:
    diagnostics = evaluate_khudaida_solver_gate(figure=2, tie_line=1, source="package")

    assert diagnostics["solver_outcome"] == "accepted"
    assert diagnostics["acceptance_gate"] == "predictive_nonlinear_solve"
    assert diagnostics["decision"] == "solver_accepts_package_fixed_tieline_feed"


def test_khudaida_package_tieline_seeded_solver_gate_accepts_known_split() -> None:
    diagnostics = evaluate_khudaida_solver_gate(figure=2, tie_line=1, source="package", seeded=True)

    assert diagnostics["source"] == "epcsaft_native_v5"
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
