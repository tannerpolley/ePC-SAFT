"""Regression API contract tests outside the hydrocarbon benchmark."""

from __future__ import annotations

import csv
from types import SimpleNamespace

import numpy as np
import pytest

import epcsaft
import epcsaft.regression as regression_module
from epcsaft import FitProblem, FitResult, create_parameter_template, fit_pure_neutral, write_fit_result
from epcsaft._types import InputError
from epcsaft.regression import (
    _debug_native_pure_neutral_objective,
    _fit_pure_neutral_least_squares_internal,
    evaluate_generic_regression_derivatives,
)
from tests.helpers.regression_cases import _methane_like_records, _minimal_neutral_metadata


def test_public_regression_surface_includes_ion_and_binary_v1():
    assert hasattr(epcsaft, "fit_pure_parameters")
    assert hasattr(epcsaft, "fit_binary_parameters")
    assert hasattr(epcsaft, "fit_liquid_electrolyte_parameters")
    assert hasattr(epcsaft, "fit_pure_neutral")
    assert hasattr(epcsaft, "fit_pure_ion")
    assert hasattr(epcsaft, "fit_binary_pair")
    assert hasattr(epcsaft, "FitParameter")
    assert hasattr(epcsaft, "BinaryInteraction")
    assert hasattr(epcsaft, "RelativePermittivityResidual")
    assert hasattr(epcsaft, "validate_regression_provenance")
    assert hasattr(epcsaft, "evaluate_reactive_electrolyte_bubble_residuals")
    assert hasattr(epcsaft, "ReactiveElectrolyteRegressionResult")
    assert hasattr(epcsaft, "TargetRow")
    assert hasattr(epcsaft, "TargetDataset")


def test_reactive_electrolyte_regression_residuals_keep_fixed_shape(monkeypatch):
    calls = []

    def fake_solve(**kwargs):
        calls.append(kwargs)
        if len(calls) == 2:
            raise epcsaft.SolutionError("synthetic bubble failure", {"best_P": 90000.0})
        return SimpleNamespace(
            success=True,
            message="converged",
            P_total=101325.0,
            y_vap={"CO2": 0.25, "H2O": 0.75},
            partial_pressures={"CO2": 25331.25},
            x_liq={"CO2": 0.1, "H2O": 0.9},
            named_reaction_residuals={"R": 0.2},
            fugacity_residual_norm=1.0e-8,
            state_failure_count=0,
        )

    monkeypatch.setattr("epcsaft.reactive_electrolyte.solve_reactive_electrolyte_bubble", fake_solve)
    result = epcsaft.evaluate_reactive_electrolyte_bubble_residuals(
        [
            {
                "row_id": "ok",
                "T": 298.15,
                "P_seed": 101325.0,
                "totals": {"carbon": 0.1, "water": 0.9},
                "initial_x": [0.1, 0.9],
                "target_partial_pressures": {"CO2": 25331.25},
                "target_x": {"CO2": 0.1},
            },
            {
                "row_id": "bad",
                "T": 298.15,
                "P_seed": 101325.0,
                "totals": {"carbon": 0.1, "water": 0.9},
                "initial_x": [0.1, 0.9],
                "target_partial_pressures": {"CO2": 25331.25},
                "target_x": {"CO2": 0.1},
            },
        ],
        species=["CO2", "H2O"],
        mixture_factory=lambda x, T, P: None,
        balances={"carbon": {"CO2": 1.0}, "water": {"H2O": 1.0}},
        reactions=[],
        vapor_species=["CO2", "H2O"],
        pressure_species=["CO2"],
        speciation_species=["CO2"],
        reaction_names=["R"],
        pressure_weight=4.0,
        speciation_weight=0.25,
        reaction_weight=9.0,
        penalty_value=8.0,
    )

    assert result.success_count == 1
    assert result.failure_count == 1
    assert result.residuals.shape == (6,)
    assert result.residual_names == (
        "ok.partial_pressure.CO2",
        "ok.x.CO2",
        "ok.reaction.R",
        "bad.partial_pressure.CO2",
        "bad.x.CO2",
        "bad.reaction.R",
    )
    assert result.residuals[0] == pytest.approx(0.0)
    assert result.residuals[1] == pytest.approx(0.0)
    assert result.residuals[2] == pytest.approx(0.6)
    assert result.residuals[3:].tolist() == pytest.approx([16.0, 4.0, 24.0])
    assert result.record_results[0]["partial_pressures"] == {"CO2": pytest.approx(25331.25)}
    assert result.record_results[0]["x_liq"] == {"CO2": pytest.approx(0.1), "H2O": pytest.approx(0.9)}
    assert result.record_results[0]["y_vap"] == {"CO2": pytest.approx(0.25), "H2O": pytest.approx(0.75)}
    assert result.record_results[0]["named_reaction_residuals"] == {"R": pytest.approx(0.2)}
    assert result.to_dict()["diagnostics"]["fixed_shape"] is True
    assert len(calls) == 2
    assert calls[0]["options"].error_mode == "result"


def test_reactive_electrolyte_regression_residuals_use_continuation_seed(monkeypatch):
    calls = []

    def fake_solve(**kwargs):
        calls.append(kwargs)
        return SimpleNamespace(
            success=True,
            message="converged",
            P_total=120000.0 if len(calls) == 1 else 121000.0,
            y_vap={"CO2": 0.2, "H2O": 0.8},
            partial_pressures={"CO2": 24000.0},
            x_liq={"CO2": 0.1, "H2O": 0.9},
            named_reaction_residuals={},
            fugacity_residual_norm=1.0e-8,
            state_failure_count=0,
        )

    monkeypatch.setattr("epcsaft.reactive_electrolyte.solve_reactive_electrolyte_bubble", fake_solve)
    options = epcsaft.ReactiveElectrolyteBubbleOptions(
        bubble_options=epcsaft.ElectrolyteBubbleOptions(initial_pressure=101325.0),
        error_mode="result",
    )
    epcsaft.evaluate_reactive_electrolyte_bubble_residuals(
        [
            {
                "row_id": "one",
                "T": 298.15,
                "P_seed": 101325.0,
                "totals": {"carbon": 0.1, "water": 0.9},
                "initial_x": [0.1, 0.9],
                "target_partial_pressures": {"CO2": 24000.0},
            },
            {
                "row_id": "two",
                "T": 298.15,
                "P_seed": 90000.0,
                "totals": {"carbon": 0.1, "water": 0.9},
                "initial_x": [0.1, 0.9],
                "target_partial_pressures": {"CO2": 24000.0},
            },
        ],
        species=["CO2", "H2O"],
        mixture_factory=lambda x, T, P: None,
        balances={"carbon": {"CO2": 1.0}, "water": {"H2O": 1.0}},
        reactions=[],
        vapor_species=["CO2", "H2O"],
        pressure_species=["CO2"],
        options=options,
        continuation="auto",
    )

    assert calls[0]["P_seed"] == pytest.approx(101325.0)
    assert calls[1]["P_seed"] == pytest.approx(120000.0)
    assert calls[1]["options"].bubble_options.initial_pressure == pytest.approx(120000.0)
    assert calls[1]["options"].bubble_options.initial_y_vap == pytest.approx({"CO2": 0.2, "H2O": 0.8})


def test_provenance_validation_rejects_indirect_dborn_without_electrostatic_data():
    with pytest.raises(InputError, match=r"d_born.*dielectric|d_born.*ion-activity"):
        epcsaft.validate_regression_provenance(
            [
                epcsaft.FitParameter(
                    "MEAH+",
                    "d_born",
                    source="mixed_reactive_vle",
                )
            ],
            strict=True,
        )

    report = epcsaft.validate_regression_provenance(
        [
            epcsaft.FitParameter(
                "MEAH+",
                "d_born",
                source="explicit_override",
                allow_without_direct_data=True,
            )
        ],
        strict=True,
    )
    assert report["warnings"] == []
    assert report["parameter_sources"]["MEAH+.d_born"] == "explicit_override"


def test_provenance_validation_gates_ion_involving_binary_interactions():
    with pytest.raises(InputError, match="same-sign ionic pair"):
        epcsaft.validate_regression_provenance(
            [
                epcsaft.BinaryInteraction(
                    ("MEACOO-", "HCO3-"),
                    parameter="k_ij",
                    source="direct_binary_vle",
                )
            ],
            species=["MEACOO-", "HCO3-"],
            charges=[-1.0, -1.0],
            strict=True,
        )

    with pytest.raises(InputError, match=r"neutral-ion.*direct"):
        epcsaft.validate_regression_provenance(
            [
                epcsaft.BinaryInteraction(
                    ("CO2", "MEACOO-"),
                    parameter="k_ij",
                    source="mixed_reactive_vle",
                )
            ],
            species=["CO2", "MEACOO-"],
            charges=[0.0, -1.0],
            strict=True,
        )

    report = epcsaft.validate_regression_provenance(
        [
            epcsaft.BinaryInteraction(
                ("MEAH+", "MEACOO-"),
                parameter="k_ij",
                source="direct_electrolyte_activity",
            )
        ],
        species=["MEAH+", "MEACOO-"],
        charges=[1.0, -1.0],
        strict=True,
    )
    assert report["warnings"] == []
    assert report["data_sources_by_parameter"]["MEAH+:MEACOO-.k_ij"] == ["direct_electrolyte_activity"]


def test_relative_permittivity_residual_builds_first_class_fit_term_record():
    residual = epcsaft.RelativePermittivityResidual(
        T=298.15,
        P=101325.0,
        composition={"H2O": 0.8, "MEA": 0.2},
        epsilon_r_exp=65.0,
        weight=2.0,
        source="dielectric_measurement",
    )

    term = residual.to_fit_term(species=["H2O", "MEA"])

    assert term.term_type == "relative_permittivity"
    assert term.weight == pytest.approx(2.0)
    assert term.residual_count == 1
    assert term.records[0]["epsilon_r_exp"] == pytest.approx(65.0)
    assert term.records[0]["x_H2O"] == pytest.approx(0.8)
    assert term.records[0]["x_MEA"] == pytest.approx(0.2)


def test_fit_pure_neutral_requires_pressure_for_density_records():
    with pytest.raises(InputError, match="experimental 'P'"):
        fit_pure_neutral(
            [{"T": 320.0, "rho": 9000.0}],
            "Toluene",
            assoc_scheme="",
            fixed_parameters=_minimal_neutral_metadata(92.141e-3),
            initial_guess={"m": 2.8, "s": 3.7, "e": 285.0},
        )


def test_fit_pure_neutral_rejects_non_phase1_targets():
    with pytest.raises(InputError, match="supports only the targets 'm', 's', and 'e'"):
        fit_pure_neutral(
            [{"T": 320.0, "P": 101325.0, "rho": 9000.0}],
            "Toluene",
            assoc_scheme="",
            fit_targets=("m", "s", "e_assoc"),
            fixed_parameters=_minimal_neutral_metadata(92.141e-3),
            initial_guess={"m": 2.8, "s": 3.7, "e_assoc": 1000.0},
        )


def test_native_pure_neutral_debug_gradient_reports_autodiff_backend():
    theta = {"m": 1.05, "s": 3.68, "e": 151.0}
    debug = _debug_native_pure_neutral_objective(
        _methane_like_records(),
        "Methane",
        assoc_scheme="",
        fixed_parameters=_minimal_neutral_metadata(16.043e-3),
        initial_guess=theta,
        x=theta,
    )
    exact = np.asarray(debug["gradient"], dtype=float)
    assert np.all(np.isfinite(exact))
    assert debug["jacobian_backend"] == "autodiff"
    assert debug["residual_evaluations"] >= 1
    assert debug["density_solves"] >= 2
    assert debug["fused_state_evaluations"] >= 2
    assert debug["callback_wall_time_s"] >= 0.0
    assert debug["jacobian_available"] is True
    assert debug["jacobian_backend"] == "autodiff"
    assert debug["jacobian_fallback_used"] is False
    assert tuple(debug["jacobian_shape"]) == (len(debug["residuals"]), 3)
    assert np.asarray(debug["jacobian_row_major"], dtype=float).shape == (len(debug["residuals"]) * 3,)
    assert debug["hessian_available"] is False
    assert debug["hessian_backend"] == "not_implemented"
    assert tuple(debug["hessian_shape"]) == (0, 0)
    assert np.asarray(debug["hessian_row_major"], dtype=float).size == 0


def test_internal_native_least_squares_backend_matches_methane_reference_band():
    result = _fit_pure_neutral_least_squares_internal(
        _methane_like_records(),
        "Methane",
        assoc_scheme="",
        fixed_parameters=_minimal_neutral_metadata(16.043e-3),
        initial_guess={"m": 1.08, "s": 3.55, "e": 155.0},
        bounds={
            "m": (0.5, 3.5),
            "s": (2.0, 5.0),
            "e": (50.0, 400.0),
        },
    )

    assert result.success, result.message
    assert result.backend == "least_squares_native"
    assert result.metrics_by_term["density"] < 0.02
    assert result.metrics_by_term["pure_vle_fugacity_balance"] < 0.02
    assert result.fitted_values["m"] == pytest.approx(1.0, rel=0.0, abs=0.05)
    assert result.fitted_values["s"] == pytest.approx(3.7039, rel=0.0, abs=0.08)
    assert result.fitted_values["e"] == pytest.approx(150.03, rel=0.0, abs=3.0)


@pytest.mark.parametrize(
    "initial_guess",
    [
        {"m": 1.08, "s": 3.55, "e": 155.0},
        {"m": 0.92, "s": 3.90, "e": 143.0},
        {"m": 1.20, "s": 3.30, "e": 170.0},
    ],
)
def test_public_pure_neutral_regression_is_robust_to_distinct_initial_guesses(initial_guess):
    result = fit_pure_neutral(
        _methane_like_records(),
        "Methane",
        assoc_scheme="",
        fixed_parameters=_minimal_neutral_metadata(16.043e-3),
        initial_guess=initial_guess,
        bounds={
            "m": (0.5, 3.5),
            "s": (2.0, 5.0),
            "e": (50.0, 400.0),
        },
    )
    assert result.success, result.message
    assert result.metrics_by_term["density"] < 0.02
    assert result.metrics_by_term["pure_vle_fugacity_balance"] < 0.02
    assert result.fitted_values["m"] == pytest.approx(1.0, rel=0.0, abs=0.06)
    assert result.fitted_values["s"] == pytest.approx(3.7039, rel=0.0, abs=0.10)
    assert result.fitted_values["e"] == pytest.approx(150.03, rel=0.0, abs=4.0)


def _minimal_nacl_records():
    return [
        {
            "T": 298.15,
            "P": 101325.0,
            "x_H2O": 0.996,
            "x_Na+": 0.002,
            "x_Cl-": 0.002,
            "osmotic_coefficient": 0.974,
            "mean_ionic_activity": 0.922,
        }
    ]


def _stub_native_generic_runner(monkeypatch):
    calls = []

    def fake_runner(
        fixed_payloads,
        native_records,
        optimization_names,
        species,
        theta0,
        lower,
        upper,
        *,
        component=None,
        pair=None,
        multistart=0,
        max_nfev=200,
    ):
        calls.append(
            {
                "fixed_payloads": fixed_payloads,
                "native_records": native_records,
                "optimization_names": tuple(optimization_names),
                "species": tuple(species),
                "theta0": np.asarray(theta0, dtype=float),
                "lower": np.asarray(lower, dtype=float),
                "upper": np.asarray(upper, dtype=float),
                "component": component,
                "pair": pair,
                "multistart": int(multistart),
                "max_nfev": int(max_nfev),
            }
        )
        metrics = {str(record["term_name"]): 0.0 for record in native_records}
        if not metrics:
            metrics = {"residual": 0.0}
        return {
            "x": np.asarray(theta0, dtype=float),
            "cost": 0.0,
            "residual_norm": 0.0,
            "initial_cost": 0.0,
            "initial_residual_norm": 0.0,
            "metrics_by_term": metrics,
            "success": True,
            "status": 1,
            "nfev": 1,
            "iterations": 0,
            "starts_tried": 1,
            "message": "stubbed native generic regression",
            "backend": "least_squares_native",
            "jacobian_available": True,
            "jacobian_backend": "stub",
            "jacobian_fallback_used": False,
            "jacobian_fallback_reason": "",
            "backend_unavailable_reason": "",
            "hessian_available": False,
            "hessian_backend": "not_implemented",
            "hessian_fallback_used": False,
            "hessian_fallback_reason": "stubbed hessian skeleton",
        }

    monkeypatch.setattr(regression_module, "_run_native_generic_least_squares", fake_runner)
    return calls


def test_fit_pure_ion_requires_composition_and_activity_or_osmotic_records():
    with pytest.raises(InputError, match="composition"):
        epcsaft.fit_pure_ion(
            [{"T": 298.15, "P": 101325.0, "osmotic_coefficient": 0.93}],
            "Na+",
            dataset="2026_Khudaida",
        )

    with pytest.raises(InputError, match=r"osmotic|mean-ionic|mean ionic"):
        epcsaft.fit_pure_ion(
            [{"T": 298.15, "P": 101325.0, "molality": 0.1}],
            "Na+",
            dataset="2026_Khudaida",
            species=["H2O", "Na+", "Cl-"],
            solvent="H2O",
        )


def test_fit_pure_ion_default_s_e_bounds_and_multistart_contract(monkeypatch):
    calls = _stub_native_generic_runner(monkeypatch)
    result = epcsaft.fit_pure_ion(
        _minimal_nacl_records(),
        "Na+",
        dataset="2026_Khudaida",
        species=["H2O", "Na+", "Cl-"],
        solvent="H2O",
        initial_guess={"s": 2.6, "e": 210.0},
        bounds={"s": (2.4, 3.2), "e": (150.0, 300.0)},
        multistart=3,
    )

    assert result.success, result.message
    assert result.backend == "least_squares_native"
    assert result.jacobian_available is True
    assert result.jacobian_backend == "stub"
    assert result.backend_unavailable_reason == ""
    assert result.hessian_available is False
    assert result.hessian_backend == "not_implemented"
    assert result.problem.mode == "pure_ion"
    assert result.problem.fit_targets == ("s", "e")
    assert result.metrics_by_term == {"osmotic_coefficient": 0.0, "mean_ionic_activity": 0.0}
    assert result.fitted_values == {"s": 2.6, "e": 210.0}
    assert len(calls) == 1
    call = calls[0]
    assert call["optimization_names"] == ("s", "e")
    assert call["component"] == "Na+"
    assert call["species"] == ("H2O", "Na+", "Cl-")
    assert call["multistart"] == 3
    np.testing.assert_allclose(call["theta0"], [2.6, 210.0])
    np.testing.assert_allclose(call["lower"], [2.4, 150.0])
    np.testing.assert_allclose(call["upper"], [3.2, 300.0])
    assert {record["term_name"] for record in call["native_records"]} == {"osmotic_coefficient", "mean_ionic_activity"}


def test_fit_pure_ion_accepts_d_born_and_born_user_options(monkeypatch):
    calls = _stub_native_generic_runner(monkeypatch)
    user_options = {
        "elec_model": {
            "rel_perm": {"rule": "empirical", "differential_mode": "autodiff"},
            "born_model": {
                "d_Born_mode": 3,
                "solvation_shell_model": True,
                "dielectric_saturation": True,
                "mu_born_model": {"differential_mode": "autodiff", "comp_dep_delta_d": True},
            },
        }
    }
    result = epcsaft.fit_pure_ion(
        _minimal_nacl_records(),
        "Na+",
        dataset="2026_Khudaida",
        species=["H2O", "Na+", "Cl-"],
        solvent="H2O",
        fit_targets=("d_born",),
        initial_guess={"d_born": 3.2},
        bounds={"d_born": (2.0, 5.0)},
        user_options=user_options,
    )

    assert result.success, result.message
    assert result.backend == "least_squares_native"
    assert result.problem.fit_targets == ("d_born",)
    assert result.metrics_by_term["osmotic_coefficient"] == 0.0
    assert len(calls) == 1
    assert calls[0]["optimization_names"] == ("d_born",)
    assert calls[0]["component"] == "Na+"
    assert calls[0]["fixed_payloads"][0]["elec_model"]["rel_perm"]["differential_mode"] == 2
    assert calls[0]["fixed_payloads"][0]["elec_model"]["born_model"]["mu_born_model"]["differential_mode"] == 2
    assert result.provenance_report["parameter_sources"]["Na+.d_born"] == "ion_activity"
    assert result.provenance_report["warnings"] == []


def test_fit_pure_ion_passes_explicit_mean_ionic_pair_label_to_native_backend(monkeypatch):
    calls = _stub_native_generic_runner(monkeypatch)
    records = [dict(record, pair_label="Na+Cl-") for record in _minimal_nacl_records()]
    result = epcsaft.fit_pure_ion(
        records,
        "Na+",
        dataset="2026_Khudaida",
        species=["H2O", "Na+", "Cl-"],
        solvent="H2O",
        initial_guess={"s": 2.6, "e": 210.0},
        bounds={"s": (2.4, 3.2), "e": (150.0, 300.0)},
    )

    assert result.success, result.message
    assert result.backend == "least_squares_native"
    assert result.metrics_by_term["mean_ionic_activity"] == 0.0
    miac_records = [record for record in calls[0]["native_records"] if record["term_name"] == "mean_ionic_activity"]
    assert len(miac_records) == 1
    assert miac_records[0]["target_index"] == 1
    assert miac_records[0]["target_index_2"] == 2


def test_fit_binary_pair_vle_kij_default_and_rejects_temperature_models(monkeypatch):
    calls = _stub_native_generic_runner(monkeypatch)
    records = [
        {"T": 330.0, "P": 101325.0, "x_H2O": 0.7, "x_Ethanol": 0.3, "y_H2O": 0.5, "y_Ethanol": 0.5},
        {"T": 340.0, "P": 101325.0, "x_H2O": 0.6, "x_Ethanol": 0.4, "y_H2O": 0.4, "y_Ethanol": 0.6},
    ]
    result = epcsaft.fit_binary_pair(
        records,
        ("H2O", "Ethanol"),
        dataset="2026_Khudaida",
        initial_guess={"k_ij": -0.02},
        bounds={"k_ij": (-0.2, 0.2)},
        multistart=2,
    )

    assert result.success, result.message
    assert result.backend == "least_squares_native"
    assert result.jacobian_available is True
    assert result.jacobian_backend == "stub"
    assert result.hessian_available is False
    assert result.problem.mode == "binary_pair"
    assert result.problem.fit_targets == ("k_ij",)
    assert set(result.fitted_values) == {"k_ij"}
    assert result.metrics_by_term == {"binary_vle_fugacity_balance": 0.0}
    assert len(calls) == 1
    assert calls[0]["optimization_names"] == ("k_ij",)
    assert calls[0]["pair"] == ("H2O", "Ethanol")
    assert calls[0]["multistart"] == 2

    with pytest.raises(InputError, match="temperature_model"):
        epcsaft.fit_binary_pair(
            records,
            ("H2O", "Ethanol"),
            dataset="2026_Khudaida",
            temperature_model="linear",
        )


def test_fit_binary_pair_can_fit_all_constant_binary_interaction_targets(monkeypatch):
    calls = _stub_native_generic_runner(monkeypatch)
    records = [
        {"T": 330.0, "P": 101325.0, "x_H2O": 0.7, "x_Ethanol": 0.3, "y_H2O": 0.5, "y_Ethanol": 0.5},
        {"T": 340.0, "P": 101325.0, "x_H2O": 0.6, "x_Ethanol": 0.4, "y_H2O": 0.4, "y_Ethanol": 0.6},
    ]

    result = epcsaft.fit_binary_pair(
        records,
        ("H2O", "Ethanol"),
        dataset="2026_Khudaida",
        fit_targets=("k_ij", "l_ij", "k_hb_ij"),
        initial_guess={"k_ij": -0.02, "l_ij": 0.01, "k_hb_ij": 0.02},
        bounds={"k_ij": (-0.2, 0.2), "l_ij": (-0.2, 0.2), "k_hb_ij": (-0.2, 0.2)},
        multistart=1,
    )

    assert result.success, result.message
    assert result.backend == "least_squares_native"
    assert result.problem.fit_targets == ("k_ij", "l_ij", "k_hb_ij")
    assert set(result.fitted_values) == {"k_ij", "l_ij", "k_hb_ij"}
    assert all(np.isfinite(value) for value in result.fitted_values.values())
    assert result.metrics_by_term == {"binary_vle_fugacity_balance": 0.0}
    assert len(calls) == 1
    assert calls[0]["optimization_names"] == ("k_ij", "l_ij", "k_hb_ij")
    assert calls[0]["pair"] == ("H2O", "Ethanol")
    assert calls[0]["multistart"] == 1
    assert result.provenance_report["parameter_sources"]["H2O:Ethanol.k_ij"] == "direct_binary_vle"


def test_fit_binary_pair_rejects_nonpositive_vle_fractions_before_native_solve():
    records = [
        {"T": 330.0, "P": 101325.0, "x_H2O": 0.7, "x_Ethanol": 0.3, "y_H2O": 0.0, "y_Ethanol": 1.0},
    ]
    with pytest.raises(InputError, match="strictly positive"):
        epcsaft.fit_binary_pair(
            records,
            ("H2O", "Ethanol"),
            dataset="2026_Khudaida",
            initial_guess={"k_ij": -0.02},
        )


def test_fit_binary_pair_rejects_ion_involving_kij_without_direct_electrolyte_provenance():
    records = [
        {
            "T": 298.15,
            "P": 101325.0,
            "x_H2O": 0.998,
            "x_Na+": 0.001,
            "x_Cl-": 0.001,
            "y_H2O": 0.998,
            "y_Na+": 0.001,
            "y_Cl-": 0.001,
        }
    ]

    with pytest.raises(InputError, match=r"opposite-sign ionic pair.*direct electrolyte"):
        epcsaft.fit_binary_pair(
            records,
            ("Na+", "Cl-"),
            dataset="2026_Khudaida",
            species=["H2O", "Na+", "Cl-"],
            initial_guess={"k_ij": 0.0},
        )


def test_write_fit_result_updates_only_target_pure_row(tmp_path):
    dataset_root = create_parameter_template(tmp_path, "pure_case", ["H2O", "Na+"])
    result = FitResult(
        problem=FitProblem(
            mode="pure_neutral",
            component="H2O",
            fit_targets=("m", "s", "e"),
            optimization_parameters=("m", "s", "e"),
        ),
        fitted_values={"m": 1.2047, "s": 2.7927, "e": 353.95},
        rendered_values={"m": 1.2047, "s": 2.7927, "e": 353.95},
        success=True,
    )

    written = write_fit_result(result, dataset_root, overwrite=False)
    assert written == [dataset_root / "pure" / "water.csv"]

    with (dataset_root / "pure" / "water.csv").open("r", encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle))

    target_row = next(row for row in rows if row["component"] == "H2O")
    other_row = next(row for row in rows if row["component"] == "Na+")
    assert target_row["m"] == "1.2047"
    assert target_row["s"] == "2.7927"
    assert target_row["e"] == "353.95"
    assert other_row["m"] == ""
    assert other_row["s"] == ""
    assert other_row["e"] == ""


def test_write_fit_result_updates_ion_row_and_binary_matrix_symmetrically(tmp_path):
    ion_root = create_parameter_template(tmp_path, "ion_case", ["H2O", "Na+"])
    ion_result = FitResult(
        problem=FitProblem(
            mode="pure_ion",
            component="Na+",
            fit_targets=("s", "e"),
            optimization_parameters=("s", "e"),
        ),
        fitted_values={"s": 2.84, "e": 231.2},
        rendered_values={"s": 2.84, "e": 231.2},
        success=True,
    )

    written = write_fit_result(ion_result, ion_root, overwrite=False)
    assert written == [ion_root / "pure" / "water.csv"]
    with (ion_root / "pure" / "water.csv").open("r", encoding="utf-8-sig", newline="") as handle:
        ion_rows = list(csv.DictReader(handle))
    ion_row = next(row for row in ion_rows if row["component"] == "Na+")
    assert ion_row["s"] == "2.84"
    assert ion_row["e"] == "231.2"

    binary_root = create_parameter_template(tmp_path, "binary_case", ["H2O", "Ethanol"])
    binary_result = FitResult(
        problem=FitProblem(
            mode="binary_pair",
            pair=("H2O", "Ethanol"),
            fit_targets=("k_ij", "l_ij", "k_hb_ij"),
            optimization_parameters=("k_ij", "l_ij", "k_hb_ij"),
        ),
        fitted_values={"k_ij": -0.06167, "l_ij": 0.0123, "k_hb_ij": -0.0345},
        rendered_values={"k_ij": -0.06167, "l_ij": 0.0123, "k_hb_ij": -0.0345},
        success=True,
    )

    written = write_fit_result(binary_result, binary_root, overwrite=False)
    assert written == [
        binary_root / "mixed" / "binary_interaction" / "k_ij.csv",
        binary_root / "mixed" / "binary_interaction" / "l_ij.csv",
        binary_root / "mixed" / "binary_interaction" / "k_hb_ij.csv",
    ]
    for filename, expected in (("k_ij.csv", "-0.06167"), ("l_ij.csv", "0.0123"), ("k_hb_ij.csv", "-0.0345")):
        with (binary_root / "mixed" / "binary_interaction" / filename).open(
            "r", encoding="utf-8-sig", newline=""
        ) as handle:
            rows = list(csv.reader(handle))
        header = rows[0]
        h2o_col = header.index("H2O")
        ethanol_col = header.index("Ethanol")
        h2o_row = next(row for row in rows if row[0] == "H2O")
        ethanol_row = next(row for row in rows if row[0] == "Ethanol")
        assert h2o_row[ethanol_col] == expected
        assert ethanol_row[h2o_col] == expected
    with pytest.raises(InputError, match="overwrite"):
        write_fit_result(binary_result, binary_root, overwrite=False)


def test_public_generic_derivative_evaluator_rejects_removed_backend_names():
    with pytest.raises(InputError, match="jacobian_backend"):
        evaluate_generic_regression_derivatives(
            fixed_payloads=[{"m": [1.0, 1.0], "s": [3.0, 3.0], "e": [200.0, 200.0]}],
            native_records=[],
            optimization_names=("k_ij",),
            species=("H2O", "Ethanol"),
            x=(-0.01,),
            jacobian_backend="finite" + "_difference",
        )


def test_public_generic_derivative_evaluator_rejects_auto_without_autodiff():
    with pytest.raises(InputError, match="backend_unavailable"):
        evaluate_generic_regression_derivatives(
            fixed_payloads=[{"m": [1.0, 1.0], "s": [3.0, 3.0], "e": [200.0, 200.0]}],
            native_records=[
                {
                    "term_name": "binary_vle_fugacity_balance",
                    "term": regression_module.NATIVE_TERM_KINDS["binary_vle_fugacity_balance"],
                    "T": 330.0,
                    "P": 101325.0,
                    "x": [0.7, 0.3],
                    "y": [0.5, 0.5],
                }
            ],
            optimization_names=("k_ij", "l_ij"),
            species=("H2O", "Ethanol"),
            pair=("H2O", "Ethanol"),
            x=(-0.01, 0.02),
        )
