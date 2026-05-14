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

def _stub_native_generic_runner(monkeypatch, *, backend="least_squares_native"):
    calls = []
    jacobian_backend = "cppad_implicit" if backend == "ceres" else "stub"

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
            "backend": backend,
            "jacobian_available": True,
            "jacobian_backend": jacobian_backend,
            "jacobian_fallback_used": False,
            "jacobian_fallback_reason": "",
            "backend_unavailable_reason": "",
            "hessian_available": False,
            "hessian_backend": "not_implemented",
            "hessian_fallback_used": False,
            "hessian_fallback_reason": "stubbed hessian skeleton",
        }

    runner_name = "_run_native_generic_ceres" if backend == "ceres" else "_run_native_generic_least_squares"
    monkeypatch.setattr(regression_module, runner_name, fake_runner)
    return calls

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

def test_fit_pure_ion_requires_composition_and_supported_records():
    with pytest.raises(InputError, match="composition"):
        epcsaft.fit_pure_ion(
            [{"T": 298.15, "P": 101325.0, "osmotic_coefficient": 0.93}],
            "Na+",
            dataset="2026_Khudaida",
        )

    with pytest.raises(InputError, match=r"density|osmotic|mean-ionic|mean ionic"):
        epcsaft.fit_pure_ion(
            [{"T": 298.15, "P": 101325.0, "molality": 0.1}],
            "Na+",
            dataset="2026_Khudaida",
            species=["H2O", "Na+", "Cl-"],
            solvent="H2O",
        )

def test_fit_pure_ion_default_s_e_bounds_and_multistart_contract(monkeypatch):
    calls = _stub_native_generic_runner(monkeypatch, backend="ceres")
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
    assert result.backend == "ceres"
    assert result.jacobian_available is True
    assert result.jacobian_backend == "cppad_implicit"
    assert result.backend_unavailable_reason == ""
    assert result.hessian_available is False
    assert result.hessian_backend == "not_implemented"
    assert result.problem.mode == "pure_ion"
    assert result.problem.fit_targets == ("s", "e")
    assert set(result.metrics_by_term) == {"osmotic_coefficient", "mean_ionic_activity"}
    assert result.fitted_values == {"s": 2.6, "e": 210.0}
    assert result.parameter_map == {"s": 2.6, "e": 210.0}
    assert result.active_bounds == []
    assert {row["row_family"] for row in result.row_diagnostics} == {"osmotic_coefficient", "mean_ionic_activity"}
    assert result.provenance_report["parameter_movement"] == {"s": 0.0, "e": 0.0}
    assert len(calls) == 1
    assert calls[0]["optimization_names"] == ("s", "e")
    assert calls[0]["component"] == "Na+"
    assert calls[0]["species"] == ("H2O", "Na+", "Cl-")
    assert calls[0]["multistart"] == 3
    assert {record["term_name"] for record in calls[0]["native_records"]} == {
        "osmotic_coefficient",
        "mean_ionic_activity",
    }

def test_fit_pure_ion_accepts_d_born_and_born_user_options(monkeypatch):
    calls = _stub_native_generic_runner(monkeypatch, backend="ceres")
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
    assert result.backend == "ceres"
    assert result.jacobian_backend == "cppad_implicit"
    assert result.problem.fit_targets == ("d_born",)
    assert "osmotic_coefficient" in result.metrics_by_term
    assert len(calls) == 1
    assert calls[0]["optimization_names"] == ("d_born",)
    assert result.provenance_report["parameter_sources"]["Na+.d_born"] == "ion_activity"
    assert result.provenance_report["warnings"] == []

def test_fit_pure_ion_passes_explicit_mean_ionic_pair_label_to_native_backend(monkeypatch):
    calls = _stub_native_generic_runner(monkeypatch, backend="ceres")
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
    assert result.backend == "ceres"
    assert "mean_ionic_activity" in result.metrics_by_term
    assert len(calls) == 1
    mean_ionic_record = next(
        record for record in calls[0]["native_records"] if record["term_name"] == "mean_ionic_activity"
    )
    assert mean_ionic_record["target_index"] == 1
    assert mean_ionic_record["target_index_2"] == 2

def test_fit_binary_pair_vle_kij_default_and_rejects_temperature_models(monkeypatch):
    calls = _stub_native_generic_runner(monkeypatch, backend="ceres")
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
    assert result.backend == "ceres"
    assert result.jacobian_available is True
    assert result.jacobian_backend == "cppad_implicit"
    assert result.hessian_available is False
    assert result.problem.mode == "binary_pair"
    assert result.problem.fit_targets == ("k_ij",)
    assert set(result.fitted_values) == {"k_ij"}
    assert result.parameter_map == result.fitted_values
    assert result.row_diagnostics == [{"row_family": "binary_vle_fugacity_balance", "metric": 0.0}]
    assert result.active_bounds == []
    assert result.metrics_by_term == {"binary_vle_fugacity_balance": 0.0}
    assert result.provenance_report["parameter_movement"] == {"k_ij": 0.0}
    assert result.provenance_report["source_summary"]["record_count"] == 2
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

def test_fit_binary_pair_rejects_unsupported_generic_binary_optimizer_targets(monkeypatch):
    calls = _stub_native_generic_runner(monkeypatch, backend="ceres")
    records = [
        {"T": 330.0, "P": 101325.0, "x_H2O": 0.7, "x_Ethanol": 0.3, "y_H2O": 0.5, "y_Ethanol": 0.5},
        {"T": 340.0, "P": 101325.0, "x_H2O": 0.6, "x_Ethanol": 0.4, "y_H2O": 0.4, "y_Ethanol": 0.6},
    ]

    with pytest.raises(InputError, match="supports only constant k_ij"):
        epcsaft.fit_binary_pair(
            records,
            ("H2O", "Ethanol"),
            dataset="2026_Khudaida",
            fit_targets=("k_ij", "l_ij", "k_hb_ij"),
            initial_guess={"k_ij": -0.02, "l_ij": 0.01, "k_hb_ij": 0.02},
            bounds={"k_ij": (-0.2, 0.2), "l_ij": (-0.2, 0.2), "k_hb_ij": (-0.2, 0.2)},
            multistart=1,
        )

    assert calls == []

def test_fit_binary_pair_rejects_native_least_squares_backend():
    records = [
        {"T": 330.0, "P": 101325.0, "x_H2O": 0.7, "x_Ethanol": 0.3, "y_H2O": 0.5, "y_Ethanol": 0.5},
        {"T": 340.0, "P": 101325.0, "x_H2O": 0.6, "x_Ethanol": 0.4, "y_H2O": 0.4, "y_Ethanol": 0.6},
    ]

    with pytest.raises(InputError, match="native analytic/CppAD/implicit"):
        epcsaft.fit_binary_pair(
            records,
            ("H2O", "Ethanol"),
            dataset="2026_Khudaida",
            initial_guess={"k_ij": -0.02},
            optimizer_backend="least_squares_native",
        )

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
