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
            "not_available_reason": "",
            "hessian_available": False,
            "hessian_backend": "not_implemented",
            "hessian_fallback_used": False,
            "hessian_fallback_reason": "stubbed hessian skeleton",
        }

    runner_name = "_run_native_generic_ceres" if backend == "ceres" else "_run_native_generic_least_squares"
    monkeypatch.setattr(regression_module, runner_name, fake_runner)
    return calls

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
