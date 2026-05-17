"""Regression API contract tests outside the hydrocarbon benchmark."""

from __future__ import annotations

import pytest

import epcsaft


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
    assert hasattr(epcsaft, "evaluate_reactive_regression_objective")
    assert hasattr(epcsaft, "ReactiveRegressionObjectiveResult")
    assert not hasattr(epcsaft, "evaluate_reactive_electrolyte" + "_bubble_residuals")
    assert not hasattr(epcsaft, "Reactive" + "ElectrolyteRegressionResult")
    assert hasattr(epcsaft, "TargetRow")
    assert hasattr(epcsaft, "TargetDataset")


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
