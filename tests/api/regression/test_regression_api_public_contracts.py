"""Regression API contract tests outside the hydrocarbon benchmark."""

from __future__ import annotations

import numpy as np
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
    assert not hasattr(epcsaft, "fit_reactive_electrolyte" + "_parameters")
    assert not hasattr(epcsaft, "ReactiveRegression" + "FitResult")
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


def test_reactive_regression_summary_uses_shared_target_family_evidence(monkeypatch):
    def fake_speciation(**kwargs):
        mix = kwargs["mixture_factory"](kwargs["initial_x"], kwargs["T"], kwargs["P"])
        sigma = float(np.asarray(mix._params["s"], dtype=float)[0])
        return epcsaft.ReactiveSpeciationResult(
            success=True,
            message="converged",
            x={"A": 0.2 + 0.01 * (sigma - 3.0), "B": 0.8},
            activity_coefficients={"A": 1.1, "B": 1.0},
            mass_balance_residuals={},
            charge_residual=0.0,
            reaction_residuals=[],
            named_reaction_residuals={},
            diagnostics={},
        )

    monkeypatch.setattr("epcsaft.reactive_speciation.solve_reactive_speciation", fake_speciation)
    params = epcsaft.ParameterSet.from_records(
        [
            epcsaft.PureRecord("A", molar_mass=10.0e-3, m=1.0, sigma=3.0, epsilon_k=200.0),
            epcsaft.PureRecord("B", molar_mass=20.0e-3, m=1.2, sigma=3.5, epsilon_k=240.0),
        ],
        metadata={"dataset": "unit-parameters"},
    )
    batch = epcsaft.ReactiveElectrolyteBatch(
        species=["A", "B"],
        rows=[
            epcsaft.ReactiveElectrolyteRow(
                row_id="row1",
                T=298.15,
                P=101325.0,
                initial_x=[0.2, 0.8],
                balances={"a_total": {"A": 1.0}, "b_total": {"B": 1.0}},
                reactions=[],
                target_speciation={"A": 0.201},
                target_activity={"A": 1.1},
                source="validation",
                split="holdout",
                mode="speciation",
            )
        ],
        balances={"a_total": {"A": 1.0}, "b_total": {"B": 1.0}},
        reactions=[],
        base_parameters=params,
        options=epcsaft.ReactiveElectrolyteBatchOptions(include_state_outputs=False),
    )

    result = epcsaft.evaluate_reactive_regression_objective(batch, parameter_map={"A.sigma": 3.1})
    summary = epcsaft.summarize_regression_result(result)

    assert summary["target_family_summaries"] == result.batch_result.diagnostics["target_family_summaries"]
    assert summary["target_family_summaries"]["speciation"]["residual_count"] == 1
    assert summary["target_family_summaries"]["activity"]["residual_count"] == 1
    assert summary["residual_block_norms"]["speciation"] == pytest.approx(
        result.batch_result.diagnostics["target_family_summaries"]["speciation"]["residual_block_norm"]
    )
