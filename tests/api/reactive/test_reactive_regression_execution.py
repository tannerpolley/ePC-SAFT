from __future__ import annotations

import numpy as np
import pytest

import epcsaft


def _tiny_base_parameters() -> dict[str, np.ndarray]:
    return {
        "m": np.asarray([1.0, 1.2], dtype=float),
        "s": np.asarray([3.0, 3.5], dtype=float),
        "e": np.asarray([200.0, 240.0], dtype=float),
    }

def test_fit_reactive_electrolyte_parameters_validates_then_requires_native_ceres(monkeypatch) -> None:
    called_solver = False

    def fail_if_called(**_kwargs):
        nonlocal called_solver
        called_solver = True
        raise AssertionError("fit route must not evaluate residual rows before native Ceres owns optimization")

    monkeypatch.setattr("epcsaft.reactive_electrolyte.solve_reactive_electrolyte_bubble", fail_if_called)
    batch = epcsaft.ReactiveElectrolyteBatch(
        species=["A", "B"],
        rows=[
            epcsaft.ReactiveElectrolyteRow(
                row_id="row1",
                T=298.15,
                P=101325.0,
                totals={"A": 0.2, "B": 0.8},
                initial_x=[0.2, 0.8],
                balances={"a_total": {"A": 1.0}, "b_total": {"B": 1.0}},
                reactions=[],
                vapor_species=["A", "B"],
                target_pressure=100200.0,
                target_speciation={"A": 0.202},
                source="train",
                split="fit",
            )
        ],
        balances={"a_total": {"A": 1.0}, "b_total": {"B": 1.0}},
        reactions=[],
        vapor_species=["A", "B"],
        base_parameters=_tiny_base_parameters(),
        options=epcsaft.ReactiveElectrolyteBatchOptions(include_state_outputs=False),
    )

    with pytest.raises(epcsaft.InputError, match="native Ceres optimizer"):
        epcsaft.fit_reactive_electrolyte_parameters(
            batch,
            initial_parameters={"A.sigma": 2.8},
            lower_bounds={"A.sigma": 2.5},
            upper_bounds={"A.sigma": 3.05},
            max_iterations=4,
            tolerance=1.0e-10,
        )
    assert called_solver is False

def test_fit_reactive_electrolyte_parameters_rejects_invalid_parameter_inputs(monkeypatch) -> None:
    monkeypatch.setattr(
        "epcsaft.reactive_electrolyte.solve_reactive_electrolyte_bubble",
        lambda **kwargs: epcsaft.ReactiveElectrolyteBubbleResult(
            success=True,
            message="converged",
            x_liq={"A": 0.2, "B": 0.8},
            activity_coefficients={"A": 1.0, "B": 1.0},
            mass_balance_residuals={"total": 0.0},
            charge_residual=0.0,
            reaction_residuals=[],
            named_reaction_residuals={},
            P_total=100000.0,
            y_vap={"A": 0.3, "B": 0.7},
            partial_pressures={"A": 30000.0},
            fugacity_residual={"A": 0.0},
            fugacity_residual_norm=1.0e-9,
            state_failure_count=0,
            penalty_residuals=[],
            diagnostics={},
        ),
    )

    batch = epcsaft.ReactiveElectrolyteBatch(
        species=["A", "B"],
        rows=[
            epcsaft.ReactiveElectrolyteRow(
                row_id="row1",
                T=298.15,
                P=101325.0,
                totals={"A": 0.2, "B": 0.8},
                initial_x=[0.2, 0.8],
                balances={"a_total": {"A": 1.0}, "b_total": {"B": 1.0}},
                reactions=[],
                vapor_species=["A", "B"],
                target_pressure=100000.0,
                target_speciation={"A": 0.2},
            )
        ],
        balances={"a_total": {"A": 1.0}, "b_total": {"B": 1.0}},
        reactions=[],
        vapor_species=["A", "B"],
        base_parameters=_tiny_base_parameters(),
        options=epcsaft.ReactiveElectrolyteBatchOptions(include_state_outputs=False),
    )

    with pytest.raises(epcsaft.InputError, match="at least one fitted parameter"):
        epcsaft.fit_reactive_electrolyte_parameters(batch, initial_parameters={})
    with pytest.raises(epcsaft.InputError, match="unknown parameters"):
        epcsaft.fit_reactive_electrolyte_parameters(
            batch,
            initial_parameters={"A.sigma": 3.0},
            lower_bounds={"B.sigma": 3.2},
        )
    with pytest.raises(epcsaft.InputError, match="inconsistent"):
        epcsaft.fit_reactive_electrolyte_parameters(
            batch,
            initial_parameters={"A.sigma": 3.0},
            lower_bounds={"A.sigma": 3.2},
            upper_bounds={"A.sigma": 3.1},
        )

def test_evaluate_reactive_regression_objective_accepts_speciation_rows(monkeypatch) -> None:
    def fake_speciation(**kwargs):
        mix = kwargs["mixture_factory"](kwargs["initial_x"], kwargs["T"], kwargs["P"])
        sigma = float(np.asarray(mix._params["s"], dtype=float)[0])
        x_a = 0.2 + 0.01 * (sigma - 3.0)
        gamma_a = 1.1 + 0.1 * (sigma - 3.0)
        return epcsaft.ReactiveSpeciationResult(
            success=True,
            message="converged",
            x={"A": x_a, "B": 1.0 - x_a},
            activity_coefficients={"A": gamma_a, "B": 1.0},
            mass_balance_residuals={"total": 0.0},
            charge_residual=0.0,
            reaction_residuals=[],
            named_reaction_residuals={},
            state_failure_count=0,
            diagnostics={},
        )

    monkeypatch.setattr("epcsaft.reactive_speciation.solve_reactive_speciation", fake_speciation)

    batch = epcsaft.ReactiveElectrolyteBatch(
        species=["A", "B"],
        rows=[
            epcsaft.ReactiveElectrolyteRow(
                row_id="row1",
                T=298.15,
                P=101325.0,
                totals={"A": 0.2, "B": 0.8},
                initial_x=[0.2, 0.8],
                balances={"a_total": {"A": 1.0}, "b_total": {"B": 1.0}},
                reactions=[],
                target_speciation={"A": 0.203},
                target_activity={"A": 1.13},
                source="validation",
                split="holdout",
                mode="speciation",
            )
        ],
        balances={"a_total": {"A": 1.0}, "b_total": {"B": 1.0}},
        reactions=[],
        base_parameters=_tiny_base_parameters(),
        options=epcsaft.ReactiveElectrolyteBatchOptions(include_state_outputs=False),
    )

    result = epcsaft.evaluate_reactive_regression_objective(
        batch,
        parameter_map={"A.sigma": 2.9},
    )

    summary = epcsaft.summarize_regression_result(result)
    assert result.batch_result.success_count == 1
    assert "validation" in summary["by_source"]
    assert "holdout" in summary["train_validation"]
    assert summary["fit_success"] is None
    assert "covariance_available" not in summary
