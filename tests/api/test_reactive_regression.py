from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

import epcsaft
import epcsaft.reactive_regression as reactive_regression


def _tiny_base_parameters() -> dict[str, np.ndarray]:
    return {
        "m": np.asarray([1.0, 1.2], dtype=float),
        "s": np.asarray([3.0, 3.5], dtype=float),
        "e": np.asarray([200.0, 240.0], dtype=float),
    }


def test_reactive_regression_context_runs_native_speciation_objective_and_jacobian() -> None:
    row = epcsaft.ReactiveElectrolyteRow(
        row_id="native-row",
        T=298.15,
        P=101325.0,
        initial_x=[0.9, 0.05, 0.025, 0.025],
        balances={
            "water": {"H2O": 1.0},
            "sodium": {"NaCl": 1.0, "Na+": 1.0},
            "chloride": {"NaCl": 1.0, "Cl-": 1.0},
        },
        totals={"water": 0.9, "sodium": 0.075, "chloride": 0.075},
        reactions=[
            epcsaft.ReactionDefinition(
                {"NaCl": -1.0, "Na+": 1.0, "Cl-": 1.0},
                np.log(1.0e-2),
                name="salt_dissociation",
                standard_state="ideal_mole_fraction",
            )
        ],
        mode="speciation",
    )
    batch = epcsaft.ReactiveElectrolyteBatch(
        species=["H2O", "NaCl", "Na+", "Cl-"],
        rows=[row],
        balances=row.balances,
        reactions=row.reactions,
        base_parameters={
            "m": np.asarray([1.2047, 1.0, 1.0, 1.0], dtype=float),
            "s": np.asarray([2.7927, 3.1, 2.8232, 2.7560], dtype=float),
            "e": np.asarray([353.95, 230.0, 230.0, 170.0], dtype=float),
            "z": np.asarray([0.0, 0.0, 1.0, -1.0], dtype=float),
            "dielc": np.asarray([78.09, 78.09, 8.0, 8.0], dtype=float),
            "d_born": np.asarray([0.0, 0.0, 3.0, 3.0], dtype=float),
        },
        options=epcsaft.ReactiveElectrolyteBatchOptions(include_state_outputs=False),
    )
    context = epcsaft.ReactiveElectrolyteRegressionContext.from_batch(
        species=batch.species,
        rows=batch.rows,
        balances=batch.balances,
        reactions=batch.reactions,
        base_parameters=batch.base_parameters,
        options=batch.options,
    )

    objective = context.evaluate_objective({"Na+.sigma": 2.8232})
    jacobian = context.finite_difference_jacobian(
        {"Na+.sigma": 2.8232},
        parameters=["Na+.sigma"],
        mode="forward",
        relative_step=1.0e-4,
        log_parameters=False,
    )

    assert objective.batch_result.success_count == 1
    assert objective.batch_result.failure_count == 0
    assert objective.residual_names == ("native-row.reaction.salt_dissociation",)
    assert objective.residuals.shape == (1,)
    assert jacobian.jacobian.shape == (1, 1)


def test_reactive_regression_context_evaluates_batch_and_reuses_warm_starts(monkeypatch) -> None:
    calls: list[dict[str, object]] = []

    def fake_solve(**kwargs):
        calls.append(
            {
                "P_seed": kwargs["P_seed"],
                "initial_x": list(kwargs["initial_x"]),
                "initial_pressure": (
                    None
                    if kwargs["options"] is None or kwargs["options"].bubble_options is None
                    else kwargs["options"].bubble_options.initial_pressure
                ),
                "initial_y_vap": (
                    None
                    if kwargs["options"] is None or kwargs["options"].bubble_options is None
                    else kwargs["options"].bubble_options.initial_y_vap
                ),
            }
        )
        return epcsaft.ReactiveElectrolyteBubbleResult(
            success=True,
            message="converged",
            x_liq={"A": 0.2, "B": 0.8},
            activity_coefficients={"A": 1.1, "B": 0.95},
            mass_balance_residuals={"total": 0.0},
            charge_residual=0.0,
            reaction_residuals=[],
            named_reaction_residuals={},
            P_total=120000.0 + 1000.0 * len(calls),
            y_vap={"A": 0.3, "B": 0.7},
            partial_pressures={"A": 30000.0},
            fugacity_residual={"A": 0.0},
            fugacity_residual_norm=1.0e-9,
            state_failure_count=0,
            penalty_residuals=[],
            diagnostics={},
        )

    monkeypatch.setattr("epcsaft.reactive_electrolyte.solve_reactive_electrolyte_bubble", fake_solve)

    batch = epcsaft.ReactiveElectrolyteBatch(
        species=["A", "B"],
        rows=[
            epcsaft.ReactiveElectrolyteRow(
                row_id="row1",
                T=298.15,
                P_seed=101325.0,
                totals={"A": 0.2, "B": 0.8},
                initial_x=[0.2, 0.8],
                balances={"a_total": {"A": 1.0}, "b_total": {"B": 1.0}},
                reactions=[],
                vapor_species=["A", "B"],
                target_partial_pressures={"A": 30000.0},
                target_speciation={"A": 0.2},
            ),
            epcsaft.ReactiveElectrolyteRow(
                row_id="row2",
                T=298.15,
                P_seed=95000.0,
                totals={"A": 0.21, "B": 0.79},
                initial_x=[0.21, 0.79],
                balances={"a_total": {"A": 1.0}, "b_total": {"B": 1.0}},
                reactions=[],
                vapor_species=["A", "B"],
                target_partial_pressures={"A": 30000.0},
                target_speciation={"A": 0.2},
            ),
        ],
        balances={"a_total": {"A": 1.0}, "b_total": {"B": 1.0}},
        reactions=[],
        vapor_species=["A", "B"],
        mixture_factory=lambda x, T, P: None,
        options=epcsaft.ReactiveElectrolyteBatchOptions(
            warm_start_rows=True,
            warm_start_objective=True,
            penalty_value=8.0,
            include_state_outputs=False,
        ),
    )
    context = epcsaft.ReactiveElectrolyteRegressionContext.from_batch(
        species=batch.species,
        rows=batch.rows,
        balances=batch.balances,
        reactions=batch.reactions,
        options=batch.options,
        vapor_species=batch.vapor_species,
        mixture_factory=batch.mixture_factory,
    )

    first = context.evaluate()
    second = context.evaluate()

    assert first.success_count == 2
    assert first.failure_count == 0
    assert first.residual_names == (
        "row1.partial_pressure.A",
        "row1.x.A",
        "row2.partial_pressure.A",
        "row2.x.A",
    )
    assert first.residuals.shape == (4,)
    assert first.row_results[0].cache_stats["warm_start_source"] == "user_initial"
    assert first.row_results[1].cache_stats["warm_start_source"] == "previous_row"
    assert second.row_results[0].cache_stats["warm_start_source"] == "objective_cache"
    assert second.row_results[1].cache_stats["warm_start_source"] == "objective_cache"
    assert calls[1]["initial_pressure"] == pytest.approx(121000.0)
    assert calls[1]["initial_y_vap"] == pytest.approx({"A": 0.3, "B": 0.7})
    assert first.cache_stats["context_cache_hits"] >= 1


def test_reactive_regression_objective_and_jacobian_are_consistent(monkeypatch) -> None:
    def fake_solve(**kwargs):
        mix = kwargs["mixture_factory"](kwargs["initial_x"], kwargs["T"], kwargs["P_seed"])
        sigma = float(np.asarray(mix._params["s"], dtype=float)[0])
        pressure = 100000.0 + 1000.0 * (sigma - 3.0)
        x_a = 0.2 + 0.01 * (sigma - 3.0)
        return epcsaft.ReactiveElectrolyteBubbleResult(
            success=True,
            message="converged",
            x_liq={"A": x_a, "B": 1.0 - x_a},
            activity_coefficients={"A": 1.0, "B": 1.0},
            mass_balance_residuals={"total": 0.0},
            charge_residual=0.0,
            reaction_residuals=[],
            named_reaction_residuals={},
            P_total=pressure,
            y_vap={"A": 0.3, "B": 0.7},
            partial_pressures={"A": 30000.0},
            fugacity_residual={"A": 0.0},
            fugacity_residual_norm=1.0e-9,
            state_failure_count=0,
            penalty_residuals=[],
            diagnostics={},
        )

    monkeypatch.setattr("epcsaft.reactive_electrolyte.solve_reactive_electrolyte_bubble", fake_solve)

    batch = epcsaft.ReactiveElectrolyteBatch(
        species=["A", "B"],
        rows=[
            epcsaft.ReactiveElectrolyteRow(
                row_id="row1",
                T=298.15,
                P_seed=101325.0,
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
    context = epcsaft.ReactiveElectrolyteRegressionContext.from_batch(
        species=batch.species,
        rows=batch.rows,
        balances=batch.balances,
        reactions=batch.reactions,
        options=batch.options,
        vapor_species=batch.vapor_species,
        base_parameters=batch.base_parameters,
    )

    objective = context.evaluate_objective({"A.sigma": 3.0})
    jacobian = context.finite_difference_jacobian(
        {"A.sigma": 3.0},
        parameters=["A.sigma"],
        mode="central",
        relative_step=1.0e-5,
        log_parameters=False,
    )

    assert objective.residuals.shape == (2,)
    assert jacobian.jacobian.shape == (2, 1)
    np.testing.assert_allclose(jacobian.gradient, jacobian.jacobian.T @ objective.residuals, rtol=1.0e-8, atol=1.0e-8)


def test_reactive_regression_reporting_helpers_write_outputs(monkeypatch, tmp_path: Path) -> None:
    def fake_solve(**kwargs):
        return epcsaft.ReactiveElectrolyteBubbleResult(
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
        )

    monkeypatch.setattr("epcsaft.reactive_electrolyte.solve_reactive_electrolyte_bubble", fake_solve)

    batch = epcsaft.ReactiveElectrolyteBatch(
        species=["A", "B"],
        rows=[
            epcsaft.ReactiveElectrolyteRow(
                row_id="row1",
                T=298.15,
                P_seed=101325.0,
                totals={"A": 0.2, "B": 0.8},
                initial_x=[0.2, 0.8],
                balances={"a_total": {"A": 1.0}, "b_total": {"B": 1.0}},
                reactions=[],
                vapor_species=["A", "B"],
                target_pressure=100000.0,
                target_speciation={"A": 0.2},
                source="train",
                split="fit",
            )
        ],
        balances={"a_total": {"A": 1.0}, "b_total": {"B": 1.0}},
        reactions=[],
        vapor_species=["A", "B"],
        mixture_factory=lambda x, T, P: None,
        options=epcsaft.ReactiveElectrolyteBatchOptions(include_state_outputs=False),
    )
    result = epcsaft.evaluate_reactive_regression_objective(batch)

    summary = tmp_path / "summary.json"
    rows_csv = tmp_path / "rows.csv"
    residuals_csv = tmp_path / "residuals.csv"
    params_csv = tmp_path / "params.csv"

    epcsaft.write_regression_summary(result, summary)
    epcsaft.write_regression_row_table(result, rows_csv)
    epcsaft.write_regression_residual_table(result, residuals_csv)
    epcsaft.write_regression_parameter_table({"A.sigma": 3.0}, params_csv, seed_map={"A.sigma": 2.9})

    assert summary.exists()
    assert rows_csv.exists()
    assert residuals_csv.exists()
    assert params_csv.exists()
    assert "objective" in summary.read_text(encoding="utf-8")


def test_fit_reactive_electrolyte_parameters_returns_fit_result_and_applies_bounds(monkeypatch, tmp_path: Path) -> None:
    def fake_solve(**kwargs):
        mix = kwargs["mixture_factory"](kwargs["initial_x"], kwargs["T"], kwargs["P_seed"])
        sigma = float(np.asarray(mix._params["s"], dtype=float)[0])
        pressure = 100000.0 + 1000.0 * (sigma - 3.0)
        x_a = 0.2 + 0.01 * (sigma - 3.0)
        return epcsaft.ReactiveElectrolyteBubbleResult(
            success=True,
            message="converged",
            x_liq={"A": x_a, "B": 1.0 - x_a},
            activity_coefficients={"A": 1.0, "B": 1.0},
            mass_balance_residuals={"total": 0.0},
            charge_residual=0.0,
            reaction_residuals=[],
            named_reaction_residuals={},
            P_total=pressure,
            y_vap={"A": 0.3, "B": 0.7},
            partial_pressures={"A": 30000.0},
            fugacity_residual={"A": 0.0},
            fugacity_residual_norm=1.0e-9,
            state_failure_count=0,
            penalty_residuals=[],
            diagnostics={},
        )

    monkeypatch.setattr("epcsaft.reactive_electrolyte.solve_reactive_electrolyte_bubble", fake_solve)

    batch = epcsaft.ReactiveElectrolyteBatch(
        species=["A", "B"],
        rows=[
            epcsaft.ReactiveElectrolyteRow(
                row_id="row1",
                T=298.15,
                P_seed=101325.0,
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

    fit = epcsaft.fit_reactive_electrolyte_parameters(
        batch,
        initial_parameters={"A.sigma": 2.8},
        lower_bounds={"A.sigma": 2.5},
        upper_bounds={"A.sigma": 3.05},
        max_iterations=4,
        tolerance=1.0e-10,
        log_parameters=False,
    )

    assert isinstance(fit, epcsaft.ReactiveRegressionFitResult)
    assert fit.parameter_map["A.sigma"] == pytest.approx(3.05)
    assert fit.active_bounds["A.sigma"] is True
    assert fit.covariance_available is True
    summary = epcsaft.summarize_regression_result(fit)
    assert summary["fit_success"] is True
    assert summary["upper_bounds"]["A.sigma"] == pytest.approx(3.05)
    assert summary["covariance_status"] == "available"

    summary_path = tmp_path / "fit_summary.json"
    rows_path = tmp_path / "fit_rows.csv"
    residuals_path = tmp_path / "fit_residuals.csv"
    params_path = tmp_path / "fit_params.csv"
    epcsaft.write_regression_summary(fit, summary_path)
    epcsaft.write_regression_row_table(fit, rows_path)
    epcsaft.write_regression_residual_table(fit, residuals_path)
    epcsaft.write_regression_parameter_table(fit, params_path)

    assert '"fit_success": true' in summary_path.read_text(encoding="utf-8")
    assert "active_bound" in params_path.read_text(encoding="utf-8")
    assert rows_path.exists()
    assert residuals_path.exists()


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
                P_seed=101325.0,
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


def test_fit_reactive_electrolyte_parameters_accepts_speciation_rows(monkeypatch) -> None:
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

    fit = epcsaft.fit_reactive_electrolyte_parameters(
        batch,
        initial_parameters={"A.sigma": 2.9},
        max_iterations=3,
        tolerance=1.0e-12,
        log_parameters=False,
    )

    summary = epcsaft.summarize_regression_result(fit)
    assert isinstance(fit, epcsaft.ReactiveRegressionFitResult)
    assert fit.objective_result.batch_result.success_count == 1
    assert "validation" in summary["by_source"]
    assert "holdout" in summary["train_validation"]
    assert isinstance(summary["fit_success"], bool)


def test_fit_reactive_electrolyte_parameters_reports_nonconverged_fit(monkeypatch) -> None:
    def fake_solve(**kwargs):
        mix = kwargs["mixture_factory"](kwargs["initial_x"], kwargs["T"], kwargs["P_seed"])
        sigma = float(np.asarray(mix._params["s"], dtype=float)[0])
        pressure = 100000.0 + 1000.0 * (sigma - 3.0)
        x_a = 0.2 + 0.01 * (sigma - 3.0)
        return epcsaft.ReactiveElectrolyteBubbleResult(
            success=True,
            message="converged",
            x_liq={"A": x_a, "B": 1.0 - x_a},
            activity_coefficients={"A": 1.0, "B": 1.0},
            mass_balance_residuals={"total": 0.0},
            charge_residual=0.0,
            reaction_residuals=[],
            named_reaction_residuals={},
            P_total=pressure,
            y_vap={"A": 0.3, "B": 0.7},
            partial_pressures={"A": 30000.0},
            fugacity_residual={"A": 0.0},
            fugacity_residual_norm=1.0e-9,
            state_failure_count=0,
            penalty_residuals=[],
            diagnostics={},
        )

    monkeypatch.setattr("epcsaft.reactive_electrolyte.solve_reactive_electrolyte_bubble", fake_solve)

    batch = epcsaft.ReactiveElectrolyteBatch(
        species=["A", "B"],
        rows=[
            epcsaft.ReactiveElectrolyteRow(
                row_id="row1",
                T=298.15,
                P_seed=101325.0,
                totals={"A": 0.2, "B": 0.8},
                initial_x=[0.2, 0.8],
                balances={"a_total": {"A": 1.0}, "b_total": {"B": 1.0}},
                reactions=[],
                vapor_species=["A", "B"],
                target_pressure=100600.0,
                target_speciation={"A": 0.206},
            )
        ],
        balances={"a_total": {"A": 1.0}, "b_total": {"B": 1.0}},
        reactions=[],
        vapor_species=["A", "B"],
        base_parameters=_tiny_base_parameters(),
        options=epcsaft.ReactiveElectrolyteBatchOptions(include_state_outputs=False),
    )

    fit = epcsaft.fit_reactive_electrolyte_parameters(
        batch,
        initial_parameters={"A.sigma": 2.8},
        max_iterations=1,
        damping=0.25,
        tolerance=1.0e-12,
        log_parameters=False,
    )

    assert fit.success is False
    assert fit.message == "maximum iterations reached"
    assert fit.iterations == 1


def test_reactive_regression_legacy_wrapper_keeps_fixed_shape(monkeypatch) -> None:
    def fake_solve(**kwargs):
        return epcsaft.ReactiveElectrolyteBubbleResult(
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
        )

    monkeypatch.setattr("epcsaft.reactive_electrolyte.solve_reactive_electrolyte_bubble", fake_solve)

    result = epcsaft.evaluate_reactive_electrolyte_bubble_residuals(
        [
            {
                "row_id": "row1",
                "T": 298.15,
                "P_seed": 101325.0,
                "totals": {"A": 0.2, "B": 0.8},
                "initial_x": [0.2, 0.8],
                "target_partial_pressures": {"A": 30000.0},
                "target_x": {"A": 0.2},
            }
        ],
        species=["A", "B"],
        mixture_factory=lambda x, T, P: None,
        balances={"a_total": {"A": 1.0}, "b_total": {"B": 1.0}},
        reactions=[],
        vapor_species=["A", "B"],
        pressure_species=["A"],
        speciation_species=["A"],
    )

    assert result.success_count == 1
    assert result.failure_count == 0
    assert result.residuals.shape == (2,)
    assert result.residual_names == ("row1.partial_pressure.A", "row1.x.A")
