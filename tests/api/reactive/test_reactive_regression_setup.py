from __future__ import annotations

from dataclasses import MISSING, fields

import numpy as np
import pytest

import epcsaft


def _tiny_base_parameters() -> dict[str, np.ndarray]:
    return {
        "m": np.asarray([1.0, 1.2], dtype=float),
        "s": np.asarray([3.0, 3.5], dtype=float),
        "e": np.asarray([200.0, 240.0], dtype=float),
    }


def test_reactive_electrolyte_regression_public_surfaces_are_current() -> None:
    row_fields = fields(epcsaft.ReactiveElectrolyteRow)
    row_names = tuple(field.name for field in row_fields)
    option_names = tuple(field.name for field in fields(epcsaft.ReactiveElectrolyteBatchOptions))
    result_names = tuple(field.name for field in fields(epcsaft.ReactiveElectrolyteRowResult))

    assert row_names == (
        "row_id",
        "T",
        "P",
        "initial_x",
        "balances",
        "totals",
        "reactions",
        "vapor_species",
        "target_pressure",
        "target_speciation",
        "target_activity",
        "target_fugacity",
        "target_density",
        "target_relative_permittivity",
        "target_partial_pressures",
        "weights",
        "source",
        "split",
        "metadata",
        "mode",
    )
    assert row_fields[2].default is MISSING
    assert option_names == ("penalty_value", "failure_residual_mode", "include_state_outputs")
    assert result_names == (
        "row_id",
        "success",
        "message",
        "composition",
        "pressure",
        "ln_fugacity",
        "activity_coefficients",
        "density",
        "relative_permittivity",
        "residuals",
        "residual_names",
        "failure_diagnostics",
        "active_bounds",
        "elapsed_seconds",
        "partial_pressures",
        "y_vap",
        "named_reaction_residuals",
        "source",
        "split",
        "metadata",
    )


def _native_mixed_pressure_speciation_batch() -> tuple[epcsaft.ReactiveElectrolyteBatch, float]:
    temperature = 298.15
    water_sigma = 2.7927 + 10.11 * np.exp(-0.01775 * temperature) - 1.417 * np.exp(-0.01146 * temperature)
    params = {
        "MW": np.asarray([18.01528e-3, 22.98e-3, 35.45e-3]),
        "m": np.asarray([1.2047, 1.0, 1.0]),
        "s": np.asarray([water_sigma, 2.8232, 2.7560]),
        "e": np.asarray([353.95, 230.0, 170.0]),
        "e_assoc": np.asarray([2425.7, 0.0, 0.0]),
        "vol_a": np.asarray([0.04509, 0.0, 0.0]),
        "assoc_scheme": ["2B", None, None],
        "z": np.asarray([0.0, 1.0, -1.0]),
        "dielc": np.asarray([78.09, 8.0, 8.0]),
        "d_born": np.asarray([0.0, 3.445, 4.1]),
        "f_solv": np.asarray([1.5, 1.0, 1.0]),
        "k_ij": np.asarray([[0.0, 0.0045, -0.25], [0.0045, 0.0, 0.317], [-0.25, 0.317, 0.0]]),
        "l_ij": np.zeros((3, 3)),
        "k_hb": np.zeros((3, 3)),
    }
    balances = {"water": {"water": 1.0}, "sodium": {"Na+": 1.0}, "chloride": {"Cl-": 1.0}}
    row = epcsaft.ReactiveElectrolyteRow(
        row_id="native-mixed",
        T=temperature,
        P=101325.0,
        totals={"water": 0.98, "sodium": 0.01, "chloride": 0.01},
        initial_x=[0.98, 0.01, 0.01],
        balances=balances,
        reactions=[],
        vapor_species=["water"],
        target_partial_pressures={"water": 3000.0},
        target_speciation={"water": 0.98},
    )
    batch = epcsaft.ReactiveElectrolyteBatch(
        species=["water", "Na+", "Cl-"],
        rows=[row],
        balances=balances,
        reactions=[],
        vapor_species=["water"],
        base_parameters=params,
        options=epcsaft.ReactiveElectrolyteBatchOptions(include_state_outputs=False),
        reactive_bubble_options=epcsaft.ReactiveElectrolyteBubbleOptions(error_mode="result"),
    )
    return batch, water_sigma

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
                standard_state="mole_fraction_activity",
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
    assert objective.batch_result.success_count == 0
    assert objective.batch_result.failure_count == 1
    assert objective.residual_names == ("native-row.reaction.salt_dissociation",)
    assert objective.residuals.shape == (1,)
    assert "nonideal reactive speciation requires a native Ipopt Gibbs/activity NLP route builder" in (
        objective.batch_result.row_results[0].message
    )
    with pytest.raises(epcsaft.InputError, match="native Ceres derivative coverage"):
        context.evaluate_derivatives({"Na+.sigma": 2.8232}, parameters=["Na+.sigma"])

def test_reactive_regression_context_evaluates_batch_with_explicit_row_inputs(monkeypatch) -> None:
    calls: list[dict[str, object]] = []

    def fake_solve(**kwargs):
        bubble_options = None if kwargs["options"] is None else kwargs["options"].bubble_options
        calls.append(
            {
                "P": kwargs["P"],
                "initial_x": list(kwargs["initial_x"]),
                "bubble_options": bubble_options,
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
                P=101325.0,
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
                P=95000.0,
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
    assert second.cache_stats["context_evaluations"] >= 2
    assert calls[1]["P"] == pytest.approx(95000.0)
    assert calls[1]["initial_x"] == pytest.approx([0.21, 0.79])
    assert calls[1]["bubble_options"] is None
    assert first.cache_stats["context_evaluations"] >= 1


def test_reactive_regression_objective_and_jacobian_are_consistent(monkeypatch) -> None:
    def fake_solve(**kwargs):
        mix = kwargs["mixture_factory"](kwargs["initial_x"], kwargs["T"], kwargs["P"])
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
    assert objective.residuals.shape == (2,)
    with pytest.raises(epcsaft.InputError, match="native Ceres derivative coverage"):
        context.evaluate_derivatives({"A.sigma": 3.0}, parameters=["A.sigma"])
