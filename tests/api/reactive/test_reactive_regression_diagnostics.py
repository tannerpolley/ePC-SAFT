from __future__ import annotations

from pathlib import Path

import numpy as np

import epcsaft


def _tiny_base_parameters() -> dict[str, np.ndarray]:
    return {
        "m": np.asarray([1.0, 1.2], dtype=float),
        "s": np.asarray([3.0, 3.5], dtype=float),
        "e": np.asarray([200.0, 240.0], dtype=float),
    }

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
        P_seed=101325.0,
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
        options=epcsaft.ReactiveElectrolyteBatchOptions(
            include_state_outputs=False,
            warm_start_rows=True,
            warm_start_objective=True,
        ),
        reactive_bubble_options=epcsaft.ReactiveElectrolyteBubbleOptions(error_mode="result"),
    )
    return batch, water_sigma

def test_reactive_regression_reporting_helpers_write_outputs(monkeypatch, tmp_path: Path) -> None:
    def fake_solve(**_kwargs):
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

def test_reactive_bubble_residual_wrapper_keeps_fixed_shape(monkeypatch) -> None:
    def fake_solve(**_kwargs):
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
