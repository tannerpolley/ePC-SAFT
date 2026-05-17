from __future__ import annotations

import numpy as np

import epcsaft


def _tiny_base_parameters() -> dict[str, np.ndarray]:
    return {
        "m": np.asarray([1.0, 1.2], dtype=float),
        "s": np.asarray([3.0, 3.5], dtype=float),
        "e": np.asarray([200.0, 240.0], dtype=float),
    }

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
