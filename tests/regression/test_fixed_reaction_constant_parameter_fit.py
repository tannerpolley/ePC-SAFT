from __future__ import annotations

import math

import numpy as np

import epcsaft


def test_fixed_reaction_constants_remain_batch_inputs_not_fit_parameters() -> None:
    reaction = epcsaft.ReactionDefinition.from_literature_constant(
        {"A": -1.0, "B": 1.0},
        log_equilibrium_constant=math.log(2.0),
        name="fixed_a_to_b",
    )
    row = epcsaft.ReactiveElectrolyteRow(
        row_id="fixed-k",
        T=298.15,
        P=101325.0,
        initial_x=[0.5, 0.5],
        balances={"total": {"A": 1.0, "B": 1.0}},
        totals={"total": 1.0},
        reactions=[reaction],
        target_speciation={"A": 0.5, "B": 0.5},
        mode="speciation",
    )
    batch = epcsaft.ReactiveElectrolyteBatch(
        species=["A", "B"],
        rows=[row],
        balances=row.balances,
        reactions=row.reactions,
        base_parameters={
            "m": np.asarray([1.0, 1.0], dtype=float),
            "s": np.asarray([3.0, 3.0], dtype=float),
            "e": np.asarray([200.0, 200.0], dtype=float),
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

    objective = context.evaluate_objective({"A.sigma": 3.0})

    assert batch.reactions[0].metadata["fitting_role"] == "fixed_input"
    assert objective.batch_result.success_count + objective.batch_result.failure_count == 1
    assert all("reaction_constant" not in name for name in objective.residual_names)
