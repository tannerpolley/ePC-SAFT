from __future__ import annotations

import numpy as np

import epcsaft


def _toy_reactive_phase_case() -> tuple[
    epcsaft.ePCSAFTMixture,
    np.ndarray,
    dict[str, object],
    epcsaft.ReactionDefinition,
]:
    mix = epcsaft.ePCSAFTMixture.from_params(
        {
            "MW": np.asarray([32.042e-3, 84.147e-3]),
            "m": np.asarray([1.5255, 2.5303]),
            "s": np.asarray([3.2300, 3.8499]),
            "e": np.asarray([188.90, 278.11]),
            "e_assoc": np.asarray([2899.5, 0.0]),
            "vol_a": np.asarray([0.035176, 0.0]),
            "assoc_scheme": ["2B", None],
            "k_ij": np.asarray([[0.0, 0.051], [0.051, 0.0]]),
            "z": np.asarray([0.0, 0.0]),
            "dielc": np.asarray([33.05, 2.02]),
        },
        species=["Methanol", "Cyclohexane"],
    )
    liq1 = np.asarray([0.11757838279937723, 0.8824216172006228])
    liq2 = np.asarray([0.7985874308392054, 0.20141256916079467])
    beta2 = 0.48813098468607985
    feed = (1.0 - beta2) * liq1 + beta2 * liq2
    reaction = epcsaft.ReactionDefinition.from_literature_constant(
        {"Methanol": -1.0, "Cyclohexane": 1.0},
        log_equilibrium_constant=-0.079259405371,
        name="methanol_to_cyclohexane",
        standard_state="mole_fraction_activity",
        source="public route smoke fixture",
    )
    return mix, feed, {"liq1": liq1, "liq2": liq2, "phase_fraction": beta2}, reaction


def test_reactive_phase_equilibrium_problem_solves_with_native_coupled_route(monkeypatch) -> None:
    mix, feed, initial_phases, reaction = _toy_reactive_phase_case()

    def fail_if_staged(*_args, **_kwargs):
        raise AssertionError("production ReactivePhaseEquilibriumProblem must not call the staged helper")

    monkeypatch.setattr(mix, "reactive_staged_equilibrium", fail_if_staged)
    problem = epcsaft.ReactivePhaseEquilibriumProblem(
        T=298.15,
        P=1.013e5,
        z=feed,
        balances={"total": {"Methanol": 1.0, "Cyclohexane": 1.0}},
        totals={"total": 1.0},
        reactions=[reaction],
        phase_kind="lle_flash",
        phase_options=epcsaft.EquilibriumOptions(max_iterations=80, tolerance=1.0e-8, min_composition=1.0e-12),
        phase_kwargs={"initial_phases": initial_phases},
    )

    result = mix.solve_equilibrium(problem)

    assert isinstance(result, epcsaft.EquilibriumResult)
    assert result.problem_kind == "reactive_phase_equilibrium"
    assert result.diagnostics["reactive_workflow_class"] == "coupled_native"
    assert result.diagnostics["staged_route_used"] is False
    assert result.diagnostics["solver_method"] == "ceres_trust_region_coupled_reactive_phase_equilibrium"
    assert result.diagnostics["jacobian_available"] is True
    assert result.diagnostics["reaction_residual_norm"] <= 1.0e-8


def test_reactive_phase_equilibrium_problem_rejects_non_lle_production_kind() -> None:
    mix, feed, _initial_phases, reaction = _toy_reactive_phase_case()
    problem = epcsaft.ReactivePhaseEquilibriumProblem(
        T=298.15,
        P=1.013e5,
        z=feed,
        balances={"total": {"Methanol": 1.0, "Cyclohexane": 1.0}},
        totals={"total": 1.0},
        reactions=[reaction],
        phase_kind="tp_flash",
    )

    with np.testing.assert_raises_regex(epcsaft.InputError, "reactive_staged_equilibrium"):
        mix.solve_equilibrium(problem)
