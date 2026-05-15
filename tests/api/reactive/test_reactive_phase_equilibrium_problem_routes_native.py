from __future__ import annotations

import math

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
            "m": np.asarray([1.0, 1.0]),
            "s": np.asarray([3.0, 3.0]),
            "e": np.asarray([200.0, 200.0]),
            "z": np.asarray([0.0, 0.0]),
            "dielc": np.asarray([2.0, 2.0]),
        },
        species=["A", "B"],
    )
    feed = np.asarray([0.4, 0.6], dtype=float)
    state = mix.state(T=298.15, P=1.0e5, x=feed, phase="liq")
    ln_phi = state.fugacity_coefficient(natural_log=True)
    log_k = -math.log(feed[0]) - float(ln_phi[0]) + math.log(feed[1]) + float(ln_phi[1])
    reaction = epcsaft.ReactionDefinition.from_literature_constant(
        {"A": -1.0, "B": 1.0},
        log_equilibrium_constant=log_k,
        name="a_to_b",
        standard_state="mole_fraction_activity",
        source="public route smoke fixture",
    )
    return mix, feed, {"liq1": feed, "liq2": feed, "phase_fraction": 0.5}, reaction


def test_reactive_phase_equilibrium_problem_solves_with_native_coupled_route(monkeypatch) -> None:
    mix, feed, initial_phases, reaction = _toy_reactive_phase_case()

    def fail_if_staged(*_args, **_kwargs):
        raise AssertionError("production ReactivePhaseEquilibriumProblem must not call the staged helper")

    monkeypatch.setattr(mix, "reactive_staged_equilibrium", fail_if_staged)
    problem = epcsaft.ReactivePhaseEquilibriumProblem(
        T=298.15,
        P=1.0e5,
        z=feed,
        balances={"total": {"A": 1.0, "B": 1.0}},
        totals={"total": 1.0},
        reactions=[reaction],
        phase_kind="lle_flash",
        phase_options=epcsaft.EquilibriumOptions(max_iterations=20, tolerance=1.0e-10, min_composition=1.0e-12),
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
        P=1.0e5,
        z=feed,
        balances={"total": {"A": 1.0, "B": 1.0}},
        totals={"total": 1.0},
        reactions=[reaction],
        phase_kind="tp_flash",
    )

    with np.testing.assert_raises_regex(epcsaft.InputError, "reactive_staged_equilibrium"):
        mix.solve_equilibrium(problem)
