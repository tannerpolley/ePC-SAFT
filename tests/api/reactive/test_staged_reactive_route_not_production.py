from __future__ import annotations

import pytest

import epcsaft
from tests.api.reactive.test_reactive_phase_equilibrium_problem_routes_native import _toy_reactive_phase_case


def test_explicit_staged_kind_remains_separate_from_production_reactive_lle(monkeypatch) -> None:
    mix, feed, initial_phases, reaction = _toy_reactive_phase_case()

    def successful_staged(**_kwargs):
        chemical = epcsaft.ReactiveSpeciationResult(
            success=True,
            message="converged",
            x={"Methanol": float(feed[0]), "Cyclohexane": float(feed[1])},
            activity_coefficients={},
            mass_balance_residuals={"total": 0.0},
            charge_residual=0.0,
            reaction_residuals=[0.0],
            named_reaction_residuals={"methanol_to_cyclohexane": 0.0},
            state_failure_count=0,
            diagnostics={"phase_equilibrium_handoff": {}},
        )
        return epcsaft.ReactiveStagedEquilibriumResult(
            success=True,
            message="converged",
            z=chemical.x,
            chemical=chemical,
            phase={},
            diagnostics={
                "reactive_workflow_class": "staged",
                "staged_feed": chemical.x,
            },
        )

    monkeypatch.setattr(mix, "reactive_staged_equilibrium", successful_staged)
    staged = mix.equilibrium(
        kind="reactive_staged",
        T=298.15,
        P=1.013e5,
        z=feed,
        balances={"total": {"Methanol": 1.0, "Cyclohexane": 1.0}},
        totals={"total": 1.0},
        reactions=[reaction],
        phase_kind="lle_flash",
        initial_phases=initial_phases,
    )

    assert isinstance(staged, epcsaft.ReactiveStagedEquilibriumResult)
    assert staged.diagnostics["reactive_workflow_class"] == "staged"
    assert set(staged.diagnostics["staged_feed"]) == {"Methanol", "Cyclohexane"}
    assert sum(staged.diagnostics["staged_feed"].values()) == pytest.approx(1.0, abs=1.0e-12)

    def fail_if_staged(*_args, **_kwargs):
        raise AssertionError("kind='reactive_lle' must bypass explicit staged compatibility")

    monkeypatch.setattr(mix, "reactive_staged_equilibrium", fail_if_staged)
    with pytest.raises(epcsaft.InputError, match="native Ipopt reactive phase-equilibrium NLP route"):
        mix.equilibrium(
            kind="reactive_lle",
            T=298.15,
            P=1.013e5,
            z=feed,
            balances={"total": {"Methanol": 1.0, "Cyclohexane": 1.0}},
            totals={"total": 1.0},
            reactions=[reaction],
            initial_phases=initial_phases,
            phase_options=epcsaft.EquilibriumOptions(max_iterations=80, tolerance=1.0e-8, min_composition=1.0e-12),
        )


def test_reactive_lle_dispatcher_does_not_accept_phase_route_controls() -> None:
    mix, feed, _initial_phases, reaction = _toy_reactive_phase_case()

    with pytest.raises(epcsaft.InputError, match="does not support phase_kwargs"):
        mix.equilibrium(
            kind="reactive_lle",
            T=298.15,
            P=1.013e5,
            z=feed,
            balances={"total": {"Methanol": 1.0, "Cyclohexane": 1.0}},
            totals={"total": 1.0},
            reactions=[reaction],
            parent_phase="liq",
        )
