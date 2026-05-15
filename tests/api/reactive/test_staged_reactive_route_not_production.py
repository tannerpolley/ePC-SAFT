from __future__ import annotations

import numpy as np
import pytest

import epcsaft
from tests.api.reactive.test_reactive_phase_equilibrium_problem_routes_native import _toy_reactive_phase_case


def test_explicit_staged_kind_remains_separate_from_production_reactive_lle(monkeypatch) -> None:
    mix, feed, initial_phases, reaction = _toy_reactive_phase_case()

    staged = mix.equilibrium(
        kind="reactive_staged",
        T=298.15,
        P=1.0e5,
        z=feed,
        balances={"total": {"A": 1.0, "B": 1.0}},
        totals={"total": 1.0},
        reactions=[reaction],
        phase_kind="lle_flash",
        initial_phases=initial_phases,
    )

    assert isinstance(staged, epcsaft.ReactiveStagedEquilibriumResult)
    assert staged.diagnostics["reactive_workflow_class"] == "staged"
    assert staged.diagnostics["staged_feed"] == pytest.approx({"A": 0.4, "B": 0.6})

    def fail_if_staged(*_args, **_kwargs):
        raise AssertionError("kind='reactive_lle' must bypass explicit staged compatibility")

    monkeypatch.setattr(mix, "reactive_staged_equilibrium", fail_if_staged)
    production = mix.equilibrium(
        kind="reactive_lle",
        T=298.15,
        P=1.0e5,
        z=feed,
        balances={"total": {"A": 1.0, "B": 1.0}},
        totals={"total": 1.0},
        reactions=[reaction],
        initial_phases=initial_phases,
        phase_options=epcsaft.EquilibriumOptions(max_iterations=20, tolerance=1.0e-10, min_composition=1.0e-12),
    )

    assert isinstance(production, epcsaft.EquilibriumResult)
    assert production.diagnostics["reactive_workflow_class"] == "coupled_native"
    assert production.diagnostics["staged_route_used"] is False
    assert production.diagnostics["solver_backend"] == "ceres"
    assert production.diagnostics["jacobian_backend"] == "cppad_implicit"
    assert "staged_feed" not in production.diagnostics


def test_reactive_lle_dispatcher_does_not_accept_phase_route_controls() -> None:
    mix, feed, _initial_phases, reaction = _toy_reactive_phase_case()

    with np.testing.assert_raises_regex(epcsaft.InputError, "does not support phase_kwargs"):
        mix.equilibrium(
            kind="reactive_lle",
            T=298.15,
            P=1.0e5,
            z=feed,
            balances={"total": {"A": 1.0, "B": 1.0}},
            totals={"total": 1.0},
            reactions=[reaction],
            parent_phase="liq",
        )
