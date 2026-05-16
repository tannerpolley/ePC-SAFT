from __future__ import annotations

import math

import numpy as np
import pytest

import epcsaft
from tests.equilibrium.core.test_vle import _assert_tp_flash_route_pending


def _toy_mixture() -> epcsaft.ePCSAFTMixture:
    return epcsaft.ePCSAFTMixture.from_params(
        {
            "m": np.asarray([1.0, 1.0]),
            "s": np.asarray([3.0, 3.0]),
            "e": np.asarray([200.0, 200.0]),
        },
        species=["A", "B"],
    )


def _fixed_literature_reaction() -> epcsaft.ReactionDefinition:
    return epcsaft.ReactionDefinition.from_literature_constant(
        {"A": -1.0, "B": 1.0},
        log_equilibrium_constant=math.log(3.0),
        name="literature_a_to_b",
        standard_state="ideal_mole_fraction",
        source="example literature table",
    )


def _successful_speciation_result() -> epcsaft.ReactiveSpeciationResult:
    return epcsaft.ReactiveSpeciationResult(
        success=True,
        message="converged",
        x={"A": 0.25, "B": 0.75},
        activity_coefficients={},
        mass_balance_residuals={"total": 0.0},
        charge_residual=0.0,
        reaction_residuals=[0.0],
        named_reaction_residuals={"literature_a_to_b": 0.0},
        state_failure_count=0,
        diagnostics={"phase_equilibrium_handoff": {}},
    )


def test_literature_reaction_constant_is_explicit_fixed_input() -> None:
    reaction = _fixed_literature_reaction()

    assert reaction.log_equilibrium_constant == pytest.approx(math.log(3.0))
    assert reaction.metadata["constant_source"] == "literature"
    assert reaction.metadata["source"] == "example literature table"
    assert reaction.metadata["fitting_role"] == "fixed_input"


def test_staged_workflow_requires_native_ipopt_phase_route_after_fixed_constant_speciation(monkeypatch) -> None:
    mix = _toy_mixture()
    monkeypatch.setattr(
        "epcsaft.reactive_staged.solve_reactive_speciation",
        lambda **kwargs: _successful_speciation_result(),
    )

    with pytest.raises(epcsaft.InputError) as excinfo:
        epcsaft.solve_reactive_staged_equilibrium(
            species=mix.species,
            mixture_factory=lambda x, T, P: mix,
            T=298.15,
            P=1.0e5,
            balances={"total": {"A": 1.0, "B": 1.0}},
            totals={"total": 1.0},
            reactions=[_fixed_literature_reaction()],
            initial_x=[0.5, 0.5],
            phase_kind="tp_flash",
        )

    _assert_tp_flash_route_pending(excinfo)


def test_reactive_phase_equilibrium_problem_is_public_generic_contract() -> None:
    assert "ReactivePhaseEquilibriumProblem" in epcsaft.__all__


def test_staged_workflow_rejects_reaction_constant_fit_as_default_role() -> None:
    mix = _toy_mixture()

    with pytest.raises(epcsaft.InputError, match="reaction-constant fitting is not a default"):
        epcsaft.solve_reactive_staged_equilibrium(
            species=mix.species,
            mixture_factory=lambda x, T, P: mix,
            T=298.15,
            P=1.0e5,
            balances={"total": {"A": 1.0, "B": 1.0}},
            totals={"total": 1.0},
            reactions=[_fixed_literature_reaction()],
            initial_x=[0.5, 0.5],
            phase_kind="tp_flash",
            workflow_options={"reaction_constant_fitting": "primary"},
        )
