from __future__ import annotations

import math

import numpy as np
import pytest

import epcsaft


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


def test_literature_reaction_constant_is_explicit_fixed_input() -> None:
    reaction = _fixed_literature_reaction()

    assert reaction.log_equilibrium_constant == pytest.approx(math.log(3.0))
    assert reaction.metadata["constant_source"] == "literature"
    assert reaction.metadata["source"] == "example literature table"
    assert reaction.metadata["fitting_role"] == "fixed_input"


def test_staged_workflow_reports_fixed_constant_boundaries() -> None:
    mix = _toy_mixture()

    result = epcsaft.solve_reactive_staged_equilibrium(
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

    assert result.success is True
    assert result.diagnostics["reactive_workflow_class"] == "staged"
    assert result.diagnostics["reaction_constant_policy"] == "fixed_literature_constants_first"
    assert result.diagnostics["reaction_constant_fitting_role"] == "secondary_optional"
    assert (
        result.diagnostics["parameter_regression_boundary"] == "fit_epcsaft_parameters_after_fixed_constant_speciation"
    )
    assert result.diagnostics["full_simultaneous_reactive_nlp"] is False
    assert result.diagnostics["derivative_policy"]["finite_difference_backend_available"] is False
    assert result.diagnostics["derivative_policy"]["unsupported_derivative_status"] == "backend_unavailable"
    assert result.chemical.diagnostics["reaction_constant_sources"]["literature_a_to_b"] == "literature"


def test_staged_workflow_rejects_reaction_constant_fit_as_default_role() -> None:
    mix = _toy_mixture()

    with pytest.raises(epcsaft.InputError, match="reaction-constant fitting is not a staged default"):
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
