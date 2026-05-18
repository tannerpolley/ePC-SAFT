from __future__ import annotations

import math

import numpy as np
import pytest

import epcsaft
from tests.api.reactive.test_reactive_speciation_options import (
    _assert_reactive_speciation_native_ipopt_dependency_required,
    _native_ipopt_compiled,
)


def _toy_mixture() -> epcsaft.ePCSAFTMixture:
    return epcsaft.ePCSAFTMixture.from_params(
        {
            "m": np.asarray([1.0, 1.0], dtype=float),
            "s": np.asarray([3.0, 3.0], dtype=float),
            "e": np.asarray([200.0, 200.0], dtype=float),
        },
        species=["A", "B"],
    )


@pytest.mark.parametrize(
    ("standard_state", "basis"),
    [
        ("ideal_mole_fraction", "mole_fraction"),
        ("mole_fraction_activity", "mole_fraction"),
        ("thermodynamic_activity", "activity"),
        ("concentration", "concentration"),
        ("molality", "molality"),
        ("apparent", "apparent"),
    ],
)
def test_reaction_constant_convention_defines_supported_bases(standard_state: str, basis: str) -> None:
    convention = epcsaft.ReactionConstantConvention(standard_state=standard_state)

    assert convention.standard_state == standard_state
    assert convention.basis == basis
    assert convention.constant_kind == "thermodynamic"
    assert convention.fitting_role == "fixed_input"


def test_reaction_constant_convention_rejects_incompatible_basis() -> None:
    with pytest.raises(epcsaft.InputError, match="basis is incompatible"):
        epcsaft.ReactionConstantConvention(standard_state="molality", basis="mole_fraction")


def test_reaction_definition_carries_explicit_convention_metadata() -> None:
    convention = epcsaft.ReactionConstantConvention(
        standard_state="thermodynamic_activity",
        source="generic table",
    )
    reaction = epcsaft.ReactionDefinition.from_literature_constant(
        {"A": -1.0, "B": 1.0},
        log_equilibrium_constant=math.log(2.0),
        name="a_to_b",
        convention=convention,
    )

    assert reaction.standard_state == "thermodynamic_activity"
    assert reaction.convention.basis == "activity"
    assert reaction.metadata["constant_source"] == "literature"
    assert reaction.metadata["fitting_role"] == "fixed_input"
    assert reaction.metadata["constant_convention"]["standard_state"] == "thermodynamic_activity"
    assert reaction.metadata["constant_convention"]["native_standard_state_code"] == 0


def test_fitted_reaction_constant_must_be_explicit_not_default() -> None:
    reaction = epcsaft.ReactionDefinition.from_fitted_constant(
        {"A": -1.0, "B": 1.0},
        log_equilibrium_constant=math.log(2.0),
        name="fitted_a_to_b",
    )

    assert reaction.metadata["constant_source"] == "fit"
    assert reaction.convention.constant_kind == "fitted"
    assert reaction.convention.fitting_role == "fitted_parameter"


def test_corrected_convention_rejects_fixed_input_fitting_role() -> None:
    with pytest.raises(epcsaft.InputError, match="non-fixed fitting_role"):
        epcsaft.ReactionConstantConvention(
            constant_kind="corrected",
            correction_terms={"regularization": 1.0},
        )


def test_reactive_speciation_auto_validates_conventions_before_native_ipopt_route() -> None:
    mix = _toy_mixture()
    reaction = epcsaft.ReactionDefinition.from_literature_constant(
        {"A": -1.0, "B": 1.0},
        log_equilibrium_constant=math.log(3.0),
        name="a_to_b",
        standard_state="mole_fraction_activity",
        source="generic fixed table",
    )

    kwargs = {
        "species": mix.species,
        "mixture_factory": lambda x, T, P: mix,
        "T": 298.15,
        "P": 1.0e5,
        "balances": {"total": {"A": 1.0, "B": 1.0}},
        "totals": {"total": 1.0},
        "reactions": [reaction],
        "initial_x": [0.5, 0.5],
        "options": epcsaft.ReactiveSpeciationOptions(tolerance=1.0e-10),
    }

    if not _native_ipopt_compiled():
        with pytest.raises(epcsaft.SolutionError) as excinfo:
            epcsaft.solve_reactive_speciation(**kwargs)
        _assert_reactive_speciation_native_ipopt_dependency_required(excinfo)
    else:
        result = epcsaft.solve_reactive_speciation(**kwargs)
        assert result.success is True
        assert result.diagnostics["problem_class"] == "homogeneous_nonideal_residual_speciation"
        assert result.diagnostics["derivative_backend"] == "cppad_implicit"

    assert reaction.metadata["constant_convention"]["standard_state"] == "mole_fraction_activity"
    assert reaction.metadata["constant_convention"]["basis"] == "mole_fraction"
    assert reaction.metadata["constant_convention"]["constant_kind"] == "thermodynamic"
    assert reaction.metadata["constant_convention"]["fitting_role"] == "fixed_input"


def test_reactive_speciation_gates_apparent_convention_before_native_ipopt_route() -> None:
    mix = _toy_mixture()
    reaction = epcsaft.ReactionDefinition.from_fitted_constant(
        {"A": -1.0, "B": 1.0},
        log_equilibrium_constant=math.log(2.0),
        name="fitted_apparent_a_to_b",
        convention=epcsaft.ReactionConstantConvention(
            standard_state="apparent",
            constant_kind="fitted",
            fitting_role="fitted_parameter",
        ),
    )

    with pytest.raises(epcsaft.InputError, match=r"native speciation requires.*apparent"):
        epcsaft.solve_reactive_speciation(
            species=mix.species,
            mixture_factory=lambda x, T, P: mix,
            T=298.15,
            P=1.0e5,
            balances={"total": {"A": 1.0, "B": 1.0}},
            totals={"total": 1.0},
            reactions=[reaction],
            initial_x=[0.5, 0.5],
            options=epcsaft.ReactiveSpeciationOptions(solver_backend="ipopt", tolerance=1.0e-10),
        )

    assert reaction.metadata["constant_convention"]["standard_state"] == "apparent"
    assert reaction.metadata["constant_convention"]["constant_kind"] == "fitted"
    assert reaction.metadata["constant_convention"]["fitting_role"] == "fitted_parameter"


def test_molality_native_route_fails_at_speciation_contract() -> None:
    mix = _toy_mixture()
    reaction = epcsaft.ReactionDefinition(
        {"A": -1.0, "B": 1.0},
        log_equilibrium_constant=math.log(3.0),
        convention=epcsaft.ReactionConstantConvention(standard_state="molality"),
    )

    with pytest.raises(epcsaft.InputError, match=r"native speciation requires.*molality"):
        epcsaft.solve_reactive_speciation(
            species=mix.species,
            mixture_factory=lambda x, T, P: mix,
            T=298.15,
            P=1.0e5,
            balances={"total": {"A": 1.0, "B": 1.0}},
            totals={"total": 1.0},
            reactions=[reaction],
            initial_x=[0.5, 0.5],
        )
