from __future__ import annotations

import math
from dataclasses import fields

import numpy as np
import pytest

import epcsaft


def _salt_speciation_mixture() -> epcsaft.ePCSAFTMixture:
    params = {
        "m": np.asarray([1.2047, 1.0, 1.0, 1.0]),
        "s": np.asarray([2.7927, 3.0, 2.8232, 2.7560]),
        "e": np.asarray([353.95, 200.0, 230.0, 170.0]),
        "z": np.asarray([0.0, 0.0, 1.0, -1.0]),
        "dielc": np.asarray([78.09, 8.0, 8.0, 8.0]),
        "d_born": np.asarray([0.0, 0.0, 3.445, 4.1]),
        "MW": np.asarray([18.01528e-3, 58.44e-3, 22.989e-3, 35.45e-3]),
    }
    return epcsaft.ePCSAFTMixture.from_params(params, species=["H2O", "NaCl", "Na+", "Cl-"])


def _assert_reactive_speciation_native_derivative_route_required(
    excinfo: pytest.ExceptionInfo[epcsaft.InputError],
) -> None:
    message = str(excinfo.value)
    assert "Native Ipopt reactive speciation currently supports ideal_mole_fraction standard states" in message
    assert "activity and concentration routes require the EOS derivative NLP blocks" in message


@pytest.mark.parametrize("standard_state", ["ideal_mole_fraction", "concentration", "mole_fraction_activity"])
def test_reaction_definition_accepts_supported_standard_states(standard_state: str) -> None:
    reaction = epcsaft.ReactionDefinition(
        stoichiometry={"NaCl": -1.0, "Na+": 1.0, "Cl-": 1.0},
        log_equilibrium_constant=0.0,
        standard_state=standard_state,
    )

    assert reaction.standard_state == standard_state

def test_reaction_definition_rejects_unknown_standard_state() -> None:
    with pytest.raises(epcsaft.InputError, match=r"ReactionDefinition\.standard_state"):
        epcsaft.ReactionDefinition(
            stoichiometry={"NaCl": -1.0, "Na+": 1.0, "Cl-": 1.0},
            log_equilibrium_constant=0.0,
            standard_state="unknown_basis",
        )

def test_reactive_speciation_options_public_surface_is_current_fields() -> None:
    assert {field.name for field in fields(epcsaft.ReactiveSpeciationOptions)} == {
        "max_iterations",
        "tolerance",
        "min_mole_fraction",
        "jacobian_backend",
        "solver_backend",
        "phase",
        "error_mode",
        "activity_output",
        "mass_tolerance",
        "charge_tolerance",
        "reaction_tolerance",
    }

@pytest.mark.parametrize(
    ("options", "message"),
    [
        (epcsaft.ReactiveSpeciationOptions(solver_backend="python_ipopt"), "solver_backend"),
        (epcsaft.ReactiveSpeciationOptions(jacobian_backend="autodiff"), "jacobian_backend"),
    ],
)
def test_reactive_speciation_rejects_invalid_optimizer_options(options, message) -> None:
    species = ["H2O", "NaCl", "Na+", "Cl-"]
    mix = _salt_speciation_mixture()

    with pytest.raises(epcsaft.InputError, match=message):
        epcsaft.solve_reactive_speciation(
            species=species,
            mixture_factory=lambda x, T, P: mix,
            T=298.15,
            P=1.0e5,
            balances={
                "water_total": {"H2O": 1.0},
                "sodium_total": {"NaCl": 1.0, "Na+": 1.0},
                "chloride_total": {"NaCl": 1.0, "Cl-": 1.0},
            },
            totals={"water_total": 0.998, "sodium_total": 0.0015, "chloride_total": 0.0015},
            reactions=[
                epcsaft.ReactionDefinition(
                    stoichiometry={"NaCl": -1.0, "Na+": 1.0, "Cl-": 1.0},
                    log_equilibrium_constant=0.0,
                )
            ],
            initial_x=[0.998, 0.001, 0.0005, 0.0005],
            options=options,
        )

def test_reactive_speciation_requested_ipopt_routes_ideal_speciation_when_compiled() -> None:
    from epcsaft import _core

    mix = epcsaft.ePCSAFTMixture.from_params(
        {
            "m": np.asarray([1.0, 1.0]),
            "s": np.asarray([3.0, 3.0]),
            "e": np.asarray([200.0, 200.0]),
        },
        species=["A", "B"],
    )
    kwargs = {
        "species": ["A", "B"],
        "mixture_factory": lambda x, T, P: mix,
        "T": 298.15,
        "P": 1.0e5,
        "balances": {"total": {"A": 1.0, "B": 1.0}},
        "totals": {"total": 1.0},
        "reactions": [
            epcsaft.ReactionDefinition(
                {"A": -1.0, "B": 1.0},
                log_equilibrium_constant=math.log(3.0),
                standard_state="ideal_mole_fraction",
            )
        ],
        "initial_x": [0.5, 0.5],
        "options": epcsaft.ReactiveSpeciationOptions(solver_backend="ipopt", tolerance=1.0e-9),
    }

    if not _core._native_ipopt_smoke()["compiled"]:
        with pytest.raises(epcsaft.SolutionError, match=r"EPCSAFT_ENABLE_IPOPT=ON"):
            epcsaft.solve_reactive_speciation(**kwargs)
        return

    result = epcsaft.solve_reactive_speciation(**kwargs)

    assert result.success is True
    assert result.x["B"] / result.x["A"] == pytest.approx(3.0, rel=1.0e-7)
    assert result.diagnostics["selected_solver_backend"] == "native_ipopt"
    assert result.diagnostics["problem_class"] == "homogeneous_ideal_gibbs_speciation"
    assert result.diagnostics["jacobian_backend"] == "analytic"

def test_reactive_speciation_requested_ipopt_handles_charged_ideal_constraint_when_compiled() -> None:
    from epcsaft import _core

    mix = epcsaft.ePCSAFTMixture.from_params(
        {
            "m": np.asarray([1.0, 1.0, 1.0]),
            "s": np.asarray([3.0, 3.0, 3.0]),
            "e": np.asarray([200.0, 200.0, 200.0]),
            "z": np.asarray([1.0, -1.0, 0.0]),
            "dielc": np.asarray([78.0, 78.0, 78.0]),
            "d_born": np.asarray([3.0, 3.0, 0.0]),
            "MW": np.asarray([20.0e-3, 20.0e-3, 40.0e-3]),
        },
        species=["C+", "A-", "N"],
    )
    kwargs = {
        "species": ["C+", "A-", "N"],
        "mixture_factory": lambda x, T, P: mix,
        "T": 298.15,
        "P": 1.0e5,
        "balances": {"formula_units": {"C+": 0.5, "A-": 0.5, "N": 1.0}},
        "totals": {"formula_units": 0.75},
        "reactions": [
            epcsaft.ReactionDefinition(
                {"C+": -1.0, "A-": -1.0, "N": 1.0},
                log_equilibrium_constant=math.log(8.0),
                standard_state="ideal_mole_fraction",
            )
        ],
        "initial_x": [0.3, 0.3, 0.4],
        "options": epcsaft.ReactiveSpeciationOptions(solver_backend="ipopt", tolerance=1.0e-9),
    }

    if not _core._native_ipopt_smoke()["compiled"]:
        with pytest.raises(epcsaft.SolutionError, match=r"EPCSAFT_ENABLE_IPOPT=ON"):
            epcsaft.solve_reactive_speciation(**kwargs)
        return

    result = epcsaft.solve_reactive_speciation(**kwargs)

    assert result.success is True
    assert result.x["C+"] == pytest.approx(0.25, rel=1.0e-7)
    assert result.x["A-"] == pytest.approx(0.25, rel=1.0e-7)
    assert result.x["N"] == pytest.approx(0.5, rel=1.0e-7)
    assert result.charge_residual == pytest.approx(0.0, abs=1.0e-10)
    assert result.diagnostics["charge_constraint_in_nlp"] is True

def test_reactive_speciation_auto_routes_ideal_speciation_to_native_ipopt_when_compiled() -> None:
    from epcsaft import _core

    mix = epcsaft.ePCSAFTMixture.from_params(
        {
            "m": np.asarray([1.0, 1.0]),
            "s": np.asarray([3.0, 3.0]),
            "e": np.asarray([200.0, 200.0]),
        },
        species=["A", "B"],
    )
    kwargs = {
        "species": ["A", "B"],
        "mixture_factory": lambda x, T, P: mix,
        "T": 298.15,
        "P": 1.0e5,
        "balances": {"total": {"A": 1.0, "B": 1.0}},
        "totals": {"total": 1.0},
        "reactions": [
            epcsaft.ReactionDefinition(
                {"A": -1.0, "B": 1.0},
                log_equilibrium_constant=math.log(3.0),
                standard_state="ideal_mole_fraction",
            )
        ],
        "initial_x": [0.5, 0.5],
    }

    if not _core._native_ipopt_smoke()["compiled"]:
        with pytest.raises(epcsaft.SolutionError, match=r"EPCSAFT_ENABLE_IPOPT=ON"):
            epcsaft.solve_reactive_speciation(**kwargs)
        return

    result = epcsaft.solve_reactive_speciation(**kwargs)

    assert result.success is True
    assert result.diagnostics["requested_solver_backend"] == "auto"
    assert result.diagnostics["selected_solver_backend"] == "native_ipopt"
    assert result.diagnostics["solver_selection_reason"] == "auto_selected_native_ipopt"


def test_reactive_speciation_auto_routes_activity_coupled_state_to_native_derivative_gate() -> None:
    mix = _salt_speciation_mixture()

    with pytest.raises(epcsaft.InputError) as excinfo:
        epcsaft.solve_reactive_speciation(
            species=["H2O", "NaCl", "Na+", "Cl-"],
            mixture_factory=lambda x, T, P: mix,
            T=298.15,
            P=1.0e5,
            balances={
                "water_total": {"H2O": 1.0},
                "sodium_total": {"NaCl": 1.0, "Na+": 1.0},
                "chloride_total": {"NaCl": 1.0, "Cl-": 1.0},
            },
            totals={"water_total": 0.998, "sodium_total": 0.0015, "chloride_total": 0.0015},
            reactions=[
                epcsaft.ReactionDefinition(
                    stoichiometry={"NaCl": -1.0, "Na+": 1.0, "Cl-": 1.0},
                    log_equilibrium_constant=0.0,
                    standard_state="mole_fraction_activity",
                )
            ],
            initial_x=[0.998, 0.001, 0.0005, 0.0005],
        )

    _assert_reactive_speciation_native_derivative_route_required(excinfo)
