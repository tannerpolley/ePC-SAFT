from __future__ import annotations

import json
import math
from dataclasses import fields

import numpy as np
import pytest

import epcsaft
import epcsaft.ipopt_backend as ipopt_backend

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

def test_reactive_speciation_options_expose_jacobian_backend_selector() -> None:
    assert "backend" not in {field.name for field in fields(epcsaft.ReactiveSpeciationOptions)}
    assert "jacobian_backend" in {field.name for field in fields(epcsaft.ReactiveSpeciationOptions)}
    assert "finite" + "_difference_step" not in {field.name for field in fields(epcsaft.ReactiveSpeciationOptions)}
    assert "solver_backend" in {field.name for field in fields(epcsaft.ReactiveSpeciationOptions)}
    assert "hessian_strategy" in {field.name for field in fields(epcsaft.ReactiveSpeciationOptions)}

@pytest.mark.parametrize(
    ("options", "message"),
    [
        (epcsaft.ReactiveSpeciationOptions(solver_backend="cyipopt"), "solver_backend"),
        (epcsaft.ReactiveSpeciationOptions(hessian_strategy="exact"), "hessian_strategy"),
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

def test_reactive_speciation_requested_ipopt_requires_cyipopt(monkeypatch) -> None:
    monkeypatch.setattr(ipopt_backend, "cyipopt_available", lambda: False)
    species = ["H2O", "NaCl", "Na+", "Cl-"]
    mix = _salt_speciation_mixture()

    with pytest.raises(epcsaft.InputError, match=r"cyipopt.*solver_backend='ipopt'"):
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
            options=epcsaft.ReactiveSpeciationOptions(solver_backend="ipopt"),
        )

def test_reactive_speciation_auto_does_not_require_cyipopt(monkeypatch) -> None:
    monkeypatch.setattr(ipopt_backend, "require_cyipopt", lambda route: (_ for _ in ()).throw(AssertionError(route)))
    mix = epcsaft.ePCSAFTMixture.from_params(
        {
            "m": np.asarray([1.0, 1.0]),
            "s": np.asarray([3.0, 3.0]),
            "e": np.asarray([200.0, 200.0]),
        },
        species=["A", "B"],
    )

    result = epcsaft.solve_reactive_speciation(
        species=["A", "B"],
        mixture_factory=lambda x, T, P: mix,
        T=298.15,
        P=1.0e5,
        balances={"total": {"A": 1.0, "B": 1.0}},
        totals={"total": 1.0},
        reactions=[
            epcsaft.ReactionDefinition(
                {"A": -1.0, "B": 1.0},
                log_equilibrium_constant=math.log(3.0),
                standard_state="ideal_mole_fraction",
            )
        ],
        initial_x=[0.5, 0.5],
    )

    assert result.success is True
    assert result.diagnostics["selected_solver_backend"] == "native"
    assert result.diagnostics["solver_selection_reason"] == "default_native"

@pytest.mark.skipif(not ipopt_backend.cyipopt_available(), reason="cyipopt is optional")
def test_reactive_speciation_ipopt_solves_easy_ideal_reaction() -> None:
    mix = epcsaft.ePCSAFTMixture.from_params(
        {
            "m": np.asarray([1.0, 1.0]),
            "s": np.asarray([3.0, 3.0]),
            "e": np.asarray([200.0, 200.0]),
        },
        species=["A", "B"],
    )

    result = epcsaft.solve_reactive_speciation(
        species=["A", "B"],
        mixture_factory=lambda x, T, P: mix,
        T=298.15,
        P=1.0e5,
        balances={"total": {"A": 1.0, "B": 1.0}},
        totals={"total": 1.0},
        reactions=[epcsaft.ReactionDefinition({"A": -1.0, "B": 1.0}, log_equilibrium_constant=math.log(3.0))],
        initial_x=[0.5, 0.5],
        options=epcsaft.ReactiveSpeciationOptions(
            solver_backend="ipopt",
            hessian_strategy="lbfgs",
            tolerance=1.0e-8,
            max_iterations=80,
        ),
    )

    assert result.success is True
    assert result.x["B"] / result.x["A"] == pytest.approx(3.0, rel=1.0e-7)
    assert result.diagnostics["backend"] == "ipopt"
    assert result.diagnostics["solver_method"] == "cyipopt_bound_min_residual"
    assert result.diagnostics["formulation"] == "bound_constrained_residual_minimization"
    assert result.diagnostics["full_constrained_nlp"] is False
    assert result.diagnostics["exact_hessian_available"] is False
    assert result.diagnostics["hessian_strategy"] == "lbfgs"
    assert result.diagnostics["hessian_kind"] == "ipopt_limited_memory"
    assert result.diagnostics["hessian_includes_second_residual_derivatives"] is False
    assert result.diagnostics["ipopt_success"] is True
    assert result.diagnostics["residual_gate_success"] is True
    assert result.diagnostics["accepted"] is True
    assert result.diagnostics["selected_solver_backend"] == "ipopt"
    assert result.diagnostics["solver_selection_reason"] == "explicit_request"
