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


def test_solve_reactive_speciation_activity_coupled_state_uses_epcsaft_activities() -> None:
    species = ["H2O", "NaCl", "Na+", "Cl-"]
    mix = _salt_speciation_mixture()
    initial_x = np.asarray([0.998, 0.001, 0.0005, 0.0005], dtype=float)
    state = mix.state(T=298.15, P=1.0e5, x=initial_x, phase="liq")
    gamma = state.activity_coefficient(species=species)
    log_k = math.log(initial_x[2] * gamma["Na+"]) + math.log(initial_x[3] * gamma["Cl-"])
    log_k -= math.log(initial_x[1] * gamma["NaCl"])

    result = epcsaft.solve_reactive_speciation(
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
                log_equilibrium_constant=log_k,
                name="salt_dissociation",
            )
        ],
        initial_x=initial_x,
        options=epcsaft.ReactiveSpeciationOptions(max_iterations=8, tolerance=1.0e-8),
    )

    assert result.success is True
    assert result.diagnostics["activity_model"] == "epcsaft_component_activity"
    assert result.diagnostics["activity_fixed_point"] is True
    assert result.diagnostics["finite_difference_backend_available"] is False
    assert result.diagnostics["activity_derivative_in_jacobian"] is False
    assert result.named_reaction_residuals["salt_dissociation"] == pytest.approx(0.0, abs=1.0e-8)


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


def test_solve_reactive_speciation_concentration_standard_state_solves_with_density_activity() -> None:
    species = ["H2O", "NaCl", "Na+", "Cl-"]
    mix = _salt_speciation_mixture()
    initial_x = np.asarray([0.998, 0.001, 0.0005, 0.0005], dtype=float)
    state = mix.state(T=298.15, P=1.0e5, x=initial_x, phase="liq")
    density = state.molar_density()
    log_k = math.log(density * initial_x[2]) + math.log(density * initial_x[3])
    log_k -= math.log(density * initial_x[1])

    result = epcsaft.solve_reactive_speciation(
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
                log_equilibrium_constant=log_k,
                name="salt_dissociation",
                standard_state="concentration",
            )
        ],
        initial_x=initial_x,
        options=epcsaft.ReactiveSpeciationOptions(max_iterations=50, tolerance=1.0e-8),
    )

    assert result.success is True
    assert result.diagnostics["activity_model"] == "concentration"
    assert result.diagnostics["density_solve_count"] > 0
    assert result.named_reaction_residuals["salt_dissociation"] == pytest.approx(0.0, abs=1.0e-8)


def test_reactive_speciation_auto_jacobian_solves_concentration_standard_state() -> None:
    species = ["H2O", "NaCl", "Na+", "Cl-"]
    mix = _salt_speciation_mixture()
    initial_x = np.asarray([0.998, 0.001, 0.0005, 0.0005], dtype=float)
    density = mix.state(T=298.15, P=1.0e5, x=initial_x, phase="liq").molar_density()
    log_k = math.log(density * initial_x[2]) + math.log(density * initial_x[3])
    log_k -= math.log(density * initial_x[1])

    result = epcsaft.solve_reactive_speciation(
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
                log_equilibrium_constant=log_k,
                name="salt_dissociation",
                standard_state="concentration",
            )
        ],
        initial_x=initial_x,
        options=epcsaft.ReactiveSpeciationOptions(max_iterations=50, tolerance=1.0e-8),
    )

    assert result.success is True
    assert result.diagnostics["jacobian_backend"] == "analytic"
    assert result.diagnostics["activity_derivative_policy"] == "not_used_by_fixed_point_outer_iteration"


def test_concentration_standard_state_can_skip_activity_output() -> None:
    species = ["H2O", "NaCl", "Na+", "Cl-"]
    mix = _salt_speciation_mixture()
    initial_x = np.asarray([0.998, 0.001, 0.0005, 0.0005], dtype=float)
    density = mix.state(T=298.15, P=1.0e5, x=initial_x, phase="liq").molar_density()
    log_k = math.log(density * initial_x[2]) + math.log(density * initial_x[3])
    log_k -= math.log(density * initial_x[1])

    result = epcsaft.solve_reactive_speciation(
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
                {"NaCl": -1.0, "Na+": 1.0, "Cl-": 1.0},
                log_equilibrium_constant=log_k,
                standard_state="concentration",
            )
        ],
        initial_x=initial_x,
        options=epcsaft.ReactiveSpeciationOptions(activity_output="never"),
    )

    assert result.success is True
    assert result.activity_coefficients == {}


def test_solve_reactive_speciation_strict_failure_reports_best_state() -> None:
    species = ["H2O", "NaCl", "Na+", "Cl-"]
    mix = _salt_speciation_mixture()

    with pytest.raises(epcsaft.SolutionError) as excinfo:
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
                    log_equilibrium_constant=100.0,
                    name="salt_dissociation",
                    standard_state="ideal_mole_fraction",
                )
            ],
            initial_x=[0.998, 0.001, 0.0005, 0.0005],
            options=epcsaft.ReactiveSpeciationOptions(
                max_iterations=0,
                tolerance=1.0e-12,
            ),
        )

    diagnostics = excinfo.value.diagnostics
    assert set(diagnostics["best_x"]) == set(species)
    assert diagnostics["best_activity_coefficients"] == {}
    assert diagnostics["named_reaction_residuals"].keys() == {"salt_dissociation"}
    assert diagnostics["native_success"] is False
    json.dumps(diagnostics, allow_nan=False)


def test_solve_reactive_speciation_best_effort_returns_nonconverged_result() -> None:
    species = ["H2O", "NaCl", "Na+", "Cl-"]
    mix = _salt_speciation_mixture()

    result = epcsaft.solve_reactive_speciation(
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
                log_equilibrium_constant=100.0,
                name="salt_dissociation",
                standard_state="ideal_mole_fraction",
            )
        ],
        initial_x=[0.998, 0.001, 0.0005, 0.0005],
        options=epcsaft.ReactiveSpeciationOptions(
            max_iterations=0,
            tolerance=1.0e-12,
            return_best_effort=True,
        ),
    )

    assert isinstance(result, epcsaft.ReactiveSpeciationResult)
    assert result.success is False
    assert set(result.x) == set(species)
    assert result.activity_coefficients == {}
    assert set(result.named_reaction_residuals) == {"salt_dissociation"}
    assert result.diagnostics["best_x"] == result.x
    assert result.diagnostics["best_activity_coefficients"] == result.activity_coefficients
    json.dumps(result.to_dict(), allow_nan=False)


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


def test_reactive_speciation_sweep_uses_continuation_and_keeps_shape() -> None:
    mix = epcsaft.ePCSAFTMixture.from_params(
        {
            "m": np.asarray([1.0, 1.0]),
            "s": np.asarray([3.0, 3.0]),
            "e": np.asarray([200.0, 200.0]),
        },
        species=["A", "B"],
    )

    results = epcsaft.solve_reactive_speciation_sweep(
        species=["A", "B"],
        mixture_factory=lambda x, T, P: mix,
        points=[
            {"T": 298.15, "P": 1.0e5, "totals": {"total": 1.0}, "initial_x": [0.5, 0.5]},
            {"T": 298.15, "P": 1.0e5, "totals": {"total": 1.0}, "initial_x": [0.9, 0.1]},
            {"T": 298.15, "P": 1.0e5, "totals": {"total": 1.0}, "initial_x": [0.9, 0.1]},
        ],
        balances={"total": {"A": 1.0, "B": 1.0}},
        reactions=[
            epcsaft.ReactionDefinition(
                {"A": -1.0, "B": 1.0},
                log_equilibrium_constant=math.log(3.0),
                standard_state="ideal_mole_fraction",
            )
        ],
        options=epcsaft.ReactiveSpeciationOptions(error_mode="result"),
        continuation="auto",
    )

    assert len(results) == 3
    assert all(isinstance(result, epcsaft.ReactiveSpeciationResult) for result in results)
    assert all("composition" in result.continuation_state for result in results)
    assert results[1].diagnostics["continuation_used"] is True
    assert results[1].diagnostics["initial_x_source"] == "previous_successful_result"


def test_reactive_speciation_sweep_returns_failed_result_shape() -> None:
    mix = epcsaft.ePCSAFTMixture.from_params(
        {
            "m": np.asarray([1.0, 1.0]),
            "s": np.asarray([3.0, 3.0]),
            "e": np.asarray([200.0, 200.0]),
        },
        species=["A", "B"],
    )

    results = epcsaft.solve_reactive_speciation_sweep(
        species=["A", "B"],
        mixture_factory=lambda x, T, P: mix,
        points=[
            {"T": 298.15, "P": 1.0e5, "totals": {"total": 1.0}, "initial_x": [0.5, 0.5]},
            {"T": 298.15, "P": 1.0e5, "totals": {"missing": 1.0}, "initial_x": [0.5, 0.5]},
        ],
        balances={"total": {"A": 1.0, "B": 1.0}},
        reactions=[
            epcsaft.ReactionDefinition(
                {"A": -1.0, "B": 1.0},
                log_equilibrium_constant=math.log(3.0),
                standard_state="ideal_mole_fraction",
            )
        ],
        options=epcsaft.ReactiveSpeciationOptions(error_mode="result"),
    )

    assert len(results) == 2
    assert results[0].success is True
    assert results[1].success is False
    assert "Missing total" in results[1].message
    assert results[1].diagnostics["structured_failure"] is True


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


@pytest.mark.parametrize(
    ("kwargs", "message"),
    [
        (
            {"reactions": [epcsaft.ReactionDefinition({"Missing": -1.0, "Na+": 1.0}, 0.0)]},
            "Unknown species 'Missing' in reaction stoichiometry",
        ),
        (
            {"totals": {"water_total": 0.998, "sodium_total": 0.0015}},
            "Missing total for balance 'chloride_total'",
        ),
        (
            {"initial_x": [0.998, 0.001, 0.0005]},
            "initial_x length must match species length",
        ),
    ],
)
def test_solve_reactive_speciation_rejects_invalid_chemistry_inputs(kwargs: dict[str, object], message: str) -> None:
    species = ["H2O", "NaCl", "Na+", "Cl-"]
    mix = _salt_speciation_mixture()
    request = {
        "species": species,
        "mixture_factory": lambda x, T, P: mix,
        "T": 298.15,
        "P": 1.0e5,
        "balances": {
            "water_total": {"H2O": 1.0},
            "sodium_total": {"NaCl": 1.0, "Na+": 1.0},
            "chloride_total": {"NaCl": 1.0, "Cl-": 1.0},
        },
        "totals": {"water_total": 0.998, "sodium_total": 0.0015, "chloride_total": 0.0015},
        "reactions": [
            epcsaft.ReactionDefinition(
                stoichiometry={"NaCl": -1.0, "Na+": 1.0, "Cl-": 1.0},
                log_equilibrium_constant=0.0,
            )
        ],
        "initial_x": [0.998, 0.001, 0.0005, 0.0005],
    }
    request.update(kwargs)

    with pytest.raises(epcsaft.InputError, match=message):
        epcsaft.solve_reactive_speciation(**request)
