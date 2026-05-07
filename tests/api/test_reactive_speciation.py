from __future__ import annotations

import json
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


def test_solve_reactive_speciation_returns_balanced_activity_coupled_state() -> None:
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
        options=epcsaft.ReactiveSpeciationOptions(
            max_iterations=8,
            tolerance=1.0e-8,
            mass_tolerance=1.0e-8,
            charge_tolerance=1.0e-8,
            reaction_tolerance=1.0e-8,
        ),
    )

    assert isinstance(result, epcsaft.ReactiveSpeciationResult)
    assert result.success is True
    assert sum(result.x.values()) == pytest.approx(1.0, abs=1.0e-10)
    assert abs(result.charge_residual) <= 1.0e-10
    assert max(abs(value) for value in result.mass_balance_residuals.values()) <= 1.0e-8
    assert max(abs(value) for value in result.reaction_residuals) <= 1.0e-8
    assert set(result.named_reaction_residuals) == {"salt_dissociation"}
    assert result.named_reaction_residuals["salt_dissociation"] == pytest.approx(result.reaction_residuals[0])
    assert set(result.activity_coefficients) == set(species)
    assert result.state_failure_count == 0
    assert result.diagnostics["solver_language"] == "c++"
    assert result.diagnostics["backend"] == "native"
    assert result.diagnostics["native_entrypoint"] == "_solve_chemical_equilibrium_native"
    assert result.diagnostics["mass_tolerance"] == pytest.approx(1.0e-8)
    assert result.diagnostics["charge_tolerance"] == pytest.approx(1.0e-8)
    assert result.diagnostics["reaction_tolerance"] == pytest.approx(1.0e-8)
    json.dumps(result.to_dict(), allow_nan=False)


@pytest.mark.parametrize("standard_state", ["ideal_mole_fraction", "concentration", "mole_fraction_activity"])
def test_reaction_definition_accepts_supported_standard_states(standard_state: str) -> None:
    reaction = epcsaft.ReactionDefinition(
        stoichiometry={"NaCl": -1.0, "Na+": 1.0, "Cl-": 1.0},
        log_equilibrium_constant=0.0,
        standard_state=standard_state,
    )

    assert reaction.standard_state == standard_state


def test_reaction_definition_rejects_unknown_standard_state() -> None:
    with pytest.raises(epcsaft.InputError, match="ReactionDefinition.standard_state"):
        epcsaft.ReactionDefinition(
            stoichiometry={"NaCl": -1.0, "Na+": 1.0, "Cl-": 1.0},
            log_equilibrium_constant=0.0,
            standard_state="molality",
        )


def test_solve_reactive_speciation_concentration_standard_state_uses_molar_density() -> None:
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
    assert result.reaction_residuals == pytest.approx([0.0], abs=1.0e-8)
    assert result.diagnostics["reaction_standard_states"] == ["concentration"]
    assert result.diagnostics["activity_basis"] == "concentration"


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
                    log_equilibrium_constant=0.0,
                    name="salt_dissociation",
                )
            ],
            initial_x=[0.998, 0.001, 0.0005, 0.0005],
            options=epcsaft.ReactiveSpeciationOptions(max_iterations=0, tolerance=1.0e-12),
        )

    diagnostics = excinfo.value.diagnostics
    assert set(diagnostics["best_x"]) == set(species)
    assert set(diagnostics["best_activity_coefficients"]) == set(species)
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
                log_equilibrium_constant=0.0,
                name="salt_dissociation",
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
    assert set(result.activity_coefficients) == set(species)
    assert set(result.named_reaction_residuals) == {"salt_dissociation"}
    assert result.diagnostics["best_x"] == result.x
    assert result.diagnostics["best_activity_coefficients"] == result.activity_coefficients
    json.dumps(result.to_dict(), allow_nan=False)


def test_reactive_speciation_options_expose_jacobian_backend_selector() -> None:
    assert "backend" not in {field.name for field in fields(epcsaft.ReactiveSpeciationOptions)}
    assert "jacobian_backend" in {field.name for field in fields(epcsaft.ReactiveSpeciationOptions)}


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
