"""Package-side generic contract smoke tests for downstream-facing APIs.

These tests prove that a local package install exposes generic public problem,
result, target, and capability contracts. They do not replace the real
downstream workflow runs required by issue #119.
"""

from __future__ import annotations

import math

import epcsaft


def _fixed_reaction(stoichiometry: dict[str, float], *, name: str) -> epcsaft.ReactionDefinition:
    return epcsaft.ReactionDefinition.from_literature_constant(
        stoichiometry,
        log_equilibrium_constant=math.log(2.0),
        name=name,
        standard_state="ideal_mole_fraction",
        source="downstream smoke fixture",
    )


def _generic_result_payload(problem_kind: str, species: tuple[str, ...]) -> dict[str, object]:
    phase = epcsaft.EquilibriumPhase(
        label="liq",
        composition=[1.0 / len(species)] * len(species),
        density=55000.0,
        temperature=298.15,
        pressure=101325.0,
        phase_fraction=1.0,
        ln_fugacity_coefficient=[0.0] * len(species),
        diagnostics={"species": list(species)},
    )
    return epcsaft.EquilibriumResult(
        backend="native",
        problem_kind=problem_kind,
        phases=(phase,),
        stable=True,
        split_detected=False,
        diagnostics={
            "derivative_policy": {
                "accepted_derivative_backends": ["analytic", "cppad", "analytic_implicit", "cppad_implicit"],
            }
        },
    ).to_dict()


def _generic_chemical_payload(species: tuple[str, ...]) -> dict[str, object]:
    x = {label: 1.0 / len(species) for label in species}
    return epcsaft.ReactiveSpeciationResult(
        success=True,
        message="converged",
        x=x,
        activity_coefficients={label: 1.0 for label in species},
        mass_balance_residuals={"total": 0.0},
        charge_residual=0.0,
        reaction_residuals=[0.0],
        named_reaction_residuals={"fixed_reaction": 0.0},
        state_failure_count=0,
        diagnostics={
            "backend": "native",
            "derivative_policy": {
                "accepted_derivative_backends": ["analytic", "cppad", "analytic_implicit", "cppad_implicit"],
            },
        },
    ).to_dict()


def test_public_api_stays_application_neutral_and_numerical_derivative_free() -> None:
    forbidden_tokens = {
        "absorption",
        "column",
        "distribution",
        "extraction",
        "lithium",
        "mea",
        "selectivity",
    }

    public_names = {name.lower() for name in epcsaft.__all__}
    assert not {name for name in public_names for token in forbidden_tokens if token in name}

    capabilities = epcsaft.capabilities()
    numerical_derivative = capabilities["derivatives"]["numerical_derivative"]
    assert numerical_derivative["available"] is False
    assert numerical_derivative["reason"] == "numerical_derivative_derivatives_forbidden"
    assert capabilities["equilibrium"]["problem_objects"]["entrypoint"] == "mixture.solve_equilibrium(problem)"


def test_mea_thermodynamics_smoke_uses_generic_reactive_problem_and_targets() -> None:
    species = ("CO2", "H2O", "Amine", "AmineH+", "HCO3-")
    problem = epcsaft.ReactiveSpeciationProblem(
        T=313.15,
        P=101325.0,
        balances={
            "carbon": {"CO2": 1.0, "HCO3-": 1.0},
            "amine": {"Amine": 1.0, "AmineH+": 1.0},
            "charge": {"AmineH+": 1.0, "HCO3-": -1.0},
        },
        totals={"carbon": 0.03, "amine": 0.12, "charge": 0.0},
        reactions=(
            _fixed_reaction(
                {"CO2": -1.0, "Amine": -1.0, "AmineH+": 1.0, "HCO3-": 1.0},
                name="fixed_reaction",
            ),
        ),
        initial_x=[0.02, 0.82, 0.10, 0.03, 0.03],
        options=epcsaft.ReactiveSpeciationOptions(jacobian_backend="cppad"),
    )
    dataset = epcsaft.TargetDataset.from_records(
        [
            {"row_family": "speciation", "T": 313.15, "target_x": {"CO2": 0.02, "HCO3-": 0.03}},
            {"row_family": "vle_partial_pressure", "T": 313.15, "target_partial_pressures": {"CO2": 2500.0}},
            {"row_family": "activity", "T": 313.15, "activity_H2O": 0.97},
        ],
        name="downstream-reactive-smoke",
        species=species,
    )

    output = _generic_chemical_payload(species)

    assert isinstance(problem, epcsaft.EquilibriumProblem)
    assert dataset.families == ("speciation", "vle_partial_pressure", "activity")
    assert output["x"]["CO2"] > 0.0


def test_lithium_extraction_smoke_uses_generic_electrolyte_lle_contracts() -> None:
    species = ("H2O", "OrganicSolvent", "Li+", "Cl-")
    problem = epcsaft.ElectrolyteLLEProblem(
        T=298.15,
        P=101325.0,
        solvent_feed={"H2O": 0.7, "OrganicSolvent": 0.3},
        salt_molality={"LiCl": 0.15},
        options=epcsaft.EquilibriumOptions(jacobian_backend="cppad"),
    )
    dataset = epcsaft.TargetDataset.from_records(
        [
            {
                "row_family": "lle_phase_composition",
                "T": 298.15,
                "P": 101325.0,
                "target_x_alpha": {"H2O": 0.93, "Li+": 0.01},
                "target_x_beta": {"OrganicSolvent": 0.88, "Li+": 0.02},
            },
            {"row_family": "mean_ionic_activity", "T": 298.15, "P": 101325.0, "mean_ionic_activity": 0.91},
            {"row_family": "regularization", "parameter": "k_ij", "target_value": 0.0},
        ],
        name="downstream-electrolyte-lle-smoke",
        species=species,
    )

    output = _generic_result_payload("electrolyte_lle", species)

    assert isinstance(problem, epcsaft.EquilibriumProblem)
    assert dataset.families == ("lle_phase_composition", "mean_ionic_activity", "regularization")
    assert output["problem_kind"] == "electrolyte_lle"
    assert output["phases"][0]["composition"]


def test_absorption_column_smoke_uses_generic_reactive_bubble_contracts() -> None:
    species = ("CO2", "H2O", "Amine", "AmineH+", "HCO3-")
    problem = epcsaft.ReactiveElectrolyteBubbleProblem(
        T=313.15,
        P=101325.0,
        balances={
            "carbon": {"CO2": 1.0, "HCO3-": 1.0},
            "amine": {"Amine": 1.0, "AmineH+": 1.0},
            "charge": {"AmineH+": 1.0, "HCO3-": -1.0},
        },
        totals={"carbon": 0.04, "amine": 0.15, "charge": 0.0},
        reactions=(
            _fixed_reaction(
                {"CO2": -1.0, "Amine": -1.0, "AmineH+": 1.0, "HCO3-": 1.0},
                name="fixed_reaction",
            ),
        ),
        initial_x=[0.03, 0.78, 0.12, 0.035, 0.035],
        vapor_species=("CO2", "H2O"),
        options=epcsaft.ReactiveElectrolyteBubbleOptions(error_mode="result"),
    )
    dataset = epcsaft.TargetDataset.from_records(
        [
            {"row_family": "vle_partial_pressure", "T": 313.15, "target_partial_pressures": {"CO2": 1800.0}},
            {"row_family": "speciation", "T": 313.15, "target_x": {"CO2": 0.03, "AmineH+": 0.035}},
            {"row_family": "fugacity", "T": 313.15, "P": 101325.0, "fugacity_CO2": 1800.0},
        ],
        name="downstream-reactive-bubble-smoke",
        species=species,
    )

    output = _generic_result_payload("reactive_electrolyte_bubble_pressure", species)

    assert isinstance(problem, epcsaft.EquilibriumProblem)
    assert dataset.families == ("vle_partial_pressure", "speciation", "fugacity")
    assert output["problem_kind"] == "reactive_electrolyte_bubble_pressure"
