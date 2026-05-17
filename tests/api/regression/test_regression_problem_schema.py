"""API contract tests for user-facing regression problem schemas."""

from __future__ import annotations

import pytest

import epcsaft
import epcsaft.regression as regression_module
from epcsaft import FitProblem, FitResult, TargetDataset, TargetRow
from epcsaft._types import InputError


def test_fit_pure_parameters_accepts_minimal_user_contract(monkeypatch):
    calls = []

    def fake_fit_pure_neutral(records, component, **kwargs):
        calls.append((records, component, kwargs))
        return FitResult(
            problem=FitProblem(
                mode="pure_neutral",
                records=tuple(records),
                component=component,
                fit_targets=tuple(kwargs["fit_targets"]),
                fixed_parameters=kwargs["fixed_parameters"],
                bounds={"m": (0.8, 1.4)},
                weights={"density": 2.0},
                loss="linear",
                solver_options={"max_nfev": 4},
                output_report=True,
            ),
            fitted_values={"m": 1.1, "s": 3.7, "e": 145.0},
            success=True,
        )

    monkeypatch.setattr(regression_module, "fit_pure_neutral", fake_fit_pure_neutral)

    result = epcsaft.fit_pure_parameters(
        species="Methane",
        data_rows=[{"T": 100.0, "P": 34375.0, "rho_sat_liq_kg_m3": 438.0}],
        parameters_to_fit=("m", "sigma", "epsilon"),
        fixed_parameters={"MW": 0.0160428, "z": 0.0},
        bounds={"m": (0.8, 1.4)},
        weights={"density": 2.0},
        loss="linear",
        solver_options={"max_nfev": 4},
        output_report=True,
    )

    assert result.success is True
    assert result.problem.fit_targets == ("m", "s", "e")
    assert result.problem.bounds == {"m": (0.8, 1.4)}
    assert result.problem.weights == {"density": 2.0}
    assert result.problem.loss == "linear"
    assert result.problem.solver_options == {"max_nfev": 4}
    assert result.problem.output_report is True
    assert calls[0][1] == "Methane"
    assert calls[0][2]["fit_targets"] == ("m", "s", "e")


def test_fit_binary_parameters_accepts_minimal_user_contract(monkeypatch):
    calls = []

    def fake_fit_binary_pair(records, pair, **kwargs):
        calls.append((records, pair, kwargs))
        return FitResult(
            problem=FitProblem(
                mode="binary_pair",
                records=tuple(records),
                pair=tuple(pair),
                dataset=str(kwargs["dataset"]),
                fit_targets=tuple(kwargs["fit_targets"]),
                bounds={"k_ij": (-0.2, 0.2)},
                weights={"binary_vle": 1.5},
                loss="soft_l1",
                solver_options={"max_nfev": 3},
            ),
            fitted_values={"k_ij": -0.01},
            success=True,
        )

    monkeypatch.setattr(regression_module, "fit_binary_pair", fake_fit_binary_pair)

    result = epcsaft.fit_binary_parameters(
        species=("H2O", "Ethanol"),
        data_rows=[{"T": 330.0, "P": 101325.0, "x_H2O": 0.7, "x_Ethanol": 0.3, "y_H2O": 0.5, "y_Ethanol": 0.5}],
        parameter_set="2026_Khudaida",
        parameters_to_fit=("k_ij",),
        bounds={"k_ij": (-0.2, 0.2)},
        weights={"binary_vle": 1.5},
        loss="soft_l1",
        solver_options={"max_nfev": 3},
    )

    assert result.success is True
    assert result.problem.mode == "binary_pair"
    assert result.problem.pair == ("H2O", "Ethanol")
    assert result.problem.fit_targets == ("k_ij",)
    assert result.problem.loss == "soft_l1"
    assert calls[0][2]["dataset"] == "2026_Khudaida"


def test_fit_liquid_electrolyte_parameters_uses_native_ceres_contract():
    result = epcsaft.fit_liquid_electrolyte_parameters(
        species=("H2O", "Na+", "Cl-"),
        data_rows=[
            {
                "T": 298.15,
                "P": 101325.0,
                "x_H2O": 0.98,
                "x_Na+": 0.01,
                "x_Cl-": 0.01,
                "osmotic_coefficient": 1.0,
            }
        ],
        parameter_set="2026_Khudaida",
        parameters_to_fit=("d_born", "f_solv"),
        fixed_parameters={"H2O": {"m": 1.0}},
        weights={"osmotic_coefficient": 1.0},
        loss="linear",
        solver_options={"max_nfev": 1},
    )

    assert result.success is True
    assert result.message.startswith("Ceres Solver Report")
    assert "without" + " optimizer" not in result.message
    assert result.backend == "ceres"
    assert result.optimizer_backend == "ceres"
    assert result.derivative_backend == "cppad_implicit"
    assert result.jacobian_backend == "cppad_implicit"
    assert result.python_objective_used is False
    assert result.problem.mode == "liquid_electrolyte"
    assert result.problem.fit_targets == ("d_born", "f_solv")
    assert result.problem.dataset == "2026_Khudaida"
    assert result.problem.weights == {"osmotic_coefficient": 1.0}
    assert {row["row_family"] for row in result.row_diagnostics} == {"osmotic_coefficient"}


def test_target_dataset_accepts_generic_row_families_and_round_trips_records():
    dataset = TargetDataset.from_records(
        [
            {"row_family": "pure_density", "row_id": "rho-1", "T": 298.15, "P": 101325.0, "rho_kg_m3": 997.0},
            {"row_family": "pure_vapor_pressure", "T": 373.15, "vapor_pressure": 101325.0},
            {
                "row_family": "binary_vle",
                "T": 330.0,
                "P": 101325.0,
                "x_H2O": 0.7,
                "x_Ethanol": 0.3,
                "y_H2O": 0.5,
                "y_Ethanol": 0.5,
            },
            {
                "row_family": "binary_lle",
                "T": 298.15,
                "P": 101325.0,
                "x_alpha_H2O": 0.8,
                "x_alpha_Ethanol": 0.2,
                "x_beta_H2O": 0.1,
                "x_beta_Ethanol": 0.9,
            },
            {"row_family": "osmotic_coefficient", "T": 298.15, "P": 101325.0, "osmotic_coefficient": 0.93},
            {"row_family": "MIAC", "T": 298.15, "mean_ionic_activity": 0.91},
            {"row_family": "relative_permittivity", "T": 298.15, "epsilon_r_exp": 78.3},
            {"row_family": "activity", "T": 298.15, "activity_H2O": 0.98},
            {"row_family": "fugacity", "T": 298.15, "fugacity_CO2": 1200.0},
            {"row_family": "speciation", "T": 298.15, "target_x": {"CO2": 0.01}},
            {"row_family": "vle_partial_pressure", "T": 313.15, "target_partial_pressures": {"CO2": 25000.0}},
            {"row_family": "lle_phase_composition", "T": 298.15, "target_x_alpha": {"H2O": 0.8}},
            {"row_family": "regularization", "parameter": "k_ij", "target_value": 0.0},
        ],
        name="generic-schema-smoke",
        species=("H2O", "Ethanol", "CO2"),
    )

    assert dataset.name == "generic-schema-smoke"
    assert dataset.species == ("H2O", "Ethanol", "CO2")
    assert dataset.families == (
        "pure_density",
        "pure_vapor_pressure",
        "binary_vle",
        "binary_lle",
        "osmotic_coefficient",
        "mean_ionic_activity",
        "relative_permittivity",
        "activity",
        "fugacity",
        "speciation",
        "vle_partial_pressure",
        "lle_phase_composition",
        "regularization",
    )
    assert dataset.rows[5].row_family == "mean_ionic_activity"
    assert dataset.to_records()[0]["row_id"] == "rho-1"
    assert dataset.to_records()[0]["rho_kg_m3"] == 997.0


def test_target_row_schema_rejects_unknown_family_missing_targets_and_bad_weights():
    with pytest.raises(InputError, match="Unsupported target row family"):
        TargetRow("lithium_extraction", {"T": 298.15})

    with pytest.raises(InputError, match="binary_vle"):
        TargetRow("binary_vle", {"T": 330.0, "P": 101325.0, "x_H2O": 0.7})

    with pytest.raises(InputError, match="weight"):
        TargetRow("regularization", {"parameter": "k_ij", "target_value": 0.0}, weight=0.0)


def test_fit_problem_can_carry_target_dataset_schema_without_optimizer_internals():
    dataset = TargetDataset(
        rows=(
            TargetRow("p_rho_t", {"T": 300.0, "P": 101325.0, "rho": 55000.0}),
            TargetRow("regularization", {"parameter": "m", "target_value": 1.0}, weight=0.1),
        ),
        name="schema-only",
    )

    problem = FitProblem(mode="generic", target_dataset=dataset)

    assert problem.target_dataset is not None
    assert problem.target_dataset.families == ("p_rho_t", "regularization")
    assert problem.target_dataset.to_records()[1]["parameter"] == "m"


@pytest.mark.parametrize(
    "function, kwargs",
    [
        (
            epcsaft.fit_pure_parameters,
            {"species": "Methane", "data_rows": [{"T": 100.0}], "parameters_to_fit": ("m",)},
        ),
        (
            epcsaft.fit_binary_parameters,
            {
                "species": ("H2O", "Ethanol"),
                "data_rows": [{"T": 330.0}],
                "parameter_set": "2026_Khudaida",
                "parameters_to_fit": ("k_ij",),
            },
        ),
        (
            epcsaft.fit_liquid_electrolyte_parameters,
            {
                "species": ("H2O", "Na+", "Cl-"),
                "data_rows": [{"T": 298.15}],
                "parameter_set": "2026_Khudaida",
                "parameters_to_fit": ("d_born",),
            },
        ),
    ],
)
def test_easy_regression_apis_reject_nonexact_derivative_solver_options(function, kwargs):
    with pytest.raises(InputError, match="Unsupported derivative"):
        function(**kwargs, solver_options={"jacobian_backend": "finite" + "_difference"})
