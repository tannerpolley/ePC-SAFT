"""API contract tests for user-facing regression problem schemas."""

from __future__ import annotations

import pytest

import epcsaft
import epcsaft.regression as regression_module
from epcsaft import FitProblem, FitResult
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


def test_fit_liquid_electrolyte_parameters_returns_stable_backend_unavailable_contract():
    result = epcsaft.fit_liquid_electrolyte_parameters(
        species=("H2O", "Na+", "Cl-"),
        data_rows=[{"T": 298.15, "P": 101325.0, "molality": 0.1, "osmotic_coefficient": 0.933}],
        parameter_set="2026_Khudaida",
        parameters_to_fit=("d_born", "f_solv"),
        fixed_parameters={"H2O": {"m": 1.0}},
        weights={"osmotic_coefficient": 1.0},
        loss="linear",
        solver_options={"max_nfev": 1},
    )

    assert result.success is False
    assert result.status == -1
    assert result.backend == "backend_unavailable"
    assert result.backend_unavailable_reason
    assert result.problem.mode == "liquid_electrolyte"
    assert result.problem.fit_targets == ("d_born", "f_solv")
    assert result.problem.dataset == "2026_Khudaida"
    assert result.problem.weights == {"osmotic_coefficient": 1.0}


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
def test_easy_regression_apis_reject_finite_difference_solver_options(function, kwargs):
    with pytest.raises(InputError, match="finite-difference"):
        function(**kwargs, solver_options={"jacobian_backend": "finite" + "_difference"})
