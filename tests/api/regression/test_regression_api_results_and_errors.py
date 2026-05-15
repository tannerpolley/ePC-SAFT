"""Regression API contract tests outside the hydrocarbon benchmark."""

from __future__ import annotations

import csv
from types import SimpleNamespace

import numpy as np
import pytest

import epcsaft
import epcsaft.regression as regression_module
from epcsaft import FitProblem, FitResult, create_parameter_template, fit_pure_neutral, write_fit_result
from epcsaft._types import InputError
from epcsaft.regression import (
    _debug_native_pure_neutral_objective,
    _fit_pure_neutral_least_squares_internal,
    evaluate_generic_regression_derivatives,
)
from tests.helpers.regression_cases import _methane_like_records, _minimal_neutral_metadata

def _minimal_nacl_records():
    return [
        {
            "T": 298.15,
            "P": 101325.0,
            "x_H2O": 0.996,
            "x_Na+": 0.002,
            "x_Cl-": 0.002,
            "osmotic_coefficient": 0.974,
            "mean_ionic_activity": 0.922,
        }
    ]

def _stub_native_generic_runner(monkeypatch, *, backend="least_squares_native"):
    calls = []
    jacobian_backend = "cppad_implicit" if backend == "ceres" else "stub"

    def fake_runner(
        fixed_payloads,
        native_records,
        optimization_names,
        species,
        theta0,
        lower,
        upper,
        *,
        component=None,
        pair=None,
        multistart=0,
        max_nfev=200,
    ):
        calls.append(
            {
                "fixed_payloads": fixed_payloads,
                "native_records": native_records,
                "optimization_names": tuple(optimization_names),
                "species": tuple(species),
                "theta0": np.asarray(theta0, dtype=float),
                "lower": np.asarray(lower, dtype=float),
                "upper": np.asarray(upper, dtype=float),
                "component": component,
                "pair": pair,
                "multistart": int(multistart),
                "max_nfev": int(max_nfev),
            }
        )
        metrics = {str(record["term_name"]): 0.0 for record in native_records}
        if not metrics:
            metrics = {"residual": 0.0}
        return {
            "x": np.asarray(theta0, dtype=float),
            "cost": 0.0,
            "residual_norm": 0.0,
            "initial_cost": 0.0,
            "initial_residual_norm": 0.0,
            "metrics_by_term": metrics,
            "success": True,
            "status": 1,
            "nfev": 1,
            "iterations": 0,
            "starts_tried": 1,
            "message": "stubbed native generic regression",
            "backend": backend,
            "jacobian_available": True,
            "jacobian_backend": jacobian_backend,
            "jacobian_fallback_used": False,
            "jacobian_fallback_reason": "",
            "not_available_reason": "",
            "hessian_available": False,
            "hessian_backend": "not_implemented",
            "hessian_fallback_used": False,
            "hessian_fallback_reason": "stubbed hessian skeleton",
        }

    runner_name = "_run_native_generic_ceres" if backend == "ceres" else "_run_native_generic_least_squares"
    monkeypatch.setattr(regression_module, runner_name, fake_runner)
    return calls

def test_provenance_validation_rejects_indirect_dborn_without_electrostatic_data():
    with pytest.raises(InputError, match=r"d_born.*dielectric|d_born.*ion-activity"):
        epcsaft.validate_regression_provenance(
            [
                epcsaft.FitParameter(
                    "MEAH+",
                    "d_born",
                    source="mixed_reactive_vle",
                )
            ],
            strict=True,
        )

    report = epcsaft.validate_regression_provenance(
        [
            epcsaft.FitParameter(
                "MEAH+",
                "d_born",
                source="explicit_override",
                allow_without_direct_data=True,
            )
        ],
        strict=True,
    )
    assert report["warnings"] == []
    assert report["parameter_sources"]["MEAH+.d_born"] == "explicit_override"

def test_provenance_validation_gates_ion_involving_binary_interactions():
    with pytest.raises(InputError, match="same-sign ionic pair"):
        epcsaft.validate_regression_provenance(
            [
                epcsaft.BinaryInteraction(
                    ("MEACOO-", "HCO3-"),
                    parameter="k_ij",
                    source="direct_binary_vle",
                )
            ],
            species=["MEACOO-", "HCO3-"],
            charges=[-1.0, -1.0],
            strict=True,
        )

    with pytest.raises(InputError, match=r"neutral-ion.*direct"):
        epcsaft.validate_regression_provenance(
            [
                epcsaft.BinaryInteraction(
                    ("CO2", "MEACOO-"),
                    parameter="k_ij",
                    source="mixed_reactive_vle",
                )
            ],
            species=["CO2", "MEACOO-"],
            charges=[0.0, -1.0],
            strict=True,
        )

    report = epcsaft.validate_regression_provenance(
        [
            epcsaft.BinaryInteraction(
                ("MEAH+", "MEACOO-"),
                parameter="k_ij",
                source="direct_electrolyte_activity",
            )
        ],
        species=["MEAH+", "MEACOO-"],
        charges=[1.0, -1.0],
        strict=True,
    )
    assert report["warnings"] == []
    assert report["data_sources_by_parameter"]["MEAH+:MEACOO-.k_ij"] == ["direct_electrolyte_activity"]

def test_write_fit_result_updates_only_target_pure_row(tmp_path):
    dataset_root = create_parameter_template(tmp_path, "pure_case", ["H2O", "Na+"])
    result = FitResult(
        problem=FitProblem(
            mode="pure_neutral",
            component="H2O",
            fit_targets=("m", "s", "e"),
            optimization_parameters=("m", "s", "e"),
        ),
        fitted_values={"m": 1.2047, "s": 2.7927, "e": 353.95},
        rendered_values={"m": 1.2047, "s": 2.7927, "e": 353.95},
        success=True,
    )

    written = write_fit_result(result, dataset_root, overwrite=False)
    assert written == [dataset_root / "pure" / "water.csv"]

    with (dataset_root / "pure" / "water.csv").open("r", encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle))

    target_row = next(row for row in rows if row["component"] == "H2O")
    other_row = next(row for row in rows if row["component"] == "Na+")
    assert target_row["m"] == "1.2047"
    assert target_row["s"] == "2.7927"
    assert target_row["e"] == "353.95"
    assert other_row["m"] == ""
    assert other_row["s"] == ""
    assert other_row["e"] == ""

def test_write_fit_result_updates_ion_row_and_binary_matrix_symmetrically(tmp_path):
    ion_root = create_parameter_template(tmp_path, "ion_case", ["H2O", "Na+"])
    ion_result = FitResult(
        problem=FitProblem(
            mode="pure_ion",
            component="Na+",
            fit_targets=("s", "e"),
            optimization_parameters=("s", "e"),
        ),
        fitted_values={"s": 2.84, "e": 231.2},
        rendered_values={"s": 2.84, "e": 231.2},
        success=True,
    )

    written = write_fit_result(ion_result, ion_root, overwrite=False)
    assert written == [ion_root / "pure" / "water.csv"]
    with (ion_root / "pure" / "water.csv").open("r", encoding="utf-8-sig", newline="") as handle:
        ion_rows = list(csv.DictReader(handle))
    ion_row = next(row for row in ion_rows if row["component"] == "Na+")
    assert ion_row["s"] == "2.84"
    assert ion_row["e"] == "231.2"

    binary_root = create_parameter_template(tmp_path, "binary_case", ["H2O", "Ethanol"])
    binary_result = FitResult(
        problem=FitProblem(
            mode="binary_pair",
            pair=("H2O", "Ethanol"),
            fit_targets=("k_ij", "l_ij", "k_hb_ij"),
            optimization_parameters=("k_ij", "l_ij", "k_hb_ij"),
        ),
        fitted_values={"k_ij": -0.06167, "l_ij": 0.0123, "k_hb_ij": -0.0345},
        rendered_values={"k_ij": -0.06167, "l_ij": 0.0123, "k_hb_ij": -0.0345},
        success=True,
    )

    written = write_fit_result(binary_result, binary_root, overwrite=False)
    assert written == [
        binary_root / "mixed" / "binary_interaction" / "k_ij.csv",
        binary_root / "mixed" / "binary_interaction" / "l_ij.csv",
        binary_root / "mixed" / "binary_interaction" / "k_hb_ij.csv",
    ]
    for filename, expected in (("k_ij.csv", "-0.06167"), ("l_ij.csv", "0.0123"), ("k_hb_ij.csv", "-0.0345")):
        with (binary_root / "mixed" / "binary_interaction" / filename).open(
            "r", encoding="utf-8-sig", newline=""
        ) as handle:
            rows = list(csv.reader(handle))
        header = rows[0]
        h2o_col = header.index("H2O")
        ethanol_col = header.index("Ethanol")
        h2o_row = next(row for row in rows if row[0] == "H2O")
        ethanol_row = next(row for row in rows if row[0] == "Ethanol")
        assert h2o_row[ethanol_col] == expected
        assert ethanol_row[h2o_col] == expected
    with pytest.raises(InputError, match="overwrite"):
        write_fit_result(binary_result, binary_root, overwrite=False)

def test_public_generic_derivative_evaluator_rejects_removed_backend_names():
    with pytest.raises(InputError, match="jacobian_backend"):
        evaluate_generic_regression_derivatives(
            fixed_payloads=[{"m": [1.0, 1.0], "s": [3.0, 3.0], "e": [200.0, 200.0]}],
            native_records=[],
            optimization_names=("k_ij",),
            species=("H2O", "Ethanol"),
            x=(-0.01,),
            jacobian_backend="finite" + "_difference",
        )

def test_public_generic_derivative_evaluator_rejects_auto_without_autodiff():
    with pytest.raises(InputError, match="not_available"):
        evaluate_generic_regression_derivatives(
            fixed_payloads=[{"m": [1.0, 1.0], "s": [3.0, 3.0], "e": [200.0, 200.0]}],
            native_records=[
                {
                    "term_name": "binary_vle_fugacity_balance",
                    "term": regression_module.NATIVE_TERM_KINDS["binary_vle_fugacity_balance"],
                    "T": 330.0,
                    "P": 101325.0,
                    "x": [0.7, 0.3],
                    "y": [0.5, 0.5],
                }
            ],
            optimization_names=("k_ij", "l_ij"),
            species=("H2O", "Ethanol"),
            pair=("H2O", "Ethanol"),
            x=(-0.01, 0.02),
        )
