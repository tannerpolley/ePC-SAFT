# -*- coding: utf-8 -*-
"""Regression API contract tests outside the hydrocarbon benchmark."""

from __future__ import annotations

import csv

import epcsaft
import numpy as np
import pytest

from epcsaft import ePCSAFTMixture
from epcsaft import FitProblem
from epcsaft import FitResult
from epcsaft import create_parameter_template
from epcsaft import fit_pure_neutral
from epcsaft import molality_to_molefraction
from epcsaft import write_fit_result
from epcsaft._types import InputError
from epcsaft.regression import _debug_native_pure_neutral_objective
from epcsaft.regression import _fit_pure_neutral_least_squares_internal
from tests.helpers.regression_cases import _methane_like_records
from tests.helpers.regression_cases import _minimal_neutral_metadata


def test_public_regression_surface_includes_ion_and_binary_v1():
    assert hasattr(epcsaft, "fit_pure_neutral")
    assert hasattr(epcsaft, "fit_pure_ion")
    assert hasattr(epcsaft, "fit_binary_pair")


def test_fit_pure_neutral_requires_pressure_for_density_records():
    with pytest.raises(InputError, match="experimental 'P'"):
        fit_pure_neutral(
            [{"T": 320.0, "rho": 9000.0}],
            "Toluene",
            assoc_scheme="",
            fixed_parameters=_minimal_neutral_metadata(92.141e-3),
            initial_guess={"m": 2.8, "s": 3.7, "e": 285.0},
        )


def test_fit_pure_neutral_rejects_non_phase1_targets():
    with pytest.raises(InputError, match="supports only the targets 'm', 's', and 'e'"):
        fit_pure_neutral(
            [{"T": 320.0, "P": 101325.0, "rho": 9000.0}],
            "Toluene",
            assoc_scheme="",
            fit_targets=("m", "s", "e_assoc"),
            fixed_parameters=_minimal_neutral_metadata(92.141e-3),
            initial_guess={"m": 2.8, "s": 3.7, "e_assoc": 1000.0},
        )


def test_native_pure_neutral_debug_gradient_matches_finite_difference():
    theta = {"m": 1.05, "s": 3.68, "e": 151.0}

    def objective_at(m: float, s: float, e: float) -> float:
        debug = _debug_native_pure_neutral_objective(
            _methane_like_records(),
            "Methane",
            assoc_scheme="",
            fixed_parameters=_minimal_neutral_metadata(16.043e-3),
            initial_guess=theta,
            x={"m": m, "s": s, "e": e},
        )
        return float(debug["objective"])

    debug = _debug_native_pure_neutral_objective(
        _methane_like_records(),
        "Methane",
        assoc_scheme="",
        fixed_parameters=_minimal_neutral_metadata(16.043e-3),
        initial_guess=theta,
        x=theta,
    )
    exact = np.asarray(debug["gradient"], dtype=float)

    eps = np.asarray([1.0e-6, 1.0e-6, 1.0e-5], dtype=float)
    fd = np.empty(3, dtype=float)
    base = np.asarray([theta["m"], theta["s"], theta["e"]], dtype=float)
    for i in range(3):
        forward = base.copy()
        backward = base.copy()
        forward[i] += eps[i]
        backward[i] -= eps[i]
        fd[i] = (
            objective_at(*forward) - objective_at(*backward)
        ) / (2.0 * eps[i])

    assert exact == pytest.approx(fd, rel=5.0e-4, abs=5.0e-6)
    assert debug["residual_evaluations"] >= 1
    assert debug["density_solves"] >= 2
    assert debug["fused_state_evaluations"] >= 2
    assert debug["callback_wall_time_s"] >= 0.0


def test_internal_native_least_squares_backend_matches_methane_reference_band():
    result = _fit_pure_neutral_least_squares_internal(
        _methane_like_records(),
        "Methane",
        assoc_scheme="",
        fixed_parameters=_minimal_neutral_metadata(16.043e-3),
        initial_guess={"m": 1.08, "s": 3.55, "e": 155.0},
        bounds={
            "m": (0.5, 3.5),
            "s": (2.0, 5.0),
            "e": (50.0, 400.0),
        },
    )

    assert result.success, result.message
    assert result.backend == "least_squares_native"
    assert result.metrics_by_term["density"] < 0.02
    assert result.metrics_by_term["pure_vle_fugacity_balance"] < 0.02
    assert result.fitted_values["m"] == pytest.approx(1.0, rel=0.0, abs=0.05)
    assert result.fitted_values["s"] == pytest.approx(3.7039, rel=0.0, abs=0.08)
    assert result.fitted_values["e"] == pytest.approx(150.03, rel=0.0, abs=3.0)



@pytest.mark.parametrize(
    "initial_guess",
    [
        {"m": 1.08, "s": 3.55, "e": 155.0},
        {"m": 0.92, "s": 3.90, "e": 143.0},
        {"m": 1.20, "s": 3.30, "e": 170.0},
    ],
)
def test_public_pure_neutral_regression_is_robust_to_distinct_initial_guesses(initial_guess):
    result = fit_pure_neutral(
        _methane_like_records(),
        "Methane",
        assoc_scheme="",
        fixed_parameters=_minimal_neutral_metadata(16.043e-3),
        initial_guess=initial_guess,
        bounds={
            "m": (0.5, 3.5),
            "s": (2.0, 5.0),
            "e": (50.0, 400.0),
        },
    )
    assert result.success, result.message
    assert result.metrics_by_term["density"] < 0.02
    assert result.metrics_by_term["pure_vle_fugacity_balance"] < 0.02
    assert result.fitted_values["m"] == pytest.approx(1.0, rel=0.0, abs=0.06)
    assert result.fitted_values["s"] == pytest.approx(3.7039, rel=0.0, abs=0.10)
    assert result.fitted_values["e"] == pytest.approx(150.03, rel=0.0, abs=4.0)



def _nacl_records(*, user_options=None):
    species = ["H2O", "Na+", "Cl-"]
    records = []
    for molality in (0.1, 0.2):
        x = molality_to_molefraction(molality, species=species, solvent="H2O")
        mixture = ePCSAFTMixture.from_dataset("2026_Khudaida", species, x, 298.15, user_options=user_options)
        state = mixture.state(298.15, x, P=101325.0, phase="liq")
        miac = state.activity_coefficient(species=species, mean_ionic_form=True, basis="molality")
        records.append(
            {
                "T": 298.15,
                "P": 101325.0,
                "molality": molality,
                "osmotic_coefficient": float(state.osmotic_coefficient()[0]),
                "mean_ionic_activity": float(miac["Na+Cl-"]),
            }
        )
    return records


def test_fit_pure_ion_requires_composition_and_activity_or_osmotic_records():
    with pytest.raises(InputError, match="composition"):
        epcsaft.fit_pure_ion(
            [{"T": 298.15, "P": 101325.0, "osmotic_coefficient": 0.93}],
            "Na+",
            dataset="2026_Khudaida",
        )

    with pytest.raises(InputError, match="osmotic|mean-ionic|mean ionic"):
        epcsaft.fit_pure_ion(
            [{"T": 298.15, "P": 101325.0, "molality": 0.1}],
            "Na+",
            dataset="2026_Khudaida",
            species=["H2O", "Na+", "Cl-"],
            solvent="H2O",
        )


def test_fit_pure_ion_default_s_e_bounds_and_deterministic_multistart():
    result = epcsaft.fit_pure_ion(
        _nacl_records(),
        "Na+",
        dataset="2026_Khudaida",
        species=["H2O", "Na+", "Cl-"],
        solvent="H2O",
        initial_guess={"s": 2.6, "e": 210.0},
        bounds={"s": (2.4, 3.2), "e": (150.0, 300.0)},
        multistart=3,
    )
    repeat = epcsaft.fit_pure_ion(
        _nacl_records(),
        "Na+",
        dataset="2026_Khudaida",
        species=["H2O", "Na+", "Cl-"],
        solvent="H2O",
        initial_guess={"s": 2.6, "e": 210.0},
        bounds={"s": (2.4, 3.2), "e": (150.0, 300.0)},
        multistart=3,
    )

    assert result.success, result.message
    assert result.backend == "least_squares_native"
    assert result.problem.mode == "pure_ion"
    assert result.problem.fit_targets == ("s", "e")
    assert result.metrics_by_term["osmotic_coefficient"] < 1.0e-3
    assert result.metrics_by_term["mean_ionic_activity"] < 2.0e-3
    assert result.fitted_values == pytest.approx(repeat.fitted_values, rel=0.0, abs=1.0e-12)


def test_fit_pure_ion_accepts_d_born_and_born_user_options():
    user_options = {
        "elec_model": {
            "rel_perm": {"rule": "empirical", "differential_mode": "numerical"},
            "born_model": {
                "d_Born_mode": 3,
                "solvation_shell_model": True,
                "dielectric_saturation": True,
                "mu_born_model": {"differential_mode": "numerical", "comp_dep_delta_d": True},
            },
        }
    }
    result = epcsaft.fit_pure_ion(
        _nacl_records(user_options=user_options),
        "Na+",
        dataset="2026_Khudaida",
        species=["H2O", "Na+", "Cl-"],
        solvent="H2O",
        fit_targets=("d_born",),
        initial_guess={"d_born": 3.2},
        bounds={"d_born": (2.0, 5.0)},
        user_options=user_options,
    )

    assert result.success, result.message
    assert result.backend == "least_squares_native"
    assert result.problem.fit_targets == ("d_born",)
    assert result.metrics_by_term["osmotic_coefficient"] < 1.0e-3


def test_fit_pure_ion_passes_explicit_mean_ionic_pair_label_to_native_backend():
    records = [dict(record, pair_label="Na+Cl-") for record in _nacl_records()]
    result = epcsaft.fit_pure_ion(
        records,
        "Na+",
        dataset="2026_Khudaida",
        species=["H2O", "Na+", "Cl-"],
        solvent="H2O",
        initial_guess={"s": 2.6, "e": 210.0},
        bounds={"s": (2.4, 3.2), "e": (150.0, 300.0)},
    )

    assert result.success, result.message
    assert result.backend == "least_squares_native"
    assert result.metrics_by_term["mean_ionic_activity"] < 2.0e-3


def test_fit_binary_pair_vle_kij_default_and_rejects_temperature_models():
    records = [
        {"T": 330.0, "P": 101325.0, "x_H2O": 0.7, "x_Ethanol": 0.3, "y_H2O": 0.5, "y_Ethanol": 0.5},
        {"T": 340.0, "P": 101325.0, "x_H2O": 0.6, "x_Ethanol": 0.4, "y_H2O": 0.4, "y_Ethanol": 0.6},
    ]
    result = epcsaft.fit_binary_pair(
        records,
        ("H2O", "Ethanol"),
        dataset="2026_Khudaida",
        initial_guess={"k_ij": -0.02},
        bounds={"k_ij": (-0.2, 0.2)},
        multistart=2,
    )

    assert result.success, result.message
    assert result.backend == "least_squares_native"
    assert result.problem.mode == "binary_pair"
    assert result.problem.fit_targets == ("k_ij",)
    assert set(result.fitted_values) == {"k_ij"}
    assert "binary_vle_fugacity_balance" in result.metrics_by_term

    with pytest.raises(InputError, match="temperature_model"):
        epcsaft.fit_binary_pair(
            records,
            ("H2O", "Ethanol"),
            dataset="2026_Khudaida",
            temperature_model="linear",
        )


def test_fit_binary_pair_can_fit_all_constant_binary_interaction_targets():
    records = [
        {"T": 330.0, "P": 101325.0, "x_H2O": 0.7, "x_Ethanol": 0.3, "y_H2O": 0.5, "y_Ethanol": 0.5},
        {"T": 340.0, "P": 101325.0, "x_H2O": 0.6, "x_Ethanol": 0.4, "y_H2O": 0.4, "y_Ethanol": 0.6},
    ]

    result = epcsaft.fit_binary_pair(
        records,
        ("H2O", "Ethanol"),
        dataset="2026_Khudaida",
        fit_targets=("k_ij", "l_ij", "k_hb_ij"),
        initial_guess={"k_ij": -0.02, "l_ij": 0.01, "k_hb_ij": 0.02},
        bounds={"k_ij": (-0.2, 0.2), "l_ij": (-0.2, 0.2), "k_hb_ij": (-0.2, 0.2)},
        multistart=1,
    )

    assert result.success, result.message
    assert result.backend == "least_squares_native"
    assert result.problem.fit_targets == ("k_ij", "l_ij", "k_hb_ij")
    assert set(result.fitted_values) == {"k_ij", "l_ij", "k_hb_ij"}
    assert all(np.isfinite(value) for value in result.fitted_values.values())
    assert "binary_vle_fugacity_balance" in result.metrics_by_term


def test_fit_binary_pair_rejects_nonpositive_vle_fractions_before_native_solve():
    records = [
        {"T": 330.0, "P": 101325.0, "x_H2O": 0.7, "x_Ethanol": 0.3, "y_H2O": 0.0, "y_Ethanol": 1.0},
    ]
    with pytest.raises(InputError, match="strictly positive"):
        epcsaft.fit_binary_pair(
            records,
            ("H2O", "Ethanol"),
            dataset="2026_Khudaida",
            initial_guess={"k_ij": -0.02},
        )


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
        with (binary_root / "mixed" / "binary_interaction" / filename).open("r", encoding="utf-8-sig", newline="") as handle:
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
