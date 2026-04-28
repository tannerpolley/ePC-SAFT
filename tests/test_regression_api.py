# -*- coding: utf-8 -*-
"""Regression API contract tests outside the hydrocarbon benchmark."""

from __future__ import annotations

import csv

import epcsaft
import numpy as np
import pytest

from epcsaft import FitProblem
from epcsaft import FitResult
from epcsaft import create_parameter_template
from epcsaft import fit_pure_neutral
from epcsaft import write_fit_result
from epcsaft._types import InputError
from epcsaft.regression import _debug_native_pure_neutral_objective
from epcsaft.regression import _fit_pure_neutral_least_squares_internal
from epcsaft.regression import fit_binary_pair
from epcsaft.regression import fit_pure_ion


def _minimal_neutral_metadata(mw: float) -> dict[str, float]:
    return {
        "MW": mw,
        "e_assoc": 0.0,
        "vol_a": 0.0,
        "z": 0.0,
        "dielc": 1.0,
        "d_born": 0.0,
        "f_solv": 1.0,
    }


def _methane_like_records() -> list[dict[str, float | str]]:
    return [
        {"T": 110.0, "P": 88130.038, "rho_sat_liq_kg_m3": 424.77725, "phase": "liq"},
        {"T": 130.0, "P": 367319.94, "rho_sat_liq_kg_m3": 394.35230, "phase": "liq"},
    ]


def test_public_regression_surface_is_neutral_only():
    assert hasattr(epcsaft, "fit_pure_neutral")
    assert not hasattr(epcsaft, "fit_pure_ion")
    assert not hasattr(epcsaft, "fit_binary_pair")


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



@pytest.mark.parametrize(
    "fn, kwargs",
    [
        (
            fit_pure_ion,
            {
                "records": [{"T": 298.15, "P": 1.0e5}],
                "component": "Na+",
                "dataset": "2012_Held",
            },
        ),
        (
            fit_binary_pair,
            {
                "records": [{"T": 360.0, "P": 1.0e5}],
                "pair": ["Benzene", "Toluene"],
                "dataset": "2012_Held",
            },
        ),
    ],
)
def test_deferred_regression_workflows_raise(fn, kwargs):
    with pytest.raises(NotImplementedError, match="Phase 1 supports only neutral-component fitting of m, s, and e"):
        fn(**kwargs)


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
