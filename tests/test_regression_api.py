# -*- coding: utf-8 -*-
"""Regression API contract tests outside the hydrocarbon benchmark."""

from __future__ import annotations

import csv

import epcsaft
import pytest

from epcsaft import FitProblem
from epcsaft import FitResult
from epcsaft import create_parameter_template
from epcsaft import fit_pure_neutral
from epcsaft import write_fit_result
from epcsaft._types import InputError
from epcsaft.regression import fit_binary_pair
from epcsaft.regression import fit_pure_ion


def _minimal_neutral_metadata(mw: float) -> dict[str, float]:
    return {
        "MW": mw,
        "e_assoc": 0.0,
        "vol_a": 0.0,
        "dipm": 0.0,
        "dip_num": 1.0,
        "z": 0.0,
        "dielc": 1.0,
        "d_born": 0.0,
        "f_solv": 1.0,
    }


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
