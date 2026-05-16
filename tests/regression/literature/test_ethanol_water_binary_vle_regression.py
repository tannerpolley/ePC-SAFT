from __future__ import annotations

import csv

import pytest

import epcsaft
from tests.helpers.binary_regression_cases import (
    ETHANOL_WATER_HELD_2012_KIJ,
    ETHANOL_WATER_JCED2021_100KPA,
    ETHANOL_WATER_PAPER_PCSAFT_KIJ_100KPA,
    ethanol_water_jced2021_vle_records,
)


def test_ethanol_water_jced2021_vle_dataset_loads_as_binary_regression_records() -> None:
    records = ethanol_water_jced2021_vle_records()
    smoke_records = ethanol_water_jced2021_vle_records(smoke_only=True)

    assert ETHANOL_WATER_JCED2021_100KPA.exists()
    assert len(records) == 26
    assert len(smoke_records) == 5
    assert {record["source_doi"] for record in records} == {"10.1021/acs.jced.0c00686"}
    assert all(record["P"] == pytest.approx(100000.0) for record in records)

    with ETHANOL_WATER_JCED2021_100KPA.open("r", encoding="utf-8-sig", newline="") as handle:
        header = next(csv.reader(handle))
    assert {"T", "P", "x_Ethanol", "x_H2O", "y_Ethanol", "y_H2O"}.issubset(header)


def test_native_binary_kij_regression_matches_real_ethanol_water_vle_reference_band() -> None:
    result = epcsaft.fit_binary_pair(
        ethanol_water_jced2021_vle_records(smoke_only=True),
        ("Ethanol", "H2O"),
        dataset="2012_Held",
        initial_guess={"k_ij": ETHANOL_WATER_HELD_2012_KIJ},
        bounds={"k_ij": (-0.15, 0.10)},
    )

    assert result.success, result.message
    assert result.backend == "ceres"
    assert result.optimizer_backend == "ceres"
    assert result.derivative_backend == "cppad_implicit"
    assert result.problem.mode == "binary_pair"
    assert result.problem.fit_targets == ("k_ij",)
    assert result.objective_final < result.objective_initial
    assert result.parameter_movement["k_ij"] == pytest.approx(
        result.fitted_values["k_ij"] - ETHANOL_WATER_HELD_2012_KIJ
    )
    assert abs(result.parameter_movement["k_ij"]) > 1.0e-5
    assert result.initial_parameters == {"k_ij": pytest.approx(ETHANOL_WATER_HELD_2012_KIJ)}
    assert result.final_parameters == pytest.approx(result.fitted_values)
    assert result.source_summaries["records"]["record_count"] == 5
    assert result.residual_block_norms["binary_vle_fugacity_balance"] == pytest.approx(
        result.metrics_by_term["binary_vle_fugacity_balance"]
    )
    assert result.metrics_by_term["binary_vle_fugacity_balance"] < 0.04
    assert abs(result.fitted_values["k_ij"] - ETHANOL_WATER_PAPER_PCSAFT_KIJ_100KPA) < 0.01
