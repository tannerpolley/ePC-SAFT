from __future__ import annotations

import json
from pathlib import Path

import pytest

import epcsaft
from tests.helpers.binary_regression_cases import ethanol_water_jced2021_vle_records


REPO_ROOT = Path(__file__).resolve().parents[2]
FIXTURE = REPO_ROOT / "tests" / "fixtures" / "literature" / "binary_kij" / "ethanol_water_jced2021_100kpa.json"
DISALLOWED_BACKENDS = {"finite_difference", "fd", "numerical_derivative", "numerical_jacobian"}


def _load_fixture() -> dict:
    with FIXTURE.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def test_literature_binary_kij_fixture_records_provenance() -> None:
    fixture = _load_fixture()

    assert fixture["family"] == "binary k_ij regression"
    for path in fixture["source"]["local_assets"]:
        assert (REPO_ROOT / path).exists()
    assert fixture["source"]["doi"] == "10.1021/acs.jced.0c00686"
    assert fixture["literature_k_ij"] == pytest.approx(-0.0269)


def test_literature_binary_kij_regression_matches_reference_band() -> None:
    fixture = _load_fixture()
    result = epcsaft.fit_binary_pair(
        ethanol_water_jced2021_vle_records(smoke_only=True),
        tuple(fixture["species"]),
        dataset=fixture["dataset"],
        initial_guess={"k_ij": fixture["initial_k_ij"]},
        bounds={"k_ij": tuple(fixture["bounds"]["k_ij"])},
        multistart=0,
    )

    assert result.success, result.message
    assert result.backend == "least_squares_native"
    assert result.problem.mode == "binary_pair"
    assert result.problem.fit_targets == ("k_ij",)
    assert result.derivative_backend.lower() not in DISALLOWED_BACKENDS
    assert result.jacobian_backend.lower() not in DISALLOWED_BACKENDS
    assert result.metrics_by_term["binary_vle_fugacity_balance"] < fixture["max_fugacity_balance"]
    assert abs(result.fitted_values["k_ij"] - fixture["literature_k_ij"]) < fixture["max_abs_kij_error"]
