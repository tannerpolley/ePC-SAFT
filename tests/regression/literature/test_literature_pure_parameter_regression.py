from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from epcsaft import ePCSAFTMixture
from epcsaft.regression import _fit_mea_co2_h2o_pure_parameter_benchmark
from tests.regression.literature.test_mea_co2_h2o_pure_parameter_benchmark import (
    SPECIES,
    _benchmark_dataset,
)

REPO_ROOT = Path(__file__).resolve().parents[3]
FIXTURE = REPO_ROOT / "tests" / "fixtures" / "literature" / "pure_neutral" / "mea_co2_h2o_benchmark.json"
DISALLOWED_BACKENDS = {"numerical_derivative", "fd", "numerical_jacobian"}
LITERATURE_USER_OPTIONS = {
    "elec_model": {
        "rel_perm": {"rule": "empirical", "differential_mode": "auto"},
        "born_model": {
            "d_Born_mode": 3,
            "solvation_shell_model": True,
            "dielectric_saturation": True,
            "mu_born_model": {"differential_mode": "auto", "comp_dep_delta_d": True},
        },
    }
}


def _load_fixture() -> dict:
    with FIXTURE.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def _benchmark_records(dataset_root: Path) -> list[dict[str, float]]:
    records = []
    for T, x in (
        (313.15, np.asarray([0.76, 0.18, 0.025, 0.018, 0.012, 0.005], dtype=float)),
        (333.15, np.asarray([0.72, 0.20, 0.035, 0.020, 0.016, 0.009], dtype=float)),
    ):
        mixture = ePCSAFTMixture.from_dataset(dataset_root, SPECIES, x, T, user_options=LITERATURE_USER_OPTIONS)
        state = mixture.state(T, x, P=101325.0, phase="liq")
        lnphi = state.fugacity_coefficient(natural_log=True)
        row = {"T": T, "P": 101325.0}
        for species, value in zip(SPECIES, x):
            row[f"x_{species}"] = float(value)
        row["rho"] = float(state.molar_density())
        row["lnphi_CO2"] = float(lnphi[SPECIES.index("CO2")])
        row["osmotic_coefficient"] = float(state.osmotic_coefficient()[0])
        records.append(row)
    return records


def test_literature_pure_parameter_fixture_records_provenance() -> None:
    fixture = _load_fixture()

    assert fixture["family"] == "pure neutral parameter regression"
    assert fixture["source"]["local_asset"].startswith("docs/papers/")
    assert (REPO_ROOT / fixture["source"]["local_asset"]).exists()
    assert fixture["validation_role"] == "fast smoke benchmark"
    assert set(fixture["components"]) == {"MEA", "MEAH+", "MEACOO-", "HCO3-"}


def test_literature_pure_parameter_regression_uses_local_benchmark_fixture(tmp_path) -> None:
    fixture = _load_fixture()
    reference_root = _benchmark_dataset(tmp_path, fixture["reference_values"], "MEA_CO2_H2O_Reference")
    records = _benchmark_records(reference_root)
    fit_root = _benchmark_dataset(tmp_path, fixture["initial_guess"], "MEA_CO2_H2O_Fit")

    results = _fit_mea_co2_h2o_pure_parameter_benchmark(
        records,
        dataset=fit_root,
        species=SPECIES,
        user_options=LITERATURE_USER_OPTIONS,
        multistart=0,
    )

    assert set(results) == set(fixture["components"])
    for component, result in results.items():
        assert result.success, result.message
        assert result.backend == "residual_score_native"
        assert result.optimizer_backend == "residual_score_native"
        assert result.derivative_backend.lower() not in DISALLOWED_BACKENDS
        assert result.jacobian_backend.lower() not in DISALLOWED_BACKENDS
        assert result.problem.fit_targets == tuple(fixture["fit_targets"][component])
        expected_mode = "pure_neutral" if component == "MEA" else "pure_ion"
        assert result.problem.mode == expected_mode
