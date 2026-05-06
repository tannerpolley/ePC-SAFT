import csv
import os
from pathlib import Path

import numpy as np
import pytest

from epcsaft import create_parameter_template
from epcsaft import ePCSAFTMixture
from epcsaft import fit_mea_co2_h2o_electrolyte
from epcsaft import write_fit_result

SPECIES = ["H2O", "MEA", "CO2", "MEAH+", "MEACOO-", "HCO3-"]

ADVANCED_USER_OPTIONS = {
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

REFERENCE_VALUES = {
    "MEA": {"m": 3.0353, "s": 3.0435, "e": 277.174, "e_assoc": 2586.3, "vol_a": 0.03747},
    "MEAH+": {"s": 2.78, "e": 260.0, "d_born": 3.30},
    "MEACOO-": {"s": 3.18, "e": 240.0, "d_born": 3.75},
    "HCO3-": {"s": 2.90, "e": 230.0, "d_born": 3.55},
}


def _write_pure_rows(
    dataset_root: Path, values: dict[str, dict[str, float]], *, fill_fitted_defaults: bool = True
) -> None:
    path = dataset_root / "pure" / "any_solvent.csv"
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        fieldnames = list(reader.fieldnames or [])
        rows = [dict(row) for row in reader]

    nonassoc = {"e_assoc": 0.0, "vol_a": 0.0, "assoc_scheme": "", "dielc": 8.0, "d_born": 0.0, "f_solv": 1.0}
    fixed = {
        "H2O": {
            "m": 2.1945,
            "s": 2.229,
            "e": 141.66,
            "e_assoc": 1804.17,
            "vol_a": 0.2039,
            "assoc_scheme": "4C",
            "z": 0,
            "dielc": 78.54,
            "d_born": 0.0,
            "f_solv": 1.0,
            "MW": 0.01801528,
        },
        "CO2": {**nonassoc, "m": 2.0729, "s": 2.7852, "e": 169.21, "z": 0, "MW": 0.04401},
        "MEA": {"assoc_scheme": "3B", "z": 0, "dielc": 8.0, "d_born": 0.0, "f_solv": 1.0, "MW": 0.06108},
        "MEAH+": {**nonassoc, "m": 1.0, "z": 1, "MW": 0.06209},
        "MEACOO-": {**nonassoc, "m": 1.0, "z": -1, "MW": 0.10509},
        "HCO3-": {**nonassoc, "m": 1.0, "z": -1, "MW": 0.0610168},
    }

    for row in rows:
        component = row["component"]
        base = fixed.get(component, {})
        if not fill_fitted_defaults and component in {"MEA", "MEAH+", "MEACOO-", "HCO3-"}:
            base = {
                key: value for key, value in base.items() if key not in {"m", "s", "e", "e_assoc", "vol_a", "d_born"}
            }
        merged = {**base, **values.get(component, {})}
        for key, value in merged.items():
            row[key] = value

    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def _benchmark_records(dataset_root: Path) -> list[dict[str, float]]:
    records = []
    for T, x in (
        (313.15, np.asarray([0.76, 0.18, 0.025, 0.018, 0.012, 0.005], dtype=float)),
        (333.15, np.asarray([0.72, 0.20, 0.035, 0.020, 0.016, 0.009], dtype=float)),
    ):
        mixture = ePCSAFTMixture.from_dataset(dataset_root, SPECIES, x, T, user_options=ADVANCED_USER_OPTIONS)
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


def _benchmark_dataset(tmp_path: Path, values: dict[str, dict[str, float]], name: str) -> Path:
    root = create_parameter_template(tmp_path, name, SPECIES)
    _write_pure_rows(root, values)
    return root


def test_mea_co2_h2o_benchmark_smoke_fits_only_pure_parameters(tmp_path):
    reference_root = _benchmark_dataset(tmp_path, REFERENCE_VALUES, "MEA_CO2_H2O_Reference")
    records = _benchmark_records(reference_root)
    fit_root = _benchmark_dataset(
        tmp_path,
        {
            "MEA": {"m": 2.8, "s": 2.9, "e": 250.0, "e_assoc": 2200.0, "vol_a": 0.025},
            "MEAH+": {"s": 2.5, "e": 220.0, "d_born": 2.8},
            "MEACOO-": {"s": 2.9, "e": 210.0, "d_born": 3.2},
            "HCO3-": {"s": 2.6, "e": 205.0, "d_born": 3.1},
        },
        "MEA_CO2_H2O_Fit",
    )
    before_binary = (fit_root / "mixed" / "binary_interaction" / "k_ij.csv").read_text(encoding="utf-8")

    results = fit_mea_co2_h2o_electrolyte(
        records,
        dataset=fit_root,
        species=SPECIES,
        user_options=ADVANCED_USER_OPTIONS,
        multistart=0,
    )

    assert set(results) == {"MEA", "MEAH+", "MEACOO-", "HCO3-"}
    assert results["MEA"].problem.fit_targets == ("m", "s", "e", "e_assoc", "vol_a")
    assert results["MEAH+"].problem.fit_targets == ("s", "e", "d_born")
    assert all(result.backend == "least_squares_native" for result in results.values())
    assert all(result.success for result in results.values())
    assert all("binary" not in term.term_type for result in results.values() for term in result.problem.terms)
    assert results["MEA"].metrics_by_term["initial_residual_norm"] == pytest.approx(results["MEA"].residual_norm)
    assert fit_root.joinpath("mixed", "binary_interaction", "k_ij.csv").read_text(encoding="utf-8") == before_binary


def test_mea_co2_h2o_benchmark_writes_only_pure_rows_and_protects_existing_cells(tmp_path):
    reference_root = _benchmark_dataset(tmp_path, REFERENCE_VALUES, "MEA_CO2_H2O_Reference")
    records = _benchmark_records(reference_root)
    fit_root = _benchmark_dataset(
        tmp_path,
        {
            "MEA": {"m": 2.9, "s": 2.95, "e": 260.0, "e_assoc": 2300.0, "vol_a": 0.03},
            "MEAH+": {"s": 2.6, "e": 230.0, "d_born": 3.0},
            "MEACOO-": {"s": 3.0, "e": 220.0, "d_born": 3.3},
            "HCO3-": {"s": 2.7, "e": 215.0, "d_born": 3.2},
        },
        "MEA_CO2_H2O_Fit",
    )
    results = fit_mea_co2_h2o_electrolyte(
        records,
        dataset=fit_root,
        species=SPECIES,
        user_options=ADVANCED_USER_OPTIONS,
    )

    empty_root = create_parameter_template(tmp_path, "MEA_CO2_H2O_Output", SPECIES)
    _write_pure_rows(empty_root, {"H2O": {}, "CO2": {}}, fill_fitted_defaults=False)
    written = []
    for result in results.values():
        written.extend(write_fit_result(result, empty_root, overwrite=False))

    assert written == [empty_root / "pure" / "any_solvent.csv"] * 4
    assert not any(path.name in {"k_ij.csv", "l_ij.csv", "k_hb_ij.csv"} for path in written)
    with pytest.raises(Exception, match="Refusing to overwrite"):
        write_fit_result(results["MEA"], empty_root, overwrite=False)


@pytest.mark.skipif(os.environ.get("EPCSAFT_RUN_MEA_REGRESSION") != "1", reason="opt-in MEA regression benchmark")
def test_mea_co2_h2o_benchmark_opt_in_real_multistart(tmp_path):
    reference_root = _benchmark_dataset(tmp_path, REFERENCE_VALUES, "MEA_CO2_H2O_Reference")
    records = _benchmark_records(reference_root)
    fit_root = _benchmark_dataset(
        tmp_path,
        {
            "MEA": {"m": 2.7, "s": 2.85, "e": 245.0, "e_assoc": 2100.0, "vol_a": 0.02},
            "MEAH+": {"s": 2.4, "e": 210.0, "d_born": 2.7},
            "MEACOO-": {"s": 2.8, "e": 205.0, "d_born": 3.1},
            "HCO3-": {"s": 2.5, "e": 200.0, "d_born": 3.0},
        },
        "MEA_CO2_H2O_Fit",
    )

    results = fit_mea_co2_h2o_electrolyte(
        records,
        dataset=fit_root,
        species=SPECIES,
        user_options=ADVANCED_USER_OPTIONS,
        multistart=1,
        max_nfev=2,
    )

    assert all(result.success for result in results.values())
    assert all(result.nfev > 1 for result in results.values())
    assert all(
        result.residual_norm <= result.metrics_by_term["initial_residual_norm"] + 1.0e-12 for result in results.values()
    )
