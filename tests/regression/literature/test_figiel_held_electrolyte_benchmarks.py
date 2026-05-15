from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
import pytest

import epcsaft
from epcsaft.parameters import molality_to_molefraction

REPO_ROOT = Path(__file__).resolve().parents[3]
FIGIEL_FIXTURE = REPO_ROOT / "tests" / "fixtures" / "literature" / "figiel_2025" / "miac_liquid_electrolyte.json"
HELD_OSMOTIC_NACL = REPO_ROOT / "data" / "reference" / "osmotic" / "water" / "NaCl.csv"
HELD_MIAC_WATER_NACL = REPO_ROOT / "data" / "reference" / "MIAC" / "water" / "water-NaCl.csv"
HELD_MIAC_WATER_METHANOL_NACL = REPO_ROOT / "data" / "reference" / "MIAC" / "water-methanol" / "water-methanol-NaCl.csv"
P_REF = 101325.0
T_REF = 298.15
SOLVENT_MW = {"water": 0.01801528, "methanol": 0.03204}


def _csv_row_by_molality(path: Path, target: float) -> dict[str, str]:
    rows = list(csv.DictReader(path.open("r", encoding="utf-8-sig", newline="")))
    key = "molality" if "molality" in rows[0] else "m (mol/kg)"
    return min(rows, key=lambda row: abs(float(row[key]) - target))


def _water_methanol_nacl_x(row: dict[str, str]) -> np.ndarray:
    x_h2o = float(row["x_H2O"])
    x_meoh = float(row["x_Methanol"])
    molality = float(row["molality"])
    solvent_mw = x_h2o * SOLVENT_MW["water"] + x_meoh * SOLVENT_MW["methanol"]
    solvent_moles = 1.0 / solvent_mw
    moles = {
        "Na+": molality,
        "Cl-": molality,
        "H2O": x_h2o * solvent_moles,
        "Methanol": x_meoh * solvent_moles,
    }
    total = sum(moles.values())
    return np.asarray([moles[name] / total for name in ("Na+", "Cl-", "H2O", "Methanol")], dtype=float)


def test_figiel_liquid_electrolyte_probe_reports_native_ceres_benchmark_evidence() -> None:
    fixture = json.loads(FIGIEL_FIXTURE.read_text(encoding="utf-8"))
    probe = fixture["regression_probe"]

    result = epcsaft.fit_liquid_electrolyte_parameters(
        species=tuple(probe["species"]),
        data_rows=probe["data_rows"],
        parameter_set=probe["parameter_set"],
        parameters_to_fit=tuple(probe["parameters_to_fit"]),
        initial_guess=probe["initial_guess"],
        bounds={key: tuple(value) for key, value in probe["bounds"].items()},
        solver_options={"optimizer_backend": "ceres", "max_nfev": 50},
    )

    assert result.success, result.message
    assert result.backend == "ceres"
    assert result.jacobian_backend == "cppad_implicit"
    assert result.python_objective_used is False
    assert result.objective_final <= result.objective_initial
    assert set(result.parameter_movement) == {"d_born", "f_solv"}
    assert any(abs(value) > 1.0e-8 for value in result.parameter_movement.values())
    assert result.problem.solver_options["target_components"] == {"d_born": "Na+", "f_solv": "H2O"}


def test_held_cameretti_aqueous_nacl_density_osmotic_miac_benchmark() -> None:
    osmotic_row = _csv_row_by_molality(HELD_OSMOTIC_NACL, 0.5)
    miac_row = _csv_row_by_molality(HELD_MIAC_WATER_NACL, 0.5)
    molality = float(osmotic_row["m (mol/kg)"])
    species = ("H2O", "Na+", "Cl-")
    x = molality_to_molefraction(molality, species=species, solvent="H2O")

    mixture = epcsaft.ePCSAFTMixture.from_dataset("2014_Held", species, x, T_REF)
    state = mixture.state(T=T_REF, x=x, P=P_REF, phase="liq")
    osmotic_calc = float(state.osmotic_coefficient()[0])
    miac_calc = state.activity_coefficient(mean_ionic_form=True, basis="molality")["Na+Cl-"]

    assert state.molar_density() > 0.0
    assert osmotic_calc == pytest.approx(float(osmotic_row["osmotic"]), abs=0.03)
    assert miac_calc == pytest.approx(float(miac_row["miac_m"]), abs=0.08)


def test_held_2012_water_methanol_nacl_density_osmotic_miac_benchmark() -> None:
    row = _csv_row_by_molality(HELD_MIAC_WATER_METHANOL_NACL, 0.0125)
    species = ("Na+", "Cl-", "H2O", "Methanol")
    x = _water_methanol_nacl_x(row)

    mixture = epcsaft.ePCSAFTMixture.from_dataset("2012_Held", species, x, T_REF)
    state = mixture.state(T=T_REF, x=x, P=P_REF, phase="liq")
    osmotic_calc = float(state.osmotic_coefficient()[0])
    miac_calc = state.activity_coefficient(mean_ionic_form=True, basis="molality")["Na+Cl-"]

    assert state.molar_density() > 0.0
    assert np.isfinite(osmotic_calc)
    assert 0.0 < osmotic_calc < 2.0
    assert miac_calc == pytest.approx(float(row["miac_m"]), abs=0.08)
