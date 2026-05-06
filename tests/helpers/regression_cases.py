from __future__ import annotations

import csv
from pathlib import Path

DATA_DIR = Path(__file__).resolve().parents[2] / "data" / "reference" / "pure_component"
REFERENCE_CSV = DATA_DIR / "hydrocarbon_basis_workbook_reference.csv"

HYDROCARBON_REFERENCE = {
    "Methane": {"m": 1.0, "s": 3.7039, "e": 150.03},
    "Ethane": {"m": 1.6069, "s": 3.5206, "e": 191.42},
    "Propane": {"m": 2.0020, "s": 3.6184, "e": 208.11},
}

HYDROCARBON_FIXED_METADATA = {
    "Methane": {"MW": 16.043e-3},
    "Ethane": {"MW": 30.070e-3},
    "Propane": {"MW": 44.097e-3},
}

REAL_DATA_TEMPERATURE_GRID = {
    "Methane": (110.0, 130.0, 150.0, 170.0),
    "Ethane": (140.0, 180.0, 220.0, 260.0),
    "Propane": (180.0, 220.0, 260.0, 300.0),
}


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


def _neutral_fixed_parameters(component: str) -> dict[str, float]:
    return {
        "MW": HYDROCARBON_FIXED_METADATA[component]["MW"],
        "e_assoc": 0.0,
        "vol_a": 0.0,
        "z": 0.0,
        "dielc": 1.0,
        "d_born": 0.0,
        "f_solv": 1.0,
    }


def _load_workbook_reference_rows() -> dict[str, dict[str, float]]:
    rows: dict[str, dict[str, float]] = {}
    with REFERENCE_CSV.open("r", encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle):
            rows[row["species"]] = {
                "m": float(row["m"]),
                "s": float(row["s_Angstrom"]),
                "e": float(row["e_over_k_K"]),
            }
    return rows


def _nist_rows_by_temperature(component: str) -> dict[float, dict[str, float]]:
    path = DATA_DIR / f"{component.lower()}_nist_saturation.csv"
    values: dict[float, dict[str, float]] = {}
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle):
            values[float(row["T_K"])] = {
                "P": float(row["p_sat_Pa"]),
                "rho_sat_liq_kg_m3": float(row["rho_sat_liq_kg_m3"]),
            }
    return values


def _real_saturation_records(component: str) -> list[dict[str, float | str]]:
    nist_rows = _nist_rows_by_temperature(component)
    records: list[dict[str, float | str]] = []
    for temperature in REAL_DATA_TEMPERATURE_GRID[component]:
        source_row = nist_rows[temperature]
        records.append(
            {
                "T": temperature,
                "P": source_row["P"],
                "rho_sat_liq_kg_m3": source_row["rho_sat_liq_kg_m3"],
                "phase": "liq",
            }
        )
    return records
