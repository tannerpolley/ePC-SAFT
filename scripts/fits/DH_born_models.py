# -*- coding: utf-8 -*-
"""Born diagnostics for MIAC with dielc rules 0-6 across born models."""

import csv
import os
from pathlib import Path

import numpy as np
import pytest

from pcsaft import pcsaft_den, pcsaft_miac_m
from data.epcsaft_properties import get_prop_dict, molality_to_molefraction


WATER_SPECIES_MAP = {
    "LiCl": ["Li+", "Cl-", "H2O-2B-Li"],
    "NaCl": ["Na+", "Cl-", "H2O-2B-Li"],
    "KCl": ["K+", "Cl-", "H2O-2B-Li"],
    "LiBr": ["Li+", "Br-", "H2O-2B-Li"],
    "NaBr": ["Na+", "Br-", "H2O-2B-Li"],
    "KBr": ["K+", "Br-", "H2O-2B-Li"],
}

WATER_SALTS = ["LiCl", "NaCl", "KCl", "LiBr", "NaBr", "KBr"]
BORN_MODELS = [0, 1, 2, 3, 4, 5]
DIELC_RULES = [0, 1, 2, 3, 4, 5, 6]

def _selected_solvent():
    solvent = os.getenv("MIAC_SOLVENT", "methanol").strip().lower()
    if solvent not in {"water", "methanol"}:
        raise ValueError(f"Unsupported MIAC_SOLVENT='{solvent}'. Supported values: water, methanol.")
    return solvent


def _selected_salt():
    salt = os.getenv("MIAC_SALT", "NaBr").strip()
    if not salt:
        raise ValueError("MIAC_SALT is empty.")
    return salt


def _axis_limits(solvent, molal_exp, gamma_exp):
    if solvent == "water":
        xmax = max(1, int(np.ceil(float(np.max(molal_exp)))))
        ymax = max(1, int(np.ceil(float(np.max(gamma_exp)))))
        return (0.0, float(xmax)), (0.0, float(ymax))
    return (0.0, 2.0), (0.0, 1.0)


def _parse_methanol_salt(salt):
    cations = ("Li", "Na", "K")
    anions = {"Cl", "Br", "I"}
    for cat in cations:
        if salt.startswith(cat):
            an = salt[len(cat) :]
            if an in anions:
                return f"{cat}+", f"{an}-"
    raise ValueError(f"Unsupported methanol salt '{salt}'. Expected forms like LiCl, NaBr, KI.")


def _species_for_combo(solvent, salt):
    if solvent == "water":
        if salt not in WATER_SPECIES_MAP:
            raise ValueError(f"Unsupported water salt '{salt}'. Supported salts: {sorted(WATER_SPECIES_MAP)}.")
        return WATER_SPECIES_MAP[salt]
    cation, anion = _parse_methanol_salt(salt)
    return [cation, anion, "Methanol"]


def _load_miac_data(solvent, salt):
    root = Path(__file__).resolve().parents[2]
    if solvent == "water":
        m_min, m_max = 0.0, 6.0
    else:
        m_min, m_max = 0.0, 2.0
    molal = []
    gamma_exp = []
    if solvent == "water":
        data_candidates = [
            root / "data" / "MIAC_m" / "water" / f"{salt}.csv",
            root / "data" / "MIAC_m" / "csv" / f"{salt}.csv",
        ]
        data_path = next((p for p in data_candidates if p.exists()), None)
        if data_path is None:
            raise FileNotFoundError(f"Water MIAC data file not found in candidates: {data_candidates}")
        with data_path.open("r", newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle)
            if "m" not in reader.fieldnames or "gamma" not in reader.fieldnames:
                raise ValueError(f"Missing required columns m/gamma in {data_path}.")
            for row in reader:
                m = float(row["m"])
                g = float(row["gamma"])
                if m_min <= m <= m_max:
                    molal.append(m)
                    gamma_exp.append(g)
    else:
        data_candidates = [
            root / "data" / "MIAC_m" / "methanol" / f"methanol-{salt}.csv",
            root / "data" / "MIAC_m" / "water_methanol" / f"water-methanol-{salt}.csv",
            root / "data" / "MIAC_m" / "water_methanol" / f"{salt}.csv",
        ]
        data_path = next((p for p in data_candidates if p.exists()), None)
        if data_path is None:
            raise FileNotFoundError(f"Methanol MIAC data file not found in candidates: {data_candidates}")
        with data_path.open("r", newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle)
            fields = set(reader.fieldnames or [])
            y_key = "gamma" if "gamma" in fields else ("miac_m" if "miac_m" in fields else None)
            if {"molality (kg/mol)"}.issubset(fields) and y_key is not None:
                for row in reader:
                    m = float(row["molality (kg/mol)"])
                    g = float(row[y_key])
                    if m_min <= m <= m_max:
                        molal.append(m)
                        gamma_exp.append(g)
            elif {"x_methanol", "molality (kg/mol)"}.issubset(fields) and y_key is not None:
                for row in reader:
                    x_meoh = float(row["x_methanol"])
                    m = float(row["molality (kg/mol)"])
                    g = float(row[y_key])
                    if abs(x_meoh - 1.0) <= 1e-12 and m_min <= m <= m_max:
                        molal.append(m)
                        gamma_exp.append(g)
            elif {"methanol_pct", "molal", "gamma_exp"}.issubset(fields):
                for row in reader:
                    pct = float(row["methanol_pct"])
                    m = float(row["molal"])
                    g = float(row["gamma_exp"])
                    if pct == 100.0 and m_min <= m <= m_max:
                        molal.append(m)
                        gamma_exp.append(g)
            else:
                raise ValueError(f"Unsupported methanol schema in {data_path}: {sorted(fields)}")

    if not molal:
        raise ValueError(f"No data in range {m_min}<=m<={m_max} for solvent={solvent}, salt={salt}.")

    molal = np.asarray(molal, dtype=float)
    gamma_exp = np.asarray(gamma_exp, dtype=float)
    if not np.all(np.isfinite(molal)) or not np.all(np.isfinite(gamma_exp)):
        raise ValueError(f"Non-finite experimental data found for solvent={solvent}, salt={salt}.")

    order = np.argsort(molal)
    return molal[order], gamma_exp[order]


def _calc_rule_curve(molal, species, dielc_rule, born_model):
    t = 298.15
    p = 101325.0
    gamma_calc = np.empty_like(molal, dtype=float)

    for i, m_salt in enumerate(molal):
        x = molality_to_molefraction(float(m_salt), species=species)
        born_model_int = int(born_model)
        born_radius_model = 5 if born_model_int >= 2 else 1
        params = get_prop_dict(
            species,
            x,
            t,
            user_options={
                "dielc_rule": int(dielc_rule),
                "born_model": born_model_int,
                "born_radius_model": born_radius_model,
                "DH_model": 1,
                "debug": False,
            },
        )
        rho = pcsaft_den(t, p, x, params, phase="liq")
        result = pcsaft_miac_m(t, rho, x, params, species=species)
        salt_key = next((k for k in result if "+" in k and "-" in k), None)
        if salt_key is None:
            raise KeyError(f"No salt key found in pcsaft_miac_m output keys={list(result)}.")
        gamma_calc[i] = float(result[salt_key])

    if not np.all(np.isfinite(gamma_calc)):
        raise ValueError(f"Non-finite MIAC values for dielc_rule={dielc_rule}, born_model={born_model}.")
    return gamma_calc


def _generate_combo_plots(solvent, salt):
    """Run born diagnostics numerically for one solvent/salt combination (no PNG outputs)."""
    species = _species_for_combo(solvent, salt)
    molal_exp, gamma_exp = _load_miac_data(solvent, salt)
    _axis_limits(solvent, molal_exp, gamma_exp)

    if solvent == "water" and salt.startswith("Li"):
        gamma_max = float(np.max(gamma_exp))
        if gamma_max > 5.0:
            pytest.skip(f"Skipping {salt}: filtered water data has gamma_max={gamma_max:.3f} (>5).")

    calc_count = 0
    for born_model in BORN_MODELS:
        for rule in DIELC_RULES:
            gamma_calc = _calc_rule_curve(molal_exp, species, rule, born_model)
            if gamma_calc.shape != molal_exp.shape:
                raise AssertionError(
                    f"Shape mismatch for solvent={solvent}, salt={salt}, born_model={born_model}, rule={rule}."
                )
            calc_count += 1

    if calc_count != len(BORN_MODELS) * len(DIELC_RULES):
        raise AssertionError(f"Expected {len(BORN_MODELS) * len(DIELC_RULES)} curve evaluations, got {calc_count}.")


def test_DH_born_models():
    """Generate diagnostics plots for selected solvent/salt combo."""
    solvent = _selected_solvent()
    salt = _selected_salt()
    _generate_combo_plots(solvent, salt)


@pytest.mark.parametrize("salt", WATER_SALTS)
def test_DH_born_models_water_batch(salt):
    """Run all water salts in one pytest run (numeric checks only)."""
    _generate_combo_plots("water", salt)
