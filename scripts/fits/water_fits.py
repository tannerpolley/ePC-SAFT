# -*- coding: utf-8 -*-
"""Water MIAC fit plots for chloride and bromide salts using dielc_rule=1."""

import csv
from pathlib import Path

import numpy as np

from pcsaft import pcsaft_den, pcsaft_miac_m
from data.epcsaft_properties import get_prop_dict, molality_to_molefraction


def _species_for_salt(salt):
    """Map each salt to a cation/anion/water species triplet."""
    mapping = {
        "LiCl": ["Li+", "Cl-", "H2O-2B-Li"],
        "NaCl": ["Na+", "Cl-", "H2O-2B-NaCl"],
        "KCl": ["K+", "Cl-", "H2O-2B-NaCl"],
        "LiBr": ["Li+", "Br-", "H2O-2B-Li"],
        "NaBr": ["Na+", "Br-", "H2O-2B-NaCl"],
        "KBr": ["K+", "Br-", "H2O-2B-NaCl"],
    }
    if salt not in mapping:
        raise ValueError(f"Unsupported salt: {salt}")
    return mapping[salt]


def _load_miac_data(salt, m_min=0.0, m_max=6.0):
    """Load molality and gamma from CSV for one salt in water."""
    root = Path(__file__).resolve().parents[2]
    candidates = [
        root / "data" / "MIAC_m" / "water" / f"{salt}.csv",
        root / "data" / "MIAC_m" / "csv" / f"{salt}.csv",
    ]
    data_path = next((path for path in candidates if path.exists()), None)
    if data_path is None:
        raise FileNotFoundError(f"Data file not found in candidates: {candidates}")

    molal, gamma_exp = [], []
    with data_path.open("r", newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            m = float(row["m"])
            g = float(row["gamma"])
            if m_min <= m <= m_max:
                molal.append(m)
                gamma_exp.append(g)

    if not molal:
        raise ValueError(f"No rows in range {m_min}<=m<={m_max} for {salt}.")
    molal = np.asarray(molal, dtype=float)
    gamma_exp = np.asarray(gamma_exp, dtype=float)
    order = np.argsort(molal)
    return molal[order], gamma_exp[order]


def _calc_miac_curve(molal, species, rule=1, t=298.15, p=101325.0):
    """Calculate model MIAC curve for one salt/species set."""
    gamma_calc = np.empty_like(molal, dtype=float)
    for i, m_salt in enumerate(molal):
        x = molality_to_molefraction(float(m_salt), species=species)
        params = get_prop_dict(species, x, t, user_options={"dielc_rule": int(rule), "born_model": 0, "debug": False})
        rho = pcsaft_den(t, p, x, params, phase="liq")
        result = pcsaft_miac_m(t, rho, x, params, species=species)
        salt_key = next(k for k in result if "+" in k and "-" in k)
        gamma_calc[i] = float(result[salt_key])

    if not np.all(np.isfinite(gamma_calc)):
        raise ValueError(f"Non-finite MIAC values for species={species}.")
    return gamma_calc


def test_water_fits():
    """Run water MIAC fit checks without writing plot files."""
    rule = 1
    salts = ["LiCl", "NaCl", "KCl", "LiBr", "NaBr", "KBr"]
    for salt in salts:
        molal, gamma_exp = _load_miac_data(salt)
        species = _species_for_salt(salt)
        gamma_calc = _calc_miac_curve(molal, species, rule=rule)
        if gamma_calc.shape != gamma_exp.shape:
            raise AssertionError(f"Shape mismatch for salt {salt}.")
        if not np.all(np.isfinite(gamma_calc)):
            raise ValueError(f"Non-finite MIAC results for salt {salt}.")
