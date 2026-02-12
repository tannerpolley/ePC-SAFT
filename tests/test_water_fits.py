# -*- coding: utf-8 -*-
"""Water MIAC fit plots for chloride and bromide salts using dielc_rule=1."""

import csv
from pathlib import Path

import matplotlib
import numpy as np

from pcsaft import pcsaft_den, pcsaft_miac_m
from data.epcsaft_properties import get_prop_dict, molality_to_molefraction

matplotlib.use("Agg")
import matplotlib.pyplot as plt


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


def _style_for_salt(salt):
    """Set marker/color style by alkali cation identity."""
    if salt.startswith("Li"):
        return "^", "orange"
    if salt.startswith("Na"):
        return "s", "green"
    if salt.startswith("K"):
        return "o", "gray"
    raise ValueError(f"Unsupported salt style: {salt}")


def _load_miac_data(salt, m_min=0.0, m_max=6.0):
    """Load molality and gamma from CSV for one salt in water."""
    root = Path(__file__).resolve().parents[1]
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


def _plot_group(salts, filename):
    """Create one MIAC fit plot for a salt group."""
    out_dir = Path(__file__).resolve().parent / "fit_plots"
    out_dir.mkdir(parents=True, exist_ok=True)
    plot_path = out_dir / filename

    fig, ax = plt.subplots(figsize=(7.5, 5.0))
    rule = 1
    for salt in salts:
        molal, gamma_exp = _load_miac_data(salt)
        species = _species_for_salt(salt)
        gamma_calc = _calc_miac_curve(molal, species, rule=rule)
        marker, color = _style_for_salt(salt)

        ax.scatter(molal, gamma_exp, marker=marker, color=color, s=34, alpha=0.9, label=f"{salt} data")
        ax.plot(molal, gamma_calc, color=color, linewidth=1.8, label=f"{salt} fit (rule {rule})")
    title = f"Water MIAC fits at 298.15 K (dielc rule {rule})",
    ax.set_xlim(0.0, 6.0)
    ax.set_ylim(0.0, 4.0)
    ax.set_yticks(np.arange(0.0, 4.0 + 0.5, 0.5))
    ax.set_xlabel("molality, m / mol kg$^{-1}$")
    ax.set_ylabel("mean ionic activity coefficient, $\\gamma_{\\pm}^{m}$")
    ax.set_title(title)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8, ncol=2)
    fig.tight_layout()
    fig.savefig(plot_path, dpi=220)
    plt.close(fig)

    if not plot_path.exists():
        raise FileNotFoundError(f"Expected plot was not written: {plot_path}")
    return plot_path


def test_water_fits():
    """Generate chloride and bromide water-fit plots with rule 1."""
    path_cl = _plot_group(
        salts=["LiCl", "NaCl", "KCl"],
        filename="test_water_fits_cl_salts.png",
    )
    path_br = _plot_group(
        salts=["LiBr", "NaBr", "KBr"],
        filename="test_water_fits_br_salts.png",
    )

    if not path_cl.exists() or not path_br.exists():
        raise FileNotFoundError("Expected water fit plots were not written.")
