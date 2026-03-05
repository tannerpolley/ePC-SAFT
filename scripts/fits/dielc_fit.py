# -*- coding: utf-8 -*-
"""Dielectric-constant fit visualization for salts in water (rules 1 and 4)."""

import csv
from pathlib import Path

import matplotlib
import numpy as np

from data.epcsaft_properties import get_prop_dict
from pcsaft import pcsaft_dielc_eval

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def _load_dielc_data():
    """Load salt dielectric data and return arrays for salt, x_ion, dielc_exp."""
    data_path = Path(__file__).resolve().parents[2] / "data" / "dielc" / "dielc_salts_in_water.csv"
    if not data_path.exists():
        raise FileNotFoundError(f"Data file not found: {data_path}")

    salts, x_ion, dielc_exp = [], [], []
    with data_path.open("r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            norm = {k.strip().lower(): v for k, v in row.items()}
            salts.append(norm["salt"].strip())
            x_ion.append(float(norm["x_ion"]))
            dielc_exp.append(float(norm["dielc"]))

    if not salts:
        raise ValueError("No rows found in dielc_salts_in_water.csv.")

    x_ion = np.asarray(x_ion, dtype=float)
    dielc_exp = np.asarray(dielc_exp, dtype=float)
    if not np.all(np.isfinite(x_ion)) or not np.all(np.isfinite(dielc_exp)):
        raise ValueError("Non-finite values found in dielectric dataset.")
    return np.asarray(salts), x_ion, dielc_exp


def _species_for_salt(salt):
    """Map salt label to species list for get_prop_dict."""
    mapping = {
        "NaCl": ["Na+", "Cl-", "H2O"],
        "NaBr": ["Na+", "Br-", "H2O"],
        "LiCl": ["Li+", "Cl-", "H2O"],
    }
    if salt not in mapping:
        raise ValueError(f"Unsupported salt in dielectric dataset: {salt}")
    return mapping[salt]


def _calc_dielc_curve(x_ion_grid, species, rule, t=298.15):
    """Calculate dielectric constant for a 1:1 electrolyte over x_ion grid."""
    dielc = np.empty_like(x_ion_grid, dtype=float)
    for i, x_ion in enumerate(x_ion_grid):
        x = np.asarray([0.5 * x_ion, 0.5 * x_ion, 1.0 - x_ion], dtype=float)
        params = get_prop_dict("bulow_2020", species, x, t, user_options={"elec_model": {"rel_perm": {"rule": int(rule)}}})
        dielc[i] = float(pcsaft_dielc_eval(x, params)[0])

    if not np.all(np.isfinite(dielc)):
        raise ValueError(f"Non-finite dielectric values for dielc_rule={rule}, species={species}.")
    return dielc


def test_dielc_fit():
    """Plot dielc data and model lines (rule 1 solid, rule 4 dashed)."""
    salts, x_ion_data, dielc_exp = _load_dielc_data()

    x_ion_min = float(np.min(x_ion_data))
    x_ion_max = float(np.max(x_ion_data))
    x_ion_grid = np.linspace(x_ion_min, x_ion_max, 200)

    species_ref = _species_for_salt("NaCl")
    dielc_rule1 = _calc_dielc_curve(x_ion_grid, species_ref, rule=1)
    dielc_rule4 = _calc_dielc_curve(x_ion_grid, species_ref, rule=4)

    out_dir = Path(__file__).resolve().parents[2] / "data" / "dielc" / "plot_fits"
    out_dir.mkdir(parents=True, exist_ok=True)
    plot_path = out_dir / "test_dielc_fit_salts_in_water_fit.png"

    fig, ax = plt.subplots(figsize=(7.5, 5.0))
    marker_map = {"NaCl": "o", "NaBr": "s", "LiCl": "^"}
    color_map = {"NaCl": "tab:blue", "NaBr": "tab:orange", "LiCl": "tab:green"}
    for salt in ["NaCl", "NaBr", "LiCl"]:
        mask = salts == salt
        if np.any(mask):
            ax.scatter(
                x_ion_data[mask],
                dielc_exp[mask],
                s=34,
                marker=marker_map[salt],
                color=color_map[salt],
                label=f"{salt} data",
                alpha=0.9,
            )

    ax.plot(x_ion_grid, dielc_rule1, color="black", linestyle="--", linewidth=1.8, label="rule 1")
    ax.plot(x_ion_grid, dielc_rule4, color="black", linestyle="--", linewidth=2.2, label="rule 4")

    ax.set_xlabel(r"$x_{\mathrm{ion}} = x_{+} + x_{-}$")
    ax.set_ylabel(r"dielectric constant, $\varepsilon_r$")
    ax.set_title("Salts in water dielectric constant fit (298.15 K)")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(plot_path, dpi=220)
    plt.close(fig)

    if not plot_path.exists():
        raise FileNotFoundError(f"Expected plot was not written: {plot_path}")
    if dielc_rule1.size == 0 or dielc_rule4.size == 0:
        raise ValueError("Model line arrays are empty.")
