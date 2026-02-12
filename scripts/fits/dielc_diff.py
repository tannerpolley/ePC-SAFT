# -*- coding: utf-8 -*-
"""Dielectric-differential visualization for salts in water (rules 1-6)."""

import csv
from pathlib import Path

import matplotlib
import numpy as np

from data.epcsaft_properties import get_prop_dict
from pcsaft import pcsaft_dielc_eval

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def _load_xion_range():
    """Load x_ion data and return min/max range from the salt dielectric CSV."""
    data_path = Path(__file__).resolve().parents[2] / "data" / "dielc" / "dielc_salts_in_water.csv"
    if not data_path.exists():
        raise FileNotFoundError(f"Data file not found: {data_path}")

    x_ion = []
    with data_path.open("r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            norm = {k.strip().lower(): v for k, v in row.items()}
            x_ion.append(float(norm["x_ion"]))

    if not x_ion:
        raise ValueError("No x_ion values found in dielc_salts_in_water.csv.")
    x_ion = np.asarray(x_ion, dtype=float)
    if not np.all(np.isfinite(x_ion)):
        raise ValueError("Non-finite x_ion values in dielectric dataset.")
    return float(np.min(x_ion)), float(np.max(x_ion))


def _calc_dielc_diff_curve(x_ion_grid, rule, t=298.15):
    """
    Compute effective d(eps_r)/d(x_ion) along x=[x_ion/2, x_ion/2, 1-x_ion]
    using the C++ dielectric evaluator.
    """
    species = ["Na+", "Cl-", "H2O-2B-NaCl"]
    diff_curve = np.empty_like(x_ion_grid, dtype=float)

    for i, x_ion in enumerate(x_ion_grid):
        x = np.asarray([0.5 * x_ion, 0.5 * x_ion, 1.0 - x_ion], dtype=float)
        params = get_prop_dict(species, x, t, user_options={"dielc_rule": int(rule)})
        deps_dx = np.asarray(pcsaft_dielc_eval(x, params)[1], dtype=float)
        # Chain-rule projection onto x_ion path:
        # d/dx_ion = 0.5*d/dx_cation + 0.5*d/dx_anion - d/dx_solvent
        diff_curve[i] = 0.5 * deps_dx[0] + 0.5 * deps_dx[1] - deps_dx[2]

    if not np.all(np.isfinite(diff_curve)):
        raise ValueError(f"Non-finite dielc_diff values for rule={rule}.")
    return diff_curve


def test_dielc_diff():
    """Plot effective dielectric differential curves for rules 1-6 over the x_ion range."""
    x_min, x_max = _load_xion_range()
    x_ion_grid = np.linspace(x_min, x_max, 200)

    rules = [1, 2, 3, 4, 5, 6]
    style_map = {1: "--", 2: "--", 3: "--", 4: "--", 5: "--", 6: "--"}
    curves = {rule: _calc_dielc_diff_curve(x_ion_grid, rule=rule) for rule in rules}

    out_dir = Path(__file__).resolve().parents[2] / "data" / "dielc" / "plot_fits"
    out_dir.mkdir(parents=True, exist_ok=True)
    plot_path = out_dir / "test_dielc_diff_salts_in_water_fit.png"

    fig, ax = plt.subplots(figsize=(7.5, 5.0))
    for rule in rules:
        ax.plot(
            x_ion_grid,
            curves[rule],
            linestyle=style_map[rule],
            linewidth=1.8,
            label=f"rule {rule}",
        )

    ax.set_xlabel(r"$x_{\mathrm{ion}} = x_{+} + x_{-}$")
    ax.set_ylabel(r"effective $d\varepsilon_r/dx_{\mathrm{ion}}$")
    ax.set_title("Dielectric differential over salt-in-water x_ion range (298.15 K)")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(plot_path, dpi=220)
    plt.close(fig)

    if not plot_path.exists():
        raise FileNotFoundError(f"Expected plot was not written: {plot_path}")
    for rule in rules:
        if curves[rule].size == 0:
            raise ValueError(f"Empty dielc_diff curve for rule {rule}.")
