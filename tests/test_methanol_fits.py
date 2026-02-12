# -*- coding: utf-8 -*-
"""Methanol fit diagnostics for molality-scale MIAC against experimental data."""

import csv
from pathlib import Path

import matplotlib
import numpy as np
import pytest

from pcsaft import pcsaft_den, pcsaft_miac_m
from data.epcsaft_properties import get_prop_dict, molality_to_molefraction

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def _load_nabr_methanol_data():
    """Load NaBr-in-methanol data for the 100% methanol set in 0-2 mol/kg."""
    root = Path(__file__).resolve().parents[1]
    data_paths = [
        root / "data" / "MIAC_m" / "methanol" / "methanol-NaBr.csv",
        root / "data" / "MIAC_m" / "water_methanol" / "water-methanol-NaBr.csv",
        root / "data" / "MIAC_m" / "water_methanol" / "NaBr.csv",
    ]
    data_path = next((p for p in data_paths if p.exists()), None)
    if data_path is None:
        raise FileNotFoundError(f"Data file not found in candidates: {data_paths}")

    molal = []
    gamma_exp = []
    with data_path.open("r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        fields = set(reader.fieldnames or [])
        y_key = "gamma" if "gamma" in fields else ("miac_m" if "miac_m" in fields else None)
        if {"molality (kg/mol)"}.issubset(fields) and y_key is not None:
            for row in reader:
                m = float(row["molality (kg/mol)"])
                g = float(row[y_key])
                if 0.0 <= m <= 2.0:
                    molal.append(m)
                    gamma_exp.append(g)
        elif {"x_methanol", "molality (kg/mol)"}.issubset(fields) and y_key is not None:
            for row in reader:
                x_meoh = float(row["x_methanol"])
                m = float(row["molality (kg/mol)"])
                g = float(row[y_key])
                if abs(x_meoh - 1.0) <= 1e-12 and 0.0 <= m <= 2.0:
                    molal.append(m)
                    gamma_exp.append(g)
        elif {"methanol_pct", "molal", "gamma_exp"}.issubset(fields):
            for row in reader:
                pct = float(row["methanol_pct"])
                m = float(row["molal"])
                g = float(row["gamma_exp"])
                if pct == 100.0 and 0.0 <= m <= 2.0:
                    molal.append(m)
                    gamma_exp.append(g)
        else:
            raise ValueError(f"Unsupported methanol NaBr schema in {data_path}: {sorted(fields)}")

    if not molal:
        raise ValueError(f"No rows matched 100% methanol and 0<=molality<=2 in {data_path}.")

    molal = np.asarray(molal, dtype=float)
    gamma_exp = np.asarray(gamma_exp, dtype=float)
    order = np.argsort(molal)
    return molal[order], gamma_exp[order]


def _calc_rule_curve(molal, rule, born_model=1, dh_model=1):
    """Calculate model-predicted mean ionic activity coefficients for one dielectric rule."""
    t = 298.15
    p = 101325.0
    species = ["Na+", "Br-", "Methanol"]

    gamma_calc = np.empty_like(molal, dtype=float)
    for i, m_salt in enumerate(molal):
        x = molality_to_molefraction(float(m_salt), species=species)
        params = get_prop_dict(
            species,
            x,
            t,
            user_options={"born_model": born_model, "DH_model": dh_model, "debug": False, "dielc_rule": rule},
        )
        rho = pcsaft_den(t, p, x, params, phase="liq")
        gamma_calc[i] = pcsaft_miac_m(t, rho, x, params, species=species)["Na+Br-"]

    if not np.all(np.isfinite(gamma_calc)):
        raise ValueError(f"Non-finite MIAC values for dielc_rule={rule}.")
    return gamma_calc


def _calc_metrics(gamma_calc, gamma_exp):
    """Return AARD (%) and RMSE for one model curve."""
    if np.any(gamma_exp == 0.0):
        raise ValueError("Experimental gamma contains zero; AARD is undefined.")
    aard = float(np.mean(np.abs((gamma_calc - gamma_exp) / gamma_exp)) * 100.0)
    rmse = float(np.sqrt(np.mean((gamma_calc - gamma_exp) ** 2)))
    return aard, rmse


def test_methanol_miac_m_fit():
    """Fit NaBr-in-methanol MIAC over 0-2 mol/kg for dielc rules 0-6 and save diagnostics."""
    molal_exp, gamma_exp = _load_nabr_methanol_data()

    out_dir = Path(__file__).resolve().parent / "fit_plots"
    out_dir.mkdir(parents=True, exist_ok=True)
    plot_path = out_dir / "test_methanol_fits_nabr_methanol_100pct.png"
    metrics_path = out_dir / "test_methanol_fits_nabr_methanol_100pct_metrics.csv"

    rules = [0, 1, 2, 3, 4, 5, 6]
    curves = {}
    metrics_rows = []

    for rule in rules:
        gamma_calc = _calc_rule_curve(molal_exp, rule)
        aard, rmse = _calc_metrics(gamma_calc, gamma_exp)
        curves[rule] = gamma_calc
        metrics_rows.append((rule, aard, rmse))

    print("\nNaBr in methanol (100%) fit metrics over 0<=m<=2 mol/kg")
    print("rule | AARD(%) | RMSE")
    for rule, aard, rmse in metrics_rows:
        print(f"{rule:>4d} | {aard:>7.3f} | {rmse:>7.5f}")

    with metrics_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["dielc_rule", "AARD_percent", "RMSE"])
        for row in metrics_rows:
            writer.writerow(row)

    fig, ax = plt.subplots(figsize=(7.5, 5.0))
    ax.scatter(molal_exp, gamma_exp, color="black", marker="o", s=35, label="exp data")
    for rule, aard, rmse in metrics_rows:
        ax.plot(molal_exp, curves[rule], linewidth=1.8, label=f"rule {rule}: AARD={aard:.2f}% RMSE={rmse:.4f}")

    ax.set_xlabel("molality, m / mol kg$^{-1}$")
    ax.set_ylabel("mean ionic activity coefficient, $\\gamma_{\\pm}^{m}$")
    ax.set_title("NaBr in methanol (100%) at 298.15 K")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(plot_path, dpi=220)
    plt.close(fig)

    if not plot_path.exists() or not metrics_path.exists():
        raise FileNotFoundError("Expected fit artifacts were not written.")


def test_methanol_dh_and_born_model_variants():
    """Exercise DH/Born model switches and verify DH_model 0/1 alias behavior."""
    molal_exp, _ = _load_nabr_methanol_data()
    m_eval = np.asarray([float(molal_exp[-1])], dtype=float)

    g_dh0_b1 = _calc_rule_curve(m_eval, rule=1, born_model=1, dh_model=0)[0]
    g_dh1_b1 = _calc_rule_curve(m_eval, rule=1, born_model=1, dh_model=1)[0]
    assert np.isfinite(g_dh0_b1)
    assert np.isfinite(g_dh1_b1)
    assert abs(g_dh0_b1 - g_dh1_b1) < 1e-12

    g_b2 = _calc_rule_curve(m_eval, rule=1, born_model=2, dh_model=1)[0]
    g_b3 = _calc_rule_curve(m_eval, rule=1, born_model=3, dh_model=1)[0]
    g_b4 = _calc_rule_curve(m_eval, rule=1, born_model=4, dh_model=1)[0]
    g_b5 = _calc_rule_curve(m_eval, rule=1, born_model=5, dh_model=1)[0]
    assert np.isfinite(g_b2)
    assert np.isfinite(g_b3)
    assert np.isfinite(g_b4)
    assert np.isfinite(g_b5)
    assert abs(g_b2 - g_b3) > 1e-10
    assert abs(g_b2 - g_b4) > 1e-10
    assert abs(g_b5 - g_b3) > 1e-10 or abs(g_b5 - g_b4) > 1e-10

    with pytest.raises(Exception, match="DH_model=2"):
        _calc_rule_curve(m_eval, rule=1, born_model=1, dh_model=2)
