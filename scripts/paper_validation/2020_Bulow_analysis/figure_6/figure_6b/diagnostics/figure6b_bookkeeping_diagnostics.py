"""Diagnose exact Figure 6b bookkeeping, including the compressibility correction."""

from __future__ import annotations

import argparse
import csv
import math
import sys
from pathlib import Path
from typing import Dict, List

import matplotlib
import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[6]
SCRIPT_DIR = Path(__file__).resolve().parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts._env import require_pcsaft_install

require_pcsaft_install()
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from figure6b_libr_ethanol_contributions import (
    P_REF,
    T_REF,
    _build_params,
    _inf_dilution_state,
    _load_exp_data,
    _molality_to_species_molefraction,
    _salt_mole_fraction_from_molality,
)
from pcsaft import pcsaft_Z, pcsaft_den, pcsaft_lnfugcoef_terms

matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUTPUT_ROOT = REPO_ROOT / "scripts" / "paper_validation" / "2020_Bulow_analysis" / "figure_6" / "figure_6b" / "diagnostics" / "output"
OUTPUT_DATA_DIR = OUTPUT_ROOT / "data"
OUTPUT_PLOTS_DIR = OUTPUT_ROOT / "plots"


def _mean_ionic_delta(terms: Dict[str, np.ndarray], terms_inf: Dict[str, np.ndarray], key: str) -> float:
    a = np.asarray(terms[key], dtype=float)
    b = np.asarray(terms_inf[key], dtype=float)
    return float(0.5 * ((a[0] - b[0]) + (a[1] - b[1])))


def run_analysis(
    miac_data_path: Path,
    output_csv: Path,
    output_plot: Path,
    grid_points: int,
) -> Dict[str, float]:
    m_exp, _, _ = _load_exp_data(miac_data_path)
    m_upper = float(np.max(m_exp))
    m_grid = np.linspace(0.0, m_upper, int(grid_points))
    x_grid = _salt_mole_fraction_from_molality(m_grid)
    params = _build_params(user_options={})

    rows: List[Dict[str, float]] = []
    max_closure = 0.0
    max_lnfug_closure = 0.0
    max_mu_closure = 0.0

    for m, x_salt in zip(m_grid, x_grid):
        m_eval = float(m) if m > 0.0 else 1e-12
        x = _molality_to_species_molefraction(m_eval)
        rho = float(pcsaft_den(T_REF, P_REF, x, params, phase="liq"))
        Z = float(pcsaft_Z(T_REF, rho, x, params))
        terms = pcsaft_lnfugcoef_terms(T_REF, rho, x, params)

        x_inf, rho_inf = _inf_dilution_state(x, rho, params)
        Z_inf = float(pcsaft_Z(T_REF, rho_inf, x_inf, params))
        terms_inf = pcsaft_lnfugcoef_terms(T_REF, rho_inf, x_inf, params)

        hc = _mean_ionic_delta(terms, terms_inf, "mu_hc")
        disp = _mean_ionic_delta(terms, terms_inf, "mu_disp")
        assoc = _mean_ionic_delta(terms, terms_inf, "mu_assoc")
        dh = _mean_ionic_delta(terms, terms_inf, "mu_ion")
        born = _mean_ionic_delta(terms, terms_inf, "mu_born")
        lnfug_hc = _mean_ionic_delta(terms, terms_inf, "lnfugcoef_hc")
        lnfug_disp = _mean_ionic_delta(terms, terms_inf, "lnfugcoef_disp")
        lnfug_assoc = _mean_ionic_delta(terms, terms_inf, "lnfugcoef_assoc")
        lnfug_dh = _mean_ionic_delta(terms, terms_inf, "lnfugcoef_ion")
        lnfug_born = _mean_ionic_delta(terms, terms_inf, "lnfugcoef_born")
        mu_total = _mean_ionic_delta(terms, terms_inf, "mu_total")
        lnfug_total = _mean_ionic_delta(terms, terms_inf, "lnfugcoef_total")
        mu_sum = hc + disp + assoc + dh + born
        lnfug_sum = lnfug_hc + lnfug_disp + lnfug_assoc + lnfug_dh + lnfug_born
        zcorr = float(math.log(Z_inf / Z))
        closure = lnfug_total - lnfug_sum
        mu_reference_gap = lnfug_total - mu_sum
        mu_closure = mu_total - mu_sum
        max_closure = max(max_closure, abs(closure))
        max_lnfug_closure = max(max_lnfug_closure, abs(mu_reference_gap))
        max_mu_closure = max(max_mu_closure, abs(mu_closure))

        rows.append(
            {
                "molality": float(m_eval),
                "x_salt": float(x_salt),
                "rho": rho,
                "rho_inf": float(rho_inf),
                "Z": Z,
                "Z_inf": Z_inf,
                "hc": hc,
                "disp": disp,
                "assoc": assoc,
                "dh": dh,
                "born": born,
                "lnfug_hc": lnfug_hc,
                "lnfug_disp": lnfug_disp,
                "lnfug_assoc": lnfug_assoc,
                "lnfug_dh": lnfug_dh,
                "lnfug_born": lnfug_born,
                "mu_sum": mu_sum,
                "lnfug_sum": lnfug_sum,
                "mu_total": mu_total,
                "lnfug_total": lnfug_total,
                "z_correction": zcorr,
                "closure_lnfug_minus_lnfug_sum": closure,
                "closure_lnfug_minus_mu_sum": mu_reference_gap,
                "closure_mu_total_minus_mu_sum": mu_closure,
            }
        )

    output_csv.parent.mkdir(parents=True, exist_ok=True)
    with output_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "molality",
                "x_salt",
                "rho",
                "rho_inf",
                "Z",
                "Z_inf",
                "hc",
                "disp",
                "assoc",
                "dh",
                "born",
                "lnfug_hc",
                "lnfug_disp",
                "lnfug_assoc",
                "lnfug_dh",
                "lnfug_born",
                "mu_sum",
                "lnfug_sum",
                "mu_total",
                "lnfug_total",
                "z_correction",
                "closure_lnfug_minus_lnfug_sum",
                "closure_lnfug_minus_mu_sum",
                "closure_mu_total_minus_mu_sum",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)

    fig, ax = plt.subplots(figsize=(8.4, 5.6))
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")
    ax.plot([r["x_salt"] for r in rows], [r["lnfug_total"] for r in rows], linewidth=2.1, color="green", label=r"$\ln(\gamma_{\pm}^{*})$ total")
    ax.plot([r["x_salt"] for r in rows], [r["lnfug_sum"] for r in rows], linewidth=1.9, color="tab:blue", linestyle="--", label=r"sum of $\ln\varphi^\alpha$ contributions")
    ax.plot([r["x_salt"] for r in rows], [r["mu_sum"] for r in rows], linewidth=1.6, color="black", linestyle=":", label=r"legacy sum of $\mu^\alpha$ contributions")
    ax.plot([r["x_salt"] for r in rows], [r["z_correction"] for r in rows], linewidth=1.6, color="tab:red", linestyle="-.", label=r"$\ln(Z_{\infty}/Z)$")
    ax.set_xlabel(r"salt mole fraction, $x_{salt}$")
    ax.set_ylabel(r"Contribution to $\ln(\gamma_{\pm}^{*})$")
    ax.set_title(r"Figure 6b bookkeeping: total vs summed $\ln\varphi^\alpha$ contributions")
    ax.grid(True, alpha=0.3, color="0.7")
    ax.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(output_plot, dpi=220)
    plt.close(fig)

    print(f"Wrote bookkeeping CSV: {output_csv}")
    print(f"Wrote bookkeeping plot: {output_plot}")
    print(f"Max |lnfug_total - lnfug_sum| = {max_closure:.6e}")
    print(f"Max |lnfug_total - mu_sum| = {max_lnfug_closure:.6e}")
    print(f"Max |mu_total - mu_sum| = {max_mu_closure:.6e}")
    return {
        "max_closure_lnfug_minus_lnfug_sum": max_closure,
        "max_closure_lnfug_minus_mu_sum": max_lnfug_closure,
        "max_closure_mu_total_minus_mu_sum": max_mu_closure,
    }


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Diagnose Figure 6b bookkeeping and -ln(Z) correction")
    parser.add_argument(
        "--miac-data",
        type=Path,
        default=REPO_ROOT / "data" / "MIAC" / "ethanol" / "ethanol-LiBr.csv",
    )
    parser.add_argument(
        "--out-csv",
        type=Path,
        default=OUTPUT_DATA_DIR / "figure6b_bookkeeping.csv",
    )
    parser.add_argument(
        "--out-plot",
        type=Path,
        default=OUTPUT_PLOTS_DIR / "figure6b_bookkeeping.png",
    )
    parser.add_argument("--grid-points", type=int, default=1201)
    return parser.parse_args()


def main() -> None:
    args = _parse_args()
    run_analysis(
        miac_data_path=Path(args.miac_data),
        output_csv=Path(args.out_csv),
        output_plot=Path(args.out_plot),
        grid_points=int(args.grid_points),
    )


if __name__ == "__main__":
    main()
