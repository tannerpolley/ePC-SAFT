"""Compare dielc_rule=1 vs dielc_rule=3 for Figure 6b-style LiBr/ethanol analysis.

This script uses the same Bulow-2020 filtered data handling as figure6b and
computes ln(MIAC) contribution curves for two dielectric rules:
- dielc_rule = 1
- dielc_rule = 3

Output is a two-panel plot:
1) Absolute overlays of selected contributions for each rule.
2) Difference curves (rule1 - rule3) to make subtle deviations visible.
"""

from __future__ import annotations

import argparse
import csv
import math
import sys
from pathlib import Path
from typing import Dict, List, Tuple

import matplotlib
import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from data.epcsaft_properties import get_prop_dict
from pcsaft import pcsaft_den, pcsaft_lnfugcoef_terms, pcsaft_p

matplotlib.use("Agg")
import matplotlib.pyplot as plt


T_REF = 298.15
P_REF = 1.0e5
SALT = "LiBr"
SOLVENT = "ethanol"
CATION = "Li+"
ANION = "Br-"
SPECIES = [CATION, ANION, "Ethanol"]
MW_ETHANOL = 46.068e-3


def _read_csv_rows(path: Path) -> Tuple[List[str], List[Dict[str, str]]]:
    with path.open("r", newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        fields = [h.strip() for h in (reader.fieldnames or []) if h and h.strip()]
        rows = []
        for row in reader:
            clean: Dict[str, str] = {}
            for k, v in row.items():
                if not k:
                    continue
                ks = k.strip()
                if not ks:
                    continue
                clean[ks] = v.strip() if isinstance(v, str) else v
            rows.append(clean)
    return fields, rows


def _salt_mole_fraction_from_molality(molality: np.ndarray) -> np.ndarray:
    m = np.asarray(molality, dtype=float)
    n_solv = 1.0 / MW_ETHANOL
    denom = m + n_solv
    return np.where(denom > 0.0, m / denom, 0.0)


def _molality_for_salt_mole_fraction(x_salt: float) -> float:
    x = float(x_salt)
    if x <= 0.0:
        return 0.0
    if x >= 1.0:
        raise ValueError("Salt mole fraction target must be < 1.")
    n_solv = 1.0 / MW_ETHANOL
    return float((x * n_solv) / (1.0 - x))


def _molality_to_species_molefraction(molality: float) -> np.ndarray:
    m = max(float(molality), 0.0)
    n_solv = 1.0 / MW_ETHANOL
    n_cat = m
    n_an = m
    n_total = n_solv + n_cat + n_an
    if n_total <= 0.0:
        raise ValueError("Computed non-positive total moles while converting molality.")
    return np.asarray([n_cat / n_total, n_an / n_total, n_solv / n_total], dtype=float)


def _load_exp_data(path: Path) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    fields, rows = _read_csv_rows(path)
    lookup = {f.lower(): f for f in fields}

    m_key = None
    for candidate in ("molality", "molality (kg/mol)", "m", "m (mol/kg)"):
        if candidate in lookup:
            m_key = lookup[candidate]
            break
    if m_key is None:
        raise ValueError(f"Missing molality column in {path}.")

    source_key = lookup.get("source")
    x_key = lookup.get("mole_fraction")
    miac_key = lookup.get("miac")
    miac_m_key = lookup.get("miac_m")

    if miac_key is None and miac_m_key is None:
        raise ValueError(f"Missing MIAC columns in {path}. Need one of 'miac' or 'miac_m'.")

    molality: List[float] = []
    x_salt: List[float] = []
    ln_miac: List[float] = []

    for row in rows:
        src = str(row.get(source_key, "") if source_key else "").strip()
        if src != "Bulow 2020":
            continue

        try:
            m_val = float(row.get(m_key, ""))
        except (TypeError, ValueError):
            continue
        if not math.isfinite(m_val):
            continue

        x_val = None
        if x_key is not None:
            try:
                x_candidate = float(row.get(x_key, ""))
                if math.isfinite(x_candidate):
                    x_val = x_candidate
            except (TypeError, ValueError):
                x_val = None
        if x_val is None:
            x_val = float(_salt_mole_fraction_from_molality(np.asarray([m_val], dtype=float))[0])

        miac_val = None
        if miac_key is not None:
            try:
                y = float(row.get(miac_key, ""))
                if math.isfinite(y):
                    miac_val = y
            except (TypeError, ValueError):
                miac_val = None
        if miac_val is None and miac_m_key is not None:
            try:
                y_m = float(row.get(miac_m_key, ""))
                if math.isfinite(y_m):
                    miac_val = y_m * (1.0 + MW_ETHANOL * m_val * 2.0)
            except (TypeError, ValueError):
                miac_val = None

        if miac_val is None or (not math.isfinite(miac_val)) or miac_val <= 0.0:
            continue

        molality.append(m_val)
        x_salt.append(x_val)
        ln_miac.append(math.log(miac_val))

    if not molality:
        raise ValueError("No Bulow 2020 rows with usable MIAC data were found.")

    m_arr = np.asarray(molality, dtype=float)
    x_arr = np.asarray(x_salt, dtype=float)
    y_arr = np.asarray(ln_miac, dtype=float)
    order = np.argsort(m_arr)
    return m_arr[order], x_arr[order], y_arr[order]


def _build_params(dielc_rule: int) -> Dict[str, object]:
    x_ref = _molality_to_species_molefraction(1e-8)
    return get_prop_dict(
        "bulow_2020",
        SPECIES,
        x_ref,
        T_REF,
        user_options={"elec_model": {"rel_perm": {"rule": int(dielc_rule)}}},
    )


def _inf_dilution_state(x: np.ndarray, rho: float, params: Dict[str, object]) -> Tuple[np.ndarray, float]:
    eps = 1e-12
    x_inf = np.full_like(x, eps)
    idx_solvent = 2
    x_inf[idx_solvent] = max(1.0 - eps * (len(x) - 1), eps)
    x_inf /= np.sum(x_inf)

    p_state = pcsaft_p(T_REF, rho, x, params)
    rho_inf = pcsaft_den(T_REF, p_state, x_inf, params, phase="liq")
    return x_inf, float(rho_inf)


def _calc_ln_miac_contributions(m_grid: np.ndarray, params: Dict[str, object]) -> Dict[str, np.ndarray]:
    out = {
        "total": np.empty_like(m_grid, dtype=float),
        "born": np.empty_like(m_grid, dtype=float),
        "dh": np.empty_like(m_grid, dtype=float),
        "assoc": np.empty_like(m_grid, dtype=float),
    }

    for idx, m in enumerate(m_grid):
        m_eval = float(m) if m > 0.0 else 1e-12
        x = _molality_to_species_molefraction(m_eval)
        rho = pcsaft_den(T_REF, P_REF, x, params, phase="liq")

        terms = pcsaft_lnfugcoef_terms(T_REF, rho, x, params)
        x_inf, rho_inf = _inf_dilution_state(x, rho, params)
        terms_inf = pcsaft_lnfugcoef_terms(T_REF, rho_inf, x_inf, params)

        def mean_ionic_delta(key: str) -> float:
            a = np.asarray(terms[key], dtype=float)
            b = np.asarray(terms_inf[key], dtype=float)
            return 0.5 * ((a[0] - b[0]) + (a[1] - b[1]))

        out["assoc"][idx] = mean_ionic_delta("mu_assoc")
        out["dh"][idx] = mean_ionic_delta("mu_ion")
        out["born"][idx] = mean_ionic_delta("mu_born")
        out["total"][idx] = mean_ionic_delta("lnfugcoef_total")

    for key, arr in out.items():
        if not np.all(np.isfinite(arr)):
            raise ValueError(f"Non-finite values in contribution curve '{key}'.")

    return out


def run_analysis(
    data_path: Path,
    output_path: Path,
    x_min: float,
    x_max: float,
    grid_points: int,
    max_molality: float | None,
) -> Path:
    m_exp, x_exp, y_exp = _load_exp_data(data_path)

    if max_molality is None:
        m_upper = max(float(np.max(m_exp)), _molality_for_salt_mole_fraction(float(x_max)))
    else:
        m_upper = float(max_molality)
    if m_upper <= 0.0:
        m_upper = float(np.max(m_exp[m_exp >= 0.0]))

    m_grid = np.linspace(0.0, m_upper, int(grid_points))
    x_grid = _salt_mole_fraction_from_molality(m_grid)

    curves_rule1 = _calc_ln_miac_contributions(m_grid, _build_params(1))
    curves_rule3 = _calc_ln_miac_contributions(m_grid, _build_params(3))

    delta = {k: curves_rule1[k] - curves_rule3[k] for k in curves_rule1.keys()}

    fig, axes = plt.subplots(1, 2, figsize=(13.2, 5.4), sharex=True)
    fig.patch.set_facecolor("white")
    ax_abs, ax_delta = axes
    ax_abs.set_facecolor("white")
    ax_delta.set_facecolor("white")

    ax_abs.scatter(
        x_exp,
        y_exp,
        color="black",
        marker="o",
        s=30,
        facecolors="black",
        linewidths=0.7,
        label="Experimental data (Bulow 2020)",
        zorder=8,
    )

    ax_abs.plot(x_grid, curves_rule1["total"], color="green", linewidth=2.0, linestyle="-", label="Total, rule 1")
    ax_abs.plot(x_grid, curves_rule3["total"], color="green", linewidth=2.0, linestyle="--", label="Total, rule 3")
    ax_abs.plot(x_grid, curves_rule1["dh"], color="tab:blue", linewidth=1.8, linestyle="-", label="DH, rule 1")
    ax_abs.plot(x_grid, curves_rule3["dh"], color="tab:blue", linewidth=1.8, linestyle="--", label="DH, rule 3")
    ax_abs.plot(x_grid, curves_rule1["assoc"], color="gray", linewidth=1.8, linestyle="-", label="Assoc, rule 1")
    ax_abs.plot(x_grid, curves_rule3["assoc"], color="gray", linewidth=1.8, linestyle="--", label="Assoc, rule 3")
    ax_abs.plot(x_grid, curves_rule1["born"], color="orange", linewidth=1.6, linestyle="-", alpha=0.85, label="Born, rule 1")
    ax_abs.plot(x_grid, curves_rule3["born"], color="orange", linewidth=1.6, linestyle="--", alpha=0.85, label="Born, rule 3")

    ax_abs.set_xlim(float(x_min), float(x_max))
    ax_abs.set_xlabel(r"salt mole fraction, $x_{salt}$")
    ax_abs.set_ylabel(r"$\ln(\gamma_{\pm}^{*})$")
    ax_abs.set_title("Absolute curves")
    ax_abs.grid(True, alpha=0.3, color="0.7")
    ax_abs.legend(fontsize=8, ncol=2)

    ax_delta.plot(x_grid, delta["total"], color="green", linewidth=1.9, linestyle="-", label="delta total (1 - 3)")
    ax_delta.plot(x_grid, delta["dh"], color="tab:blue", linewidth=1.8, linestyle="-", label="delta DH (1 - 3)")
    ax_delta.plot(x_grid, delta["assoc"], color="gray", linewidth=1.8, linestyle="-", label="delta assoc (1 - 3)")
    ax_delta.plot(x_grid, delta["born"], color="orange", linewidth=1.6, linestyle="-", label="delta Born (1 - 3)")
    ax_delta.axhline(0.0, color="black", linewidth=1.0, linestyle="--")

    max_abs_delta = max(float(np.max(np.abs(v))) for v in delta.values())
    pad = max(0.05, 0.15 * max_abs_delta)
    y_lim = max_abs_delta + pad
    ax_delta.set_ylim(-y_lim, y_lim)
    ax_delta.set_xlim(float(x_min), float(x_max))
    ax_delta.set_xlabel(r"salt mole fraction, $x_{salt}$")
    ax_delta.set_ylabel("delta ln contribution")
    ax_delta.set_title("Rule-1 minus Rule-3")
    ax_delta.grid(True, alpha=0.3, color="0.7")
    ax_delta.legend(fontsize=8)

    for ax in axes:
        ax.tick_params(colors="black")
        for spine in ax.spines.values():
            spine.set_color("black")
            spine.set_linewidth(1.0)

    fig.suptitle("LiBr in ethanol at 298.15 K and 1 bar: dielectric-rule comparison", fontsize=12)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.95))
    fig.savefig(output_path, dpi=220)
    plt.close(fig)

    if not output_path.exists():
        raise FileNotFoundError(f"Expected plot was not written: {output_path}")

    print(f"Loaded Bulow-2020 rows: {len(m_exp)} from {data_path}")
    print(f"Molality grid points: {len(m_grid)} (0 to {m_upper:.6g} mol/kg)")
    for key in ("total", "dh", "assoc", "born"):
        print(
            f"max |delta {key}| = {float(np.max(np.abs(delta[key]))):.6e} "
            f"(rule1 - rule3)"
        )
    print(f"Wrote: {output_path}")
    return output_path


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Figure 6b dielc_rule comparison for LiBr/ethanol")
    parser.add_argument(
        "--data",
        type=Path,
        default=REPO_ROOT / "data" / "MIAC" / "ethanol" / "ethanol-LiBr.csv",
        help="Input CSV with columns including molality, source, and miac/miac_m.",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=REPO_ROOT / "scripts" / "2020_model_analysis" / "output" / "figure6b_libr_ethanol_dielc_rule_comparison.png",
        help="Output PNG path.",
    )
    parser.add_argument("--x-min", type=float, default=0.0)
    parser.add_argument("--x-max", type=float, default=0.2)
    parser.add_argument("--grid-points", type=int, default=1201)
    parser.add_argument(
        "--max-molality",
        type=float,
        default=None,
        help="Optional upper molality bound. Default uses max(data molality, x-max mapping).",
    )
    return parser.parse_args()


def main() -> None:
    args = _parse_args()
    run_analysis(
        data_path=Path(args.data),
        output_path=Path(args.out),
        x_min=float(args.x_min),
        x_max=float(args.x_max),
        grid_points=int(args.grid_points),
        max_molality=None if args.max_molality is None else float(args.max_molality),
    )


if __name__ == "__main__":
    main()
