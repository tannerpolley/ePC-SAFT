"""Figure 6b-style DH contribution comparison for LiBr in ethanol with d_ion variants."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Dict

import matplotlib
import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from figure6b_libr_ethanol_contributions import (
    AXIS_LABEL_SIZE,
    AXIS_TICK_SIZE,
    T_REF,
    _build_params,
    _calc_ln_miac_contributions,
    _load_exp_data,
    _molality_for_salt_mole_fraction,
    _salt_mole_fraction_from_molality,
)

matplotlib.use("Agg")
import matplotlib.pyplot as plt

VARIANTS: Dict[str, Dict[str, object]] = {
    "dion_0": {
        "label": r"2020 with $d_{ion}=0$",
        "color": "black",
        "linestyle": "--",
        "linewidth": 1.9,
        "user_options": {"elec_model": {"rel_perm": {"rule": 7}, "DH_model": {"d_ion_mode": 0}}},
    },
    "dion_1": {
        "label": r"2020 baseline, $d_{ion}=1$",
        "color": "green",
        "linestyle": "-",
        "linewidth": 2.1,
        "user_options": {"elec_model": {"rel_perm": {"rule": 7}}},
    },
    "dion_2": {
        "label": r"2020 with $d_{ion}=2$",
        "color": "tab:blue",
        "linestyle": "-.",
        "linewidth": 1.9,
        "user_options": {"elec_model": {"rel_perm": {"rule": 7}, "DH_model": {"d_ion_mode": 2}}},
    },
}


def run_analysis(
    data_path: Path,
    output_path: Path,
    x_min: float,
    x_max: float,
    y_min: float,
    y_max: float,
    grid_points: int,
    max_molality: float | None,
) -> Path:
    m_exp, x_exp, y_exp = _load_exp_data(data_path)

    m_upper = float(np.max(m_exp)) if max_molality is None else float(max_molality)
    if m_upper <= 0.0:
        m_upper = float(np.max(m_exp[m_exp >= 0.0]))

    m_grid = np.linspace(0.0, m_upper, int(grid_points))
    x_grid = _salt_mole_fraction_from_molality(m_grid)

    curves = {}
    for key, cfg in VARIANTS.items():
        params = _build_params(user_options=dict(cfg["user_options"]))
        curves[key] = _calc_ln_miac_contributions(m_grid, params)["dh"]

    fig, ax = plt.subplots(figsize=(8.2, 5.6))
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")

    ax.scatter(
        x_exp,
        y_exp,
        color="black",
        marker="o",
        s=34,
        facecolors="black",
        linewidths=0.8,
        label="Experimental data (Bulow 2020)",
        zorder=7,
    )

    for key, cfg in VARIANTS.items():
        ax.plot(
            x_grid,
            curves[key],
            color=str(cfg["color"]),
            linestyle=str(cfg["linestyle"]),
            linewidth=float(cfg["linewidth"]),
            label=str(cfg["label"]),
            zorder=5,
        )

    ax.set_xlim(float(x_min), float(x_max))
    ax.set_ylim(float(y_min), float(y_max))
    ax.set_xlabel(r"salt mole fraction, $x_{salt}$", fontsize=AXIS_LABEL_SIZE)
    ax.set_ylabel(r"DH contribution to $\ln(\gamma_{\pm}^{*})$", fontsize=AXIS_LABEL_SIZE)
    ax.set_title(r"LiBr in ethanol at 298.15 K and 1 bar (Figure 6b-style DH, rule 1a, $d_{ion}$ variants)")
    ax.grid(True, alpha=0.3, color="0.7")
    ax.tick_params(colors="black", labelsize=AXIS_TICK_SIZE)
    for spine in ax.spines.values():
        spine.set_color("black")
        spine.set_linewidth(1.0)

    legend = ax.legend(fontsize=8)
    frame = legend.get_frame()
    frame.set_facecolor("white")
    frame.set_edgecolor("black")
    frame.set_alpha(1.0)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(output_path, dpi=220)
    plt.close(fig)

    print(f"Loaded Bulow-2020 rows: {len(m_exp)} from {data_path}")
    print(f"Molality grid points: {len(m_grid)} (0 to {m_upper:.6g} mol/kg)")
    for key, cfg in VARIANTS.items():
        print(f"- {cfg['label']}")
    print(f"Wrote: {output_path}")
    return output_path


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Figure 6b-style LiBr/ethanol DH d_ion analysis")
    parser.add_argument(
        "--data",
        type=Path,
        default=REPO_ROOT / "data" / "MIAC" / "ethanol" / "ethanol-LiBr.csv",
        help="Input CSV with columns including molality, source, and miac/miac_m.",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=REPO_ROOT / "scripts" / "2020_model_analysis" / "output" / "figure6b_libr_ethanol_d_ion_variants.png",
        help="Output PNG path.",
    )
    parser.add_argument("--x-min", type=float, default=0.0)
    parser.add_argument("--x-max", type=float, default=0.2)
    parser.add_argument("--y-min", type=float, default=-3.0)
    parser.add_argument("--y-max", type=float, default=4.0)
    parser.add_argument("--grid-points", type=int, default=1201)
    parser.add_argument(
        "--max-molality",
        type=float,
        default=None,
        help="Optional upper molality bound for dense curve generation. Default: max Bulow-2020 molality.",
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = _parse_args()
    run_analysis(
        data_path=Path(args.data),
        output_path=Path(args.out),
        x_min=float(args.x_min),
        x_max=float(args.x_max),
        y_min=float(args.y_min),
        y_max=float(args.y_max),
        grid_points=int(args.grid_points),
        max_molality=None if args.max_molality is None else float(args.max_molality),
    )
