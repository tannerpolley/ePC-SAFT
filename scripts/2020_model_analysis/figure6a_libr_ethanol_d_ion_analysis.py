"""Figure 6a-style MIAC analysis for LiBr in ethanol with d_ion variants."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Dict, List, Tuple

import matplotlib
import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from figure6a_libr_ethanol_analysis import (
    _calc_curve,
    _color_for_source,
    _high_outlier_mask,
    _load_exp_data,
    _molality_for_salt_mole_fraction,
    _salt_mole_fraction_from_molality,
)

matplotlib.use("Agg")
import matplotlib.pyplot as plt


VARIANTS: List[Dict[str, object]] = [
    {
        "label": r"2020 with $d_{ion}=0$",
        "dataset": "bulow_2020",
        "color": "black",
        "lw": 2.0,
        "linestyle": "--",
        "zorder": 3,
        "user_options": {"elec_model": {"DH_model": {"d_ion_mode": 0}}},
    },
    {
        "label": r"2020 baseline, $d_{ion}=1$",
        "dataset": "bulow_2020",
        "color": "green",
        "lw": 2.2,
        "linestyle": "-",
        "zorder": 5,
        "user_options": {},
    },
    {
        "label": r"2020 with $d_{ion}=2$",
        "dataset": "bulow_2020",
        "color": "tab:blue",
        "lw": 2.0,
        "linestyle": "-.",
        "zorder": 4,
        "user_options": {"elec_model": {"DH_model": {"d_ion_mode": 2}}},
    },
]


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
    m_exp_raw, x_exp_raw, y_exp_raw, source_exp_raw = _load_exp_data(data_path)

    keep = _high_outlier_mask(y_exp_raw)
    removed = int(np.count_nonzero(~keep))
    if np.any(keep):
        m_exp = m_exp_raw[keep]
        x_exp = x_exp_raw[keep]
        y_exp = y_exp_raw[keep]
        source_exp = source_exp_raw[keep]
    else:
        m_exp = m_exp_raw
        x_exp = x_exp_raw
        y_exp = y_exp_raw
        source_exp = source_exp_raw
        removed = 0

    if removed > 0:
        print(f"[outlier-filter] d_ion figure removed {removed} high outlier experimental point(s).")

    m_upper = float(np.max(m_exp)) if max_molality is None else float(max_molality)
    if m_upper <= 0.0:
        m_upper = float(np.max(m_exp[m_exp >= 0.0]))
    m_upper = max(m_upper, _molality_for_salt_mole_fraction(x_max))

    m_grid = np.linspace(0.0, m_upper, int(grid_points))
    x_grid = _salt_mole_fraction_from_molality(m_grid)

    curves: List[Tuple[Dict[str, object], np.ndarray]] = []
    for cfg in VARIANTS:
        y_curve = _calc_curve(m_grid, str(cfg["dataset"]), user_options=dict(cfg["user_options"]))
        curves.append((cfg, y_curve))

    fig, ax = plt.subplots(figsize=(8.2, 5.6))
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")

    unique_sources = []
    for src in source_exp.tolist():
        if src not in unique_sources:
            unique_sources.append(src)

    if len(unique_sources) == 1 and unique_sources[0] == "Unspecified source":
        ax.scatter(
            x_exp,
            y_exp,
            color="black",
            marker="o",
            s=38,
            facecolors="none",
            linewidths=1.1,
            label="Experimental data",
            zorder=7,
        )
    else:
        for src in unique_sources:
            mask = source_exp == src
            color = _color_for_source(str(src))
            lbl = f"Experimental data - {src}" if src != "Unspecified source" else "Experimental data - Unspecified source"
            ax.scatter(
                x_exp[mask],
                y_exp[mask],
                color=color,
                marker="o",
                s=40,
                facecolors="none",
                linewidths=1.15,
                label=lbl,
                zorder=7,
            )

    for cfg, y_curve in curves:
        ax.plot(
            x_grid,
            y_curve,
            color=str(cfg["color"]),
            linewidth=float(cfg["lw"]),
            linestyle=str(cfg.get("linestyle", "-")),
            zorder=float(cfg.get("zorder", 3)),
            label=str(cfg["label"]),
        )

    ax.set_xlim(float(x_min), float(x_max))
    ax.set_ylim(float(y_min), float(y_max))
    ax.set_xlabel(r"salt mole fraction, $x_{salt}$")
    ax.set_ylabel(r"mean ionic activity coefficient, $\gamma_{\pm}^{*}$")
    ax.set_title(r"LiBr in ethanol at 298.15 K and 1 bar (Figure 6a-style, $d_{ion}$ variants)")
    ax.grid(True, alpha=0.3, color="0.7")
    ax.tick_params(colors="black")
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

    if not output_path.exists():
        raise FileNotFoundError(f"Expected plot was not written: {output_path}")

    print(f"Loaded rows: {len(m_exp)} from {data_path}")
    src_counts = {s: int(np.sum(source_exp == s)) for s in np.unique(source_exp)}
    print(f"Data sources: {src_counts}")
    print(f"Molality grid points: {len(m_grid)} (0 to {m_upper:.6g} mol/kg)")
    print("Curves:")
    for cfg, _ in curves:
        print(f"- {cfg['label']} [{cfg['dataset']}]")
    print(f"Wrote: {output_path}")
    return output_path


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Figure 6a-style LiBr/ethanol MIAC d_ion analysis")
    parser.add_argument(
        "--data",
        type=Path,
        default=REPO_ROOT / "data" / "MIAC" / "ethanol" / "ethanol-LiBr.csv",
        help="Input CSV with columns including molality, mole_fraction, miac.",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=REPO_ROOT / "scripts" / "2020_model_analysis" / "output" / "figure6a_libr_ethanol_d_ion_variants.png",
        help="Output PNG path.",
    )
    parser.add_argument("--x-min", type=float, default=0.0)
    parser.add_argument("--x-max", type=float, default=0.2)
    parser.add_argument("--y-min", type=float, default=0.0)
    parser.add_argument("--y-max", type=float, default=4.0)
    parser.add_argument("--grid-points", type=int, default=1201)
    parser.add_argument(
        "--max-molality",
        type=float,
        default=None,
        help="Optional upper molality bound for dense curve generation. Default: max experimental molality.",
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
