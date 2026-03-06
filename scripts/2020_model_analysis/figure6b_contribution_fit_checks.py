"""Generate one Figure 6b fit-comparison plot per contribution.

Each plot overlays the actual PC-SAFT 2020 model curve used in
figure6b_libr_ethanol_2020_contributions.png against the digitized contribution
points from Figure6b_curves.csv.
"""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path
from typing import Dict, Tuple

import matplotlib
import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPT_DIR = Path(__file__).resolve().parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from figure6b_libr_ethanol_contributions import (
    AXIS_LABEL_SIZE,
    AXIS_TICK_SIZE,
    _build_params,
    _calc_ln_miac_contributions,
    _load_exp_data,
    _salt_mole_fraction_from_molality,
)
from figure6b_digitized_reference_replica import (
    SERIES_STYLES,
    _load_digitized_curves,
)

matplotlib.use("Agg")
import matplotlib.pyplot as plt

CONTRIBUTION_ORDER = ["born", "dh", "hc", "disp", "assoc", "total"]
FILE_LABELS = {
    "born": "born",
    "dh": "dh",
    "hc": "hard_chain",
    "disp": "dispersion",
    "assoc": "association",
    "total": "total",
}
TITLE_LABELS = {
    "born": "Born",
    "dh": "Debye-Huckel",
    "hc": "Hard-chain",
    "disp": "Dispersion",
    "assoc": "Association",
    "total": "Total",
}


def _metric_summary(x_data: np.ndarray, y_data: np.ndarray, x_model: np.ndarray, y_model: np.ndarray) -> Tuple[float, float, float]:
    y_interp = np.interp(x_data, x_model, y_model)
    delta = y_interp - y_data
    rmse = float(np.sqrt(np.mean(delta * delta)))
    mae = float(np.mean(np.abs(delta)))
    max_abs = float(np.max(np.abs(delta)))
    return rmse, mae, max_abs


def _y_limits(y_data: np.ndarray, y_model: np.ndarray) -> Tuple[float, float]:
    y_all = np.concatenate([np.asarray(y_data, dtype=float), np.asarray(y_model, dtype=float)])
    y_min = float(np.min(y_all))
    y_max = float(np.max(y_all))
    span = max(y_max - y_min, 0.2)
    pad = 0.08 * span
    return y_min - pad, y_max + pad


def run_analysis(
    miac_data_path: Path,
    digitized_path: Path,
    output_dir: Path,
    grid_points: int,
    x_min: float,
    x_max: float,
) -> Dict[str, Dict[str, float]]:
    m_exp, _, _ = _load_exp_data(miac_data_path)
    m_upper = float(np.max(m_exp))
    m_grid = np.linspace(0.0, m_upper, int(grid_points))
    x_grid = _salt_mole_fraction_from_molality(m_grid)

    params = _build_params(user_options={})
    model_curves = _calc_ln_miac_contributions(m_grid, params)
    digitized = _load_digitized_curves(digitized_path)

    output_dir.mkdir(parents=True, exist_ok=True)
    results: Dict[str, Dict[str, float]] = {}

    for name in CONTRIBUTION_ORDER:
        x_data, y_data = digitized[name]
        y_model = np.asarray(model_curves[name], dtype=float)
        rmse, mae, max_abs = _metric_summary(x_data, y_data, x_grid, y_model)
        results[name] = {
            "rmse": rmse,
            "mae": mae,
            "max_abs": max_abs,
            "n_points": float(len(x_data)),
        }

        style = SERIES_STYLES[name]
        fig, ax = plt.subplots(figsize=(7.6, 5.2))
        fig.patch.set_facecolor("white")
        ax.set_facecolor("white")

        ax.scatter(
            x_data,
            y_data,
            color="black",
            marker="o",
            s=34,
            facecolors="black",
            linewidths=0.8,
            label="Digitized paper data",
            zorder=7,
        )
        ax.plot(
            x_grid,
            y_model,
            color=str(style["color"]),
            linestyle=str(style["linestyle"]),
            linewidth=float(style["linewidth"]),
            label="PC-SAFT 2020 model",
            zorder=5,
        )

        y_lo, y_hi = _y_limits(y_data, y_model)
        ax.set_xlim(float(x_min), float(x_max))
        ax.set_ylim(y_lo, y_hi)
        ax.set_xlabel(r"salt mole fraction, $x_{salt}$", fontsize=AXIS_LABEL_SIZE)
        ax.set_ylabel(r"Contribution to $\ln(\gamma_{\pm}^{*})$", fontsize=AXIS_LABEL_SIZE)
        ax.set_title(f"Figure 6b fit check: {TITLE_LABELS[name]}")
        ax.grid(True, alpha=0.3, color="0.7")
        ax.tick_params(colors="black", labelsize=AXIS_TICK_SIZE)
        for spine in ax.spines.values():
            spine.set_color("black")
            spine.set_linewidth(1.0)

        metric_text = f"RMSE = {rmse:.4f}\nMAE = {mae:.4f}\nMax |Δ| = {max_abs:.4f}"
        ax.text(
            0.98,
            0.97,
            metric_text,
            transform=ax.transAxes,
            ha="right",
            va="top",
            fontsize=9,
            bbox={"facecolor": "white", "edgecolor": "black", "alpha": 1.0, "boxstyle": "round,pad=0.25"},
        )

        legend = ax.legend(fontsize=9)
        frame = legend.get_frame()
        frame.set_facecolor("white")
        frame.set_edgecolor("black")
        frame.set_alpha(1.0)

        out_path = output_dir / f"figure6b_fit_{FILE_LABELS[name]}.png"
        fig.tight_layout()
        fig.savefig(out_path, dpi=220)
        plt.close(fig)
        print(f"{name}: rmse={rmse:.6f}, mae={mae:.6f}, max_abs={max_abs:.6f}, wrote={out_path}")

    return results


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate per-contribution Figure 6b model-vs-digitized fit plots")
    parser.add_argument(
        "--miac-data",
        type=Path,
        default=REPO_ROOT / "data" / "MIAC" / "ethanol" / "ethanol-LiBr.csv",
    )
    parser.add_argument(
        "--digitized",
        type=Path,
        default=Path(r"C:\Users\Tanner\Downloads\Figure6b_curves.csv"),
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=REPO_ROOT / "scripts" / "2020_model_analysis" / "output" / "figure6b_fit_checks",
    )
    parser.add_argument("--grid-points", type=int, default=1201)
    parser.add_argument("--x-min", type=float, default=0.0)
    parser.add_argument("--x-max", type=float, default=0.2)
    return parser.parse_args()


def main() -> None:
    args = _parse_args()
    run_analysis(
        miac_data_path=Path(args.miac_data),
        digitized_path=Path(args.digitized),
        output_dir=Path(args.out_dir),
        grid_points=int(args.grid_points),
        x_min=float(args.x_min),
        x_max=float(args.x_max),
    )


if __name__ == "__main__":
    main()
