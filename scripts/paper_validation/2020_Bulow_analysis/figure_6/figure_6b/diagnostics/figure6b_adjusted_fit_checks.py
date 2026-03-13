"""Generate Figure 6b fit-check plots from bookkeeping curves plus z-correction weights."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Dict, List, Tuple

import matplotlib
import numpy as np

from figure6b_digitized_reference_replica import SERIES_STYLES, _load_digitized_curves

matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO_ROOT = Path(__file__).resolve().parents[6]
OUTPUT_ROOT = REPO_ROOT / "scripts" / "paper_validation" / "2020_Bulow_analysis" / "figure_6" / "figure_6b" / "diagnostics" / "output"
OUTPUT_DATA_DIR = OUTPUT_ROOT / "data"
OUTPUT_PLOTS_DIR = OUTPUT_ROOT / "plots"

CONTRIBUTIONS = ["born", "dh", "hc", "disp", "assoc"]
FILE_LABELS = {
    "born": "born",
    "dh": "dh",
    "hc": "hard_chain",
    "disp": "dispersion",
    "assoc": "association",
}


def _read_dict_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def _interp(x_grid: np.ndarray, y_grid: np.ndarray, x: np.ndarray) -> np.ndarray:
    return np.interp(np.asarray(x, dtype=float), np.asarray(x_grid, dtype=float), np.asarray(y_grid, dtype=float))


def run_analysis(bookkeeping_csv: Path, digitized_csv: Path, weights_csv: Path, output_dir: Path) -> None:
    book_rows = _read_dict_rows(bookkeeping_csv)
    x_model = np.asarray([float(r["x_salt"]) for r in book_rows], dtype=float)
    zcorr = np.asarray([float(r["z_correction"]) for r in book_rows], dtype=float)
    model = {key: np.asarray([float(r[key]) for r in book_rows], dtype=float) for key in CONTRIBUTIONS}

    weight_rows = _read_dict_rows(weights_csv)
    weights = {str(r["contribution"]).strip(): float(r["best_scalar_z_weight"]) for r in weight_rows}

    digitized = _load_digitized_curves(digitized_csv)
    output_dir.mkdir(parents=True, exist_ok=True)

    for key in CONTRIBUTIONS:
        x_data, y_data = digitized[key]
        y_model = model[key] + weights[key] * zcorr
        y_interp = _interp(x_model, y_model, x_data)
        delta = y_interp - np.asarray(y_data, dtype=float)
        rmse = float(np.sqrt(np.mean(delta * delta)))
        mae = float(np.mean(np.abs(delta)))

        style = SERIES_STYLES[key]
        fig, ax = plt.subplots(figsize=(7.6, 5.2))
        fig.patch.set_facecolor("white")
        ax.set_facecolor("white")
        ax.scatter(x_data, y_data, color="black", s=34, linewidths=0.8, label="Digitized paper data", zorder=6)
        ax.plot(
            x_model,
            y_model,
            color=str(style["color"]),
            linestyle=str(style["linestyle"]),
            linewidth=float(style["linewidth"]),
            label=f"Adjusted model ({key})",
            zorder=4,
        )
        ax.set_xlim(0.0, 0.2)
        ax.set_xlabel(r"salt mole fraction, $x_{salt}$")
        ax.set_ylabel(r"Contribution to $\ln(\gamma_{\pm}^{*})$")
        ax.set_title(f"Adjusted Figure 6b fit check: {key}")
        ax.grid(True, alpha=0.3, color="0.7")
        ax.legend(fontsize=8)
        ax.text(
            0.98,
            0.97,
            f"w = {weights[key]:.4f}\nRMSE = {rmse:.4f}\nMAE = {mae:.4f}",
            transform=ax.transAxes,
            ha="right",
            va="top",
            fontsize=9,
            bbox={"facecolor": "white", "edgecolor": "black", "boxstyle": "round,pad=0.25"},
        )
        fig.tight_layout()
        fig.savefig(output_dir / f"figure6b_adjusted_fit_{FILE_LABELS[key]}.png", dpi=220)
        plt.close(fig)

    print(f"Wrote adjusted fit checks to {output_dir}")


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate adjusted Figure 6b fit-check plots from z-correction weights")
    parser.add_argument(
        "--bookkeeping-csv",
        type=Path,
        default=OUTPUT_DATA_DIR / "figure6b_bookkeeping.csv",
    )
    parser.add_argument(
        "--digitized-csv",
        type=Path,
        default=Path(r"C:\Users\Tanner\Documents\git\PC-SAFT\scripts\paper_validation\2020_Bulow_analysis\figure_6\figure_6b\data\Figure6b_curves.csv"),
    )
    parser.add_argument(
        "--weights-csv",
        type=Path,
        default=OUTPUT_DATA_DIR / "figure6b_accounting_best_scalar_zfits.csv",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=OUTPUT_PLOTS_DIR / "figure6b_fit_checks_best_scalar_zfits",
    )
    return parser.parse_args()


def main() -> None:
    args = _parse_args()
    run_analysis(
        bookkeeping_csv=Path(args.bookkeeping_csv),
        digitized_csv=Path(args.digitized_csv),
        weights_csv=Path(args.weights_csv),
        output_dir=Path(args.out_dir),
    )


if __name__ == "__main__":
    main()
