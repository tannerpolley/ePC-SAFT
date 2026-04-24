"""Plot a digitized replica of Figure 6b from curve CSV exports."""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path
from typing import Dict, List, Tuple

import matplotlib
import numpy as np
from scipy.interpolate import PchipInterpolator

matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO_ROOT = Path(__file__).resolve().parents[6]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts.plot_outputs import paper_validation_path

OUTPUT_ROOT = paper_validation_path(Path(__file__).resolve().parent, "output")
OUTPUT_PLOTS_DIR = OUTPUT_ROOT / "plots"

AXIS_LABEL_SIZE = 12
AXIS_TICK_SIZE = 11

SERIES_STYLES = {
    "data": {"label": "Experimental data", "kind": "scatter", "color": "black"},
    "born": {"label": "Born contribution", "color": "orange", "linestyle": "-", "linewidth": 1.9},
    "dh": {"label": "DH contribution", "color": "tab:blue", "linestyle": "-", "linewidth": 1.9},
    "hc": {"label": "Hard-chain contribution", "color": "gray", "linestyle": "-", "linewidth": 1.8},
    "disp": {"label": "Dispersion contribution", "color": "gray", "linestyle": "--", "linewidth": 1.8},
    "assoc": {"label": "Association contribution", "color": "gray", "linestyle": "-.", "linewidth": 1.8},
    "total": {"label": "Total (2020)", "color": "green", "linestyle": "-", "linewidth": 2.1},
}


def _normalize_series_name(name: str) -> str:
    s = str(name or "").strip().lower()
    aliases = {
        "born": "born",
        "hc": "hc",
        "hard-chain": "hc",
        "hard_chain": "hc",
        "disp": "disp",
        "dispersion": "disp",
        "assoc": "assoc",
        "association": "assoc",
        "dh": "dh",
        "total": "total",
        "data": "data",
    }
    if s not in aliases:
        raise ValueError(f"Unknown Figure 6b series name: {name}")
    return aliases[s]


def _load_digitized_curves(path: Path) -> Dict[str, Tuple[np.ndarray, np.ndarray]]:
    with path.open("r", newline="", encoding="utf-8-sig") as handle:
        reader = csv.reader(handle)
        rows = list(reader)

    if not rows:
        raise ValueError(f"No content found in {path}")

    header = [str(cell).strip() for cell in rows[0]]
    if len(header) % 2 != 0:
        raise ValueError("Expected x/y column pairs in the digitized Figure 6b CSV.")

    series: Dict[str, List[Tuple[float, float]]] = {}
    for col in range(0, len(header), 2):
        y_name = _normalize_series_name(header[col + 1])
        pairs: List[Tuple[float, float]] = []
        for row in rows[1:]:
            x_raw = row[col].strip() if col < len(row) and row[col] is not None else ""
            y_raw = row[col + 1].strip() if col + 1 < len(row) and row[col + 1] is not None else ""
            if not x_raw or not y_raw:
                continue
            x_val = float(x_raw)
            y_val = float(y_raw)
            if np.isfinite(x_val) and np.isfinite(y_val):
                pairs.append((x_val, y_val))
        if not pairs:
            continue
        pairs.sort(key=lambda item: item[0])
        x_arr = np.asarray([p[0] for p in pairs], dtype=float)
        y_arr = np.asarray([p[1] for p in pairs], dtype=float)
        series[y_name] = (x_arr, y_arr)

    missing = [name for name in SERIES_STYLES if name not in series]
    if missing:
        raise ValueError(f"Digitized Figure 6b CSV is missing series: {missing}")
    return series


def _interp_dense_curve(x: np.ndarray, y: np.ndarray, points: int) -> Tuple[np.ndarray, np.ndarray]:
    x_dense = np.linspace(float(np.min(x)), float(np.max(x)), int(points))
    if len(x) < 3 or np.any(np.diff(x) <= 0.0):
        y_dense = np.interp(x_dense, x, y)
        return x_dense, y_dense
    interp = PchipInterpolator(x, y)
    y_dense = interp(x_dense)
    return x_dense, y_dense


def run_analysis(
    data_path: Path,
    output_path: Path,
    x_min: float,
    x_max: float,
    y_min: float,
    y_max: float,
    interp_points: int,
) -> Path:
    series = _load_digitized_curves(data_path)

    fig, ax = plt.subplots(figsize=(8.2, 5.6))
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")

    x_data, y_data = series["data"]
    ax.scatter(
        x_data,
        y_data,
        color="black",
        marker="o",
        s=34,
        facecolors="black",
        linewidths=0.8,
        label="Experimental data",
        zorder=7,
    )

    for name in ["born", "dh", "hc", "disp", "assoc", "total"]:
        x, y = series[name]
        x_dense, y_dense = _interp_dense_curve(x, y, interp_points)
        style = SERIES_STYLES[name]
        ax.plot(
            x_dense,
            y_dense,
            color=str(style["color"]),
            linestyle=str(style["linestyle"]),
            linewidth=float(style["linewidth"]),
            label=str(style["label"]),
            zorder=5,
        )

    ax.set_xlim(float(x_min), float(x_max))
    ax.set_ylim(float(y_min), float(y_max))
    ax.set_xlabel(r"salt mole fraction, $x_{salt}$", fontsize=AXIS_LABEL_SIZE)
    ax.set_ylabel(r"$\ln(\gamma_{\pm}^{*})$", fontsize=AXIS_LABEL_SIZE)
    ax.set_title("Digitized Figure 6b replica")
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

    print(f"Loaded digitized series from {data_path}")
    for name, (x, _) in series.items():
        print(f"- {name}: {len(x)} points")
    print(f"Wrote: {output_path}")
    return output_path


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Plot a digitized Figure 6b replica from CSV curves")
    parser.add_argument(
        "--data",
        type=Path,
        default=Path(r"C:\Users\Tanner\Documents\git\ePC-SAFT\scripts\paper_validation\2020_Bulow_analysis\figure_6\figure_6b\data\Figure6b_curves.csv"),
        help="Digitized Figure 6b CSV with repeated x/y column pairs.",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=OUTPUT_PLOTS_DIR / "figure6b_digitized_reference_replica.png",
        help="Output PNG path.",
    )
    parser.add_argument("--x-min", type=float, default=0.0)
    parser.add_argument("--x-max", type=float, default=0.2)
    parser.add_argument("--y-min", type=float, default=-3.0)
    parser.add_argument("--y-max", type=float, default=4.0)
    parser.add_argument("--interp-points", type=int, default=4000)
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
        interp_points=int(args.interp_points),
    )

