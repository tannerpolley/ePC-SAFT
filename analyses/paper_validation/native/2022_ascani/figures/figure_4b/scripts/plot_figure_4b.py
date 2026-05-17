from __future__ import annotations

import csv
import math
import sys
from pathlib import Path

for _candidate in Path(__file__).resolve().parents:
    if (_candidate / "scripts" / "plot_outputs.py").is_file():
        if str(_candidate) not in sys.path:
            sys.path.insert(0, str(_candidate))
        break
else:
    raise ModuleNotFoundError("Could not locate repo root containing scripts/plot_outputs.py")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from scripts.plot_outputs import figure_output_path, save_plot_figure

SQRT3_OVER_2 = math.sqrt(3.0) / 2.0


def _read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _style() -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Serif",
            "font.size": 10,
            "axes.linewidth": 1.0,
            "axes.grid": False,
            "legend.frameon": False,
            "xtick.direction": "in",
            "ytick.direction": "in",
        }
    )


def _draw_axes(ax: plt.Axes) -> None:
    triangle = np.asarray([[0.0, 0.0], [1.0, 0.0], [0.5, SQRT3_OVER_2], [0.0, 0.0]], dtype=float)
    ax.plot(triangle[:, 0], triangle[:, 1], color="black", linewidth=1.25, zorder=3)
    for frac in np.linspace(0.1, 0.9, 9):
        ax.plot([0.5 * frac, 1.0 - 0.5 * frac], [SQRT3_OVER_2 * frac] * 2, color="0.86", linewidth=0.55, zorder=0)
        ax.plot([frac, 0.5 + 0.5 * frac], [0.0, SQRT3_OVER_2 * (1.0 - frac)], color="0.86", linewidth=0.55, zorder=0)
        ax.plot([1.0 - frac, 0.5 * (1.0 - frac)], [0.0, SQRT3_OVER_2 * (1.0 - frac)], color="0.86", linewidth=0.55, zorder=0)
    ax.text(-0.04, -0.04, r"$H_2O$", ha="right", va="top", fontsize=12)
    ax.text(1.04, -0.04, "1-butanol", ha="left", va="top", fontsize=12)
    ax.text(0.5, SQRT3_OVER_2 + 0.04, "total salt", ha="center", va="bottom", fontsize=12)
    ax.set_aspect("equal")
    ax.set_xlim(-0.08, 1.08)
    ax.set_ylim(-0.08, SQRT3_OVER_2 + 0.08)
    ax.axis("off")


def _plot_series(ax: plt.Axes, rows: list[dict[str, str]], series: str, *, color: str, marker: str, label: str, linestyle: str) -> None:
    phase_order = ["aq", "org"] if series != "current_feed" else ["feed"]
    pts = []
    for phase in phase_order:
        match = next(row for row in rows if row["series"] == series and row["phase"] == phase)
        pts.append((float(match["x_plot"]), float(match["y_plot"])))
    arr = np.asarray(pts, dtype=float)
    if arr.shape[0] > 1:
        ax.plot(arr[:, 0], arr[:, 1], color=color, linestyle=linestyle, linewidth=1.6, zorder=5, label=label)
    else:
        ax.scatter(arr[:, 0], arr[:, 1], color=color, marker=marker, s=55, zorder=7, label=label)
        return
    ax.scatter(arr[:, 0], arr[:, 1], facecolor="white", edgecolor=color, marker=marker, s=62, linewidth=1.2, zorder=6)


def main() -> None:
    _style()
    data_path = figure_output_path(__file__, "data/figure_4b_phase_diagram.csv")
    rows = _read_rows(data_path)
    fig, ax = plt.subplots(figsize=(6.2, 5.6))
    _draw_axes(ax)
    _plot_series(ax, rows, "paper_case2", color="0.45", marker="o", label="Ascani 2022 Case 2", linestyle=":")
    _plot_series(ax, rows, "current_ipopt", color="black", marker="s", label="current native Ipopt", linestyle="-")
    _plot_series(ax, rows, "current_feed", color="#2ca25f", marker="^", label="current feed", linestyle="")
    ax.legend(loc="upper right", fontsize=8)
    ax.set_title("Ascani 2022 Figure 4b-style phase split", fontsize=13)
    output = figure_output_path(__file__, "figure_4b.png")
    save_plot_figure(fig, output, dpi=300, svg_companion=True)
    plt.close(fig)
    print(f"[write] {output}")


if __name__ == "__main__":
    main()
