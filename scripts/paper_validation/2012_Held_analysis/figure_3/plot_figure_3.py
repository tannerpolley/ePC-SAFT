from __future__ import annotations

from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import _common as common

DATA_PATH = Path(__file__).with_name("data.csv")
OUTPUT = Path(__file__).with_name("figure_3.png")

SERIES = [
    {
        "csv_name": "H2O",
        "solvent_system": "water",
        "label": "water",
        "marker": "*",
        "color": "#1f5aa6",
        "size": 58,
        "facecolors": "#1f5aa6",
    },
    {
        "csv_name": "methanol",
        "solvent_system": "methanol",
        "label": "methanol",
        "marker": "s",
        "color": "#17632d",
        "size": 28,
        "facecolors": "none",
    },
    {
        "csv_name": "ethanol",
        "solvent_system": "ethanol",
        "label": "ethanol",
        "marker": "o",
        "color": "#b55d09",
        "size": 28,
        "facecolors": "none",
    },
]


def load_series_rows() -> dict[str, list[dict[str, float | str]]]:
    _, rows = common.read_csv_rows(DATA_PATH)
    grouped: dict[str, list[dict[str, float | str]]] = {series["csv_name"]: [] for series in SERIES}
    for row in rows:
        solvent = str(row.get("solvent", "")).strip()
        if solvent not in grouped:
            continue
        m = common.parse_float(row.get("m"))
        gamma = common.parse_float(row.get("miac"))
        if m is None or gamma is None:
            continue
        grouped[solvent].append({"m": m, "gamma": gamma, "solvent": solvent})
    for values in grouped.values():
        values.sort(key=lambda item: float(item["m"]))
    return grouped


def main() -> None:
    common.configure_style()
    grouped = load_series_rows()

    fig, ax = plt.subplots(figsize=(6.9, 4.4))
    ymin = np.inf
    ymax = -np.inf
    legend_handles = []

    for series in SERIES:
        rows = grouped[series["csv_name"]]
        if not rows:
            raise ValueError(f"No usable rows found for {series['csv_name']} in {DATA_PATH}.")

        m_data = np.asarray([float(row["m"]) for row in rows], dtype=float)
        y_data = np.asarray([float(row["gamma"]) for row in rows], dtype=float)

        grid, y_model = common.solvent_activity_curve(
            dataset="held_2012",
            salt="NaI",
            solvent_system=series["solvent_system"],
            m_max=1.6,
            points=500,
        )

        ax.plot(grid, y_model, color=series["color"], linewidth=2.0)
        ax.scatter(
            m_data,
            y_data,
            s=series["size"],
            marker=series["marker"],
            facecolors=series["facecolors"],
            edgecolors=series["color"],
            linewidths=1.0,
            zorder=5,
        )

        ymin = min(ymin, float(np.min(y_data)), float(np.min(y_model)))
        ymax = max(ymax, float(np.max(y_data)), float(np.max(y_model)))
        legend_handles.append(
            Line2D(
                [0],
                [0],
                color=series["color"],
                marker=series["marker"],
                markerfacecolor=series["facecolors"],
                markeredgecolor=series["color"],
                linewidth=2.0,
                markersize=8,
                label=series["label"],
            )
        )

    ax.set_xlim(0.0, 1.6)
    ax.set_ylim(max(0.97, ymin - 0.002), min(1.02, ymax + 0.002))
    ax.set_xlabel(r"molality, $m$ / mol kg$^{-1}$")
    ax.set_ylabel(r"$\gamma_{solv}$ / -")
    ax.set_title("2012 Fig. 3 style: NaI solvent activity coefficients (298.15 K)")
    ax.grid(True, alpha=0.25)
    ax.legend(handles=legend_handles, loc="best", fontsize=9)

    common.save_figure(fig, OUTPUT)
    plt.close(fig)


if __name__ == "__main__":
    main()
