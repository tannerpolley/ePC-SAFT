from __future__ import annotations

from pathlib import Path
import sys

import numpy as np
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import _common as common

OUTPUT = Path(__file__).with_name("figure_6.png")
PANELS = [
    ("a)", "K+", "methanol", common.REPO_ROOT / "data" / "G_trans" / "water" / "methanol" / "K.csv", r"$x_{MeOH}$ / -", (0.0, 12.5)),
    ("b)", "Br-", "methanol", common.REPO_ROOT / "data" / "G_trans" / "water" / "methanol" / "Br.csv", r"$x_{MeOH}$ / -", (0.0, 12.0)),
    ("c)", "Na+", "ethanol", common.REPO_ROOT / "data" / "G_trans" / "water" / "ethanol" / "Na.csv", r"$x_{EtOH}$ / -", (0.0, 20.0)),
    ("d)", "Cl-", "ethanol", common.REPO_ROOT / "data" / "G_trans" / "water" / "ethanol" / "Cl.csv", r"$x_{EtOH}$ / -", (0.0, 20.0)),
]


def _load_xy(path: Path):
    fields, rows = common.read_csv_rows(path)
    x_key = fields[0]
    y_key = fields[1]
    xs, ys = [], []
    for row in rows:
        x = common.parse_float(row.get(x_key))
        y = common.parse_float(row.get(y_key))
        if x is None or y is None:
            continue
        xs.append(x)
        ys.append(y)
    return np.asarray(xs, dtype=float), np.asarray(ys, dtype=float)


def main() -> None:
    common.configure_style()
    fig, axes = plt.subplots(2, 2, figsize=(6.6, 5.8))
    for ax, (label, ion, organic, csv_path, xlabel, ylim) in zip(axes.flat, PANELS):
        x_data, y_data = _load_xy(csv_path)
        x_grid = np.linspace(0.0, 1.0, 401)
        y_model = common.transfer_curve("figiel_2025", ion, organic, x_grid)
        ax.plot(x_grid, y_model, color="black", linewidth=1.5)
        ax.scatter(x_data, y_data, s=24, facecolor=common.LIGHT_GRAY, edgecolor=common.GRAY_COLOR, linewidth=0.8)
        common.panel_label(ax, label)
        ax.set_xlim(0.0, 1.0)
        ax.set_ylim(*ylim)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(r"$\Delta G_i^{trans,\infty,x}$ / kJ mol$^{-1}$")
    common.save_figure(fig, OUTPUT)
    plt.close(fig)


if __name__ == "__main__":
    main()
