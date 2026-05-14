from __future__ import annotations

from functools import lru_cache
from pathlib import Path
import sys

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ANALYSIS_ROOT = Path(__file__).resolve().parents[3]
if str(ANALYSIS_ROOT) not in sys.path:
    sys.path.insert(0, str(ANALYSIS_ROOT))

from shared import _common as common
from shared import figure_data

OUTPUT = Path(__file__).with_name("figure_6.png")
OUTPUT_METHANOL = Path(__file__).with_name("figure_6_methanol.png")
OUTPUT_ETHANOL = Path(__file__).with_name("figure_6_ethanol.png")
OUTPUTS = [
    Path(__file__).with_name("figure_6a.png"),
    Path(__file__).with_name("figure_6b.png"),
    Path(__file__).with_name("figure_6c.png"),
    Path(__file__).with_name("figure_6d.png"),
]
PANELS = [
    ("a)", "K+", "methanol", r"$x_{MeOH}$ / -", (0.0, 12.5)),
    ("b)", "Br-", "methanol", r"$x_{MeOH}$ / -", (0.0, 12.0)),
    ("c)", "Na+", "ethanol", r"$x_{EtOH}$ / -", (0.0, 20.0)),
    ("d)", "Cl-", "ethanol", r"$x_{EtOH}$ / -", (0.0, 30.0)),
]


TITLE_ALL = "Gibbs energies of transfer at infinite dilution from water to water + organic solvent systems\n$\\Delta G_i^{trans,\\infty}$ of different ions at 298.15 K and 1 bar."
TITLE_METHANOL = "Gibbs energies of transfer at infinite dilution from water to water + $\\mathbf{methanol}$ systems\n$\\Delta G_i^{trans,\\infty}$ of different ions at 298.15 K and 1 bar."
TITLE_ETHANOL = "Gibbs energies of transfer at infinite dilution from water to water + $\\mathbf{ethanol}$ systems\n$\\Delta G_i^{trans,\\infty}$ of different ions at 298.15 K and 1 bar."


@lru_cache(maxsize=1)
def _payload():
    return tuple(figure_data.read_payload("figure_6"))


def _plot_panel(ax, label, ion, organic, xlabel, ylim) -> None:
    rows = list(_payload())
    x_data, y_data = figure_data.xy(figure_data.select_rows(rows, panel_id=label, series_id=f"data_{ion}_{organic}"))
    x_grid, y_model = figure_data.xy(figure_data.select_rows(rows, panel_id=label, series_id=f"model_{ion}_{organic}"))
    ax.plot(x_grid, y_model, color="black", linewidth=1.5)
    ax.scatter(x_data, y_data, s=24, facecolor=common.LIGHT_GRAY, edgecolor=common.GRAY_COLOR, linewidth=0.8)
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(*ylim)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(r"$\Delta G_i^{trans,\infty,x}$ / kJ mol$^{-1}$")


def _legend_handles():
    return [
        Line2D([0], [0], color="black", linewidth=1.5, label="ePC-SAFT fit"),
        Line2D(
            [0],
            [0],
            marker="o",
            linestyle="None",
            markersize=5.5,
            markerfacecolor=common.LIGHT_GRAY,
            markeredgecolor=common.GRAY_COLOR,
            label="Literature data",
        ),
    ]


def _save_combined(
    panels, output: Path, title: str, shape: tuple[int, int], figsize: tuple[float, float], top: float, legend_y: float
):
    common.configure_style()
    fig, axes = plt.subplots(*shape, figsize=figsize)
    axes_flat = np.atleast_1d(axes).ravel()
    for ax, cfg in zip(axes_flat, panels):
        _plot_panel(ax, *cfg)
    handles = _legend_handles()
    fig.legend(handles=handles, loc="upper center", ncol=2, bbox_to_anchor=(0.5, legend_y), fontsize=9, frameon=False)
    fig.suptitle(title, fontsize=11, y=0.995)
    fig.subplots_adjust(left=0.10, right=0.98, bottom=0.14, top=top, wspace=0.24, hspace=0.28)
    common.save_figure(fig, output)
    plt.close(fig)


def main() -> None:
    _save_combined(PANELS, OUTPUT, TITLE_ALL, (2, 2), (7.6, 6.3), 0.87, 0.945)
    _save_combined(PANELS[:2], OUTPUT_METHANOL, TITLE_METHANOL, (1, 2), (7.6, 3.8), 0.72, 0.88)
    _save_combined(PANELS[2:], OUTPUT_ETHANOL, TITLE_ETHANOL, (1, 2), (7.6, 3.8), 0.72, 0.88)

    for cfg, out in zip(PANELS, OUTPUTS):
        common.save_panel_figure(lambda ax, cfg=cfg: _plot_panel(ax, *cfg), out, figsize=(3.8, 3.3))


if __name__ == "__main__":
    main()
