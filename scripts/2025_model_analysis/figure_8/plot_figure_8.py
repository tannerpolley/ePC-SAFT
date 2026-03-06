from __future__ import annotations

from pathlib import Path
import sys

import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import _common as common

OUTPUT = Path(__file__).with_name("figure_8.png")
PANELS = [
    ("a)", "LiBr", 5.0, 3.0),
    ("b)", "NaI", 1.5, 1.125),
    ("c)", "NaBr", 1.5, 1.125),
]


def _plot_panel(ax, label, salt, m_max, y_max):
    methanol_data = common.read_miac_dataset(common.REPO_ROOT / "data" / "MIAC" / "methanol" / f"methanol-{salt}.csv", "methanol")
    ethanol_data = common.read_miac_dataset(common.REPO_ROOT / "data" / "MIAC" / "ethanol" / f"ethanol-{salt}.csv", "ethanol")
    m_grid_meoh, y_meoh = common.mean_ionic_activity_curve("figiel_2025", salt, "methanol", {"methanol": 1.0}, m_max, points=600)
    m_grid_etoh, y_etoh = common.mean_ionic_activity_curve("figiel_2025", salt, "ethanol", {"ethanol": 1.0}, m_max, points=600)
    ax.scatter([r["molality"] for r in methanol_data], [r["miac_m"] for r in methanol_data], s=24, facecolor="none", edgecolor=common.GRAY_COLOR, linewidth=0.9)
    ax.scatter([r["molality"] for r in ethanol_data], [r["miac_m"] for r in ethanol_data], s=24, marker="s", facecolor=common.GREEN_COLOR, edgecolor=common.GREEN_COLOR, linewidth=0.8)
    ax.plot(m_grid_meoh, y_meoh, color=common.GRAY_COLOR, linewidth=1.5)
    ax.plot(m_grid_etoh, y_etoh, color="black", linewidth=1.5)
    common.panel_label(ax, label)
    ax.set_xlim(0.0, m_max)
    ax.set_ylim(0.0, y_max)
    ax.set_title(salt, fontsize=10)
    ax.set_xlabel(r"$\bar{m}_{salt}$ / mol kg$^{-1}$")
    ax.set_ylabel(r"$\gamma_{\pm}^{m,*}$ / -")


def main() -> None:
    common.configure_style()
    fig = plt.figure(figsize=(6.6, 6.3))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.0])
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[1, :])
    for ax, cfg in zip([ax1, ax2, ax3], PANELS):
        _plot_panel(ax, *cfg)
    handles = [
        plt.Line2D([0], [0], marker="o", linestyle="None", markerfacecolor="none", markeredgecolor=common.GRAY_COLOR, color=common.GRAY_COLOR, label="Methanol data"),
        plt.Line2D([0], [0], marker="s", linestyle="None", markerfacecolor=common.GREEN_COLOR, markeredgecolor=common.GREEN_COLOR, color=common.GREEN_COLOR, label="Ethanol data"),
        plt.Line2D([0], [0], color=common.GRAY_COLOR, linewidth=1.5, label="Methanol model"),
        plt.Line2D([0], [0], color="black", linewidth=1.5, label="Ethanol model"),
    ]
    fig.legend(handles=handles, loc="upper center", ncol=2, bbox_to_anchor=(0.5, 1.02), fontsize=9)
    common.save_figure(fig, OUTPUT)
    plt.close(fig)


if __name__ == "__main__":
    main()
