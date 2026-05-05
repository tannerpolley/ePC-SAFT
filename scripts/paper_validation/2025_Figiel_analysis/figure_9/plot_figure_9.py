from __future__ import annotations

from functools import lru_cache
from pathlib import Path
import sys

import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import _common as common
import figure_data

OUTPUT = Path(__file__).with_name("figure_9.png")
OUTPUT_40_METHANOL = Path(__file__).with_name("figure_9_methanol_40wt.png")
OUTPUT_40_ETHANOL = Path(__file__).with_name("figure_9_ethanol_40wt.png")
OUTPUT_40_COMBINED = Path(__file__).with_name("figure_9_40wt_combined.png")
OUTPUTS = [
    Path(__file__).with_name("figure_9a.png"),
    Path(__file__).with_name("figure_9b.png"),
    Path(__file__).with_name("figure_9c.png"),
    Path(__file__).with_name("figure_9d.png"),
]
PANELS = figure_data.FIG9_PANELS
DEFAULT_TARGETS = [
    (
        0.8,
        common.ORGANIC_COLOR,
        "^",
        common.ORGANIC_COLOR,
        "-",
        r"Model, $w_{org}^{salt-free}=0.8$",
        r"$w_{org}^{salt-free}=0.8$",
    ),
    (
        0.4,
        common.GREEN_COLOR,
        "s",
        common.GREEN_COLOR,
        "-",
        r"Model, $w_{org}^{salt-free}=0.4$",
        r"$w_{org}^{salt-free}=0.4$",
    ),
]
SALT_STYLES = {
    "NaBr": {"color": common.ORGANIC_COLOR, "marker": "^"},
    "NaCl": {"color": common.BLUE_COLOR, "marker": "s"},
}


@lru_cache(maxsize=1)
def _payload():
    return tuple(figure_data.read_payload("figure_9"))


def _series_xy(series_id: str, panel_id: str):
    return figure_data.xy(figure_data.select_rows(list(_payload()), panel_id=panel_id, series_id=series_id))


def _plot_panel(ax, label: str, salt: str, solvent_system: str, m_max: float, include_legend: bool = False) -> None:
    for target_w, marker_color, marker, line_color, line_style, line_label, data_label_text in DEFAULT_TARGETS:
        target_token = f"{target_w:.1f}"
        m_grid, y_model = _series_xy(f"model_{salt}_{solvent_system}_w{target_token}", label)
        ax.plot(
            m_grid,
            y_model,
            color=line_color,
            linestyle=line_style,
            linewidth=2.2,
            zorder=5,
            label=line_label if include_legend else None,
        )

        data_m, data_y = _series_xy(f"data_{salt}_{solvent_system}_w{target_token}", label)
        if data_m:
            ax.scatter(
                data_m,
                data_y,
                s=24,
                marker=marker,
                facecolor="none",
                edgecolor=marker_color,
                linewidth=1.0,
                zorder=6,
                label=data_label_text if include_legend else None,
            )

    organic = [s for s in solvent_system.split("-") if s != "water"][0]
    organic_label = "methanol" if organic == "methanol" else "ethanol"
    ax.set_xlim(0.0, m_max)
    ax.set_ylim(0.0, 1.125)
    ax.set_title(f"{salt} in aqueous {organic_label}", fontsize=10, pad=8)
    ax.set_xlabel(r"$\bar{m}_{salt}$ / mol kg$^{-1}$", labelpad=4)
    ax.set_ylabel(r"$\gamma_{\pm}^{m,*}$ of salt / -", labelpad=4)
    if include_legend:
        ax.legend(loc="upper right", fontsize=7.5, frameon=False)


def _plot_40wt_solvent(
    ax, solvent_system: str, organic_label: str, salts: list[str], x_max: float, y_max: float
) -> None:
    panel_id = f"40wt_{solvent_system}"
    for salt in salts:
        style = SALT_STYLES[salt]
        data_m, data_y = _series_xy(f"data_{salt}_{solvent_system}_40wt", panel_id)
        m_grid, y_model = _series_xy(f"model_{salt}_{solvent_system}_40wt", panel_id)
        ax.plot(m_grid, y_model, color=style["color"], linewidth=2.2, label=f"{salt} fit")
        ax.scatter(
            data_m,
            data_y,
            s=28,
            marker=style["marker"],
            facecolor="none",
            edgecolor=style["color"],
            linewidth=1.0,
            label=f"{salt} data",
        )

    ax.set_xlim(0.0, x_max)
    ax.set_ylim(0.0, y_max)
    panel_title = (
        "40 wt% organic in aqueous methanol" if organic_label == "methanol" else "40 wt% organic in aqueous ethanol"
    )
    ax.set_title(panel_title, fontsize=10, pad=8, fontweight="bold")
    ax.set_xlabel(r"$\bar{m}_{salt}$ / mol kg$^{-1}$")
    ax.set_ylabel(r"$\gamma_{\pm}^{m,*}$ of salt / -")
    ax.legend(loc="upper right", fontsize=8, frameon=False)


def main() -> None:
    common.configure_style()
    fig, axes = plt.subplots(2, 2, figsize=(9.3, 8.1))
    for ax, cfg in zip(axes.flat, PANELS):
        _plot_panel(ax, *cfg, include_legend=False)
    handles = [
        plt.Line2D([0], [0], color=common.ORGANIC_COLOR, linewidth=2.2, label=r"Model, $w_{org}^{salt-free}=0.8$"),
        plt.Line2D(
            [0],
            [0],
            marker="^",
            linestyle="None",
            markerfacecolor="none",
            markeredgecolor=common.ORGANIC_COLOR,
            color=common.ORGANIC_COLOR,
            label=r"$w_{org}^{salt-free}=0.8$",
        ),
        plt.Line2D([0], [0], color=common.GREEN_COLOR, linewidth=2.2, label=r"Model, $w_{org}^{salt-free}=0.4$"),
        plt.Line2D(
            [0],
            [0],
            marker="s",
            linestyle="None",
            markerfacecolor="none",
            markeredgecolor=common.GREEN_COLOR,
            color=common.GREEN_COLOR,
            label=r"$w_{org}^{salt-free}=0.4$",
        ),
    ]
    fig.legend(handles=handles, loc="upper center", ncol=2, bbox_to_anchor=(0.5, 0.965), fontsize=9)
    fig.suptitle(
        "Molality-based salt mean ionic activity coefficients $\\gamma_{\\pm}^{m,*}$ in aqueous-organic solvent mixtures",
        fontsize=12,
        y=0.995,
    )
    fig.subplots_adjust(left=0.08, right=0.98, bottom=0.08, top=0.86, wspace=0.20, hspace=0.28)
    common.save_figure(fig, OUTPUT)
    plt.close(fig)

    for cfg, out in zip(PANELS, OUTPUTS):
        common.save_panel_figure(
            lambda ax, cfg=cfg: _plot_panel(ax, *cfg, include_legend=True), out, figsize=(4.2, 3.8)
        )

    common.save_panel_figure(
        lambda ax: _plot_40wt_solvent(ax, "water-methanol", "methanol", ["NaBr", "NaCl"], 3.3, 1.0),
        OUTPUT_40_METHANOL,
        figsize=(5.1, 3.9),
    )
    common.save_panel_figure(
        lambda ax: _plot_40wt_solvent(ax, "water-ethanol", "ethanol", ["NaBr", "NaCl"], 4.4, 1.0),
        OUTPUT_40_ETHANOL,
        figsize=(5.1, 3.9),
    )

    fig40, axes40 = plt.subplots(1, 2, figsize=(10.2, 4.1))
    _plot_40wt_solvent(axes40[0], "water-methanol", "methanol", ["NaBr", "NaCl"], 3.3, 1.0)
    _plot_40wt_solvent(axes40[1], "water-ethanol", "ethanol", ["NaBr", "NaCl"], 4.4, 1.0)
    fig40.suptitle(
        "Molality-based salt mean ionic activity coefficients at 40 wt% organic solvent\nin mixed organic-aqueous solvent systems at 298.15 K and 1 bar.",
        fontsize=11,
        y=0.995,
    )
    fig40.subplots_adjust(left=0.08, right=0.98, bottom=0.16, top=0.80, wspace=0.22)
    common.save_figure(fig40, OUTPUT_40_COMBINED)
    plt.close(fig40)


if __name__ == "__main__":
    main()
