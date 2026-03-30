from __future__ import annotations

from pathlib import Path
import sys

import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import _common as common

OUTPUT = Path(__file__).with_name("figure_2d.png")


PANELS = [
    {
        "label": "a)",
        "solvent_system": "ethanol",
        "salt": "NaI",
        "data_path": common.REPO_ROOT / "data" / "MIAC" / "ethanol" / "ethanol-NaI.csv",
        "title": "NaI in ethanol",
        "color": "#b55d09",
    },
    {
        "label": "b)",
        "solvent_system": "methanol",
        "salt": "NaI",
        "data_path": common.REPO_ROOT / "data" / "MIAC" / "methanol" / "methanol-NaI.csv",
        "title": "NaI in methanol",
        "color": "#17632d",
    },
]


def _plot_panel(ax, panel: dict) -> None:
    entries = common.read_miac_dataset(panel["data_path"], panel["solvent_system"])
    if not entries:
        raise ValueError(f"No usable rows in {panel['data_path']}")

    m_data = [float(row["molality"]) for row in entries]
    y_data = [float(row["miac_m"]) for row in entries]
    m_max = max(m_data)

    grid, y_model = common.mean_ionic_activity_curve(
        dataset="2012_Held",
        salt=panel["salt"],
        solvent_system=panel["solvent_system"],
        comp={panel["solvent_system"]: 1.0},
        m_max=m_max,
        points=220,
    )

    ax.plot(grid, y_model, color=panel["color"], linewidth=2.0, label="ePC-SAFT (2012_Held)")
    ax.scatter(
        m_data,
        y_data,
        s=24,
        marker="o",
        facecolors="none",
        edgecolors="black",
        linewidths=0.9,
        label="experimental",
        zorder=5,
    )

    ax.set_xlim(0.0, 2.0)
    ax.set_ylim(0.0, 1.0)
    ax.set_xlabel(r"molality, $m$ / mol kg$^{-1}$")
    ax.set_title(panel["title"], fontsize=10)
    ax.grid(True, alpha=0.25)
    ax.text(0.03, 0.96, panel["label"], transform=ax.transAxes, ha="left", va="top", fontsize=11)


def main() -> None:
    common.configure_style()

    fig, axes = plt.subplots(1, 2, figsize=(7.8, 3.5), sharey=True)
    for ax, panel in zip(axes, PANELS):
        _plot_panel(ax, panel)

    axes[0].set_ylabel(r"$\gamma_{\pm}^{m,*}$ / -")
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=2, bbox_to_anchor=(0.5, 1.07), fontsize=9)

    common.save_figure(fig, OUTPUT)
    plt.close(fig)


if __name__ == "__main__":
    main()

