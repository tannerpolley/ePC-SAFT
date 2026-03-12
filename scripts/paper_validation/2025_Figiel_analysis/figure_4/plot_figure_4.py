from __future__ import annotations

from pathlib import Path
import sys

import numpy as np
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import _common as common

OUTPUT = Path(__file__).with_name("figure_4.png")


def _safe_model(dataset: str, ion: str) -> float:
    try:
        return -common.gsolv_ion(dataset, ion, "water", {"water": 1.0})
    except Exception as exc:
        print(f"[figure_4] skipping {dataset} {ion}: {exc}")
        return float("nan")


def _bar_panel(ax, ions, title, ylim):
    literature = common.literature_gsolv_water()
    x = np.arange(len(ions), dtype=float)
    width = 0.22
    lit_vals = np.array([-literature.get(ion, np.nan) for ion in ions], dtype=float)
    figiel_vals = np.array([_safe_model("figiel_2025", ion) for ion in ions], dtype=float)
    bulow_vals = np.array([_safe_model("bulow_2020", ion) for ion in ions], dtype=float)

    mask_lit = np.isfinite(lit_vals)
    mask_figiel = np.isfinite(figiel_vals)
    mask_bulow = np.isfinite(bulow_vals)

    ax.bar(x[mask_lit] - width, lit_vals[mask_lit], width=width, color=common.LIGHT_GRAY, edgecolor="black", linewidth=0.8, label="Literature")
    ax.bar(x[mask_figiel], figiel_vals[mask_figiel], width=width, color=common.BLUE_COLOR, edgecolor="black", linewidth=0.8, label="This work")
    ax.bar(x[mask_bulow] + width, bulow_vals[mask_bulow], width=width, color=common.BROWN_COLOR, edgecolor="black", linewidth=0.8, label="ePC-SAFT advanced")

    labels = [r"$H^+$", r"$Li^+$", r"$Na^+$", r"$K^+$"] if ions[0].endswith("+") else [r"$Cl^-$", r"$Br^-$", r"$I^-$"]
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_title(title, fontsize=10)
    ax.set_ylim(*ylim)
    ax.set_ylabel(r"$-\Delta G_i^{solv,\infty,x}$ / kJ mol$^{-1}$")


def main() -> None:
    common.configure_style()
    fig, axes = plt.subplots(1, 2, figsize=(7.4, 3.4))
    _bar_panel(axes[0], ["H+", "Li+", "Na+", "K+"], "Cations", (0.0, 1200.0))
    _bar_panel(axes[1], ["Cl-", "Br-", "I-"], "Anions", (0.0, 800.0))
    common.panel_label(axes[0], "a)")
    common.panel_label(axes[1], "b)")
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=3, bbox_to_anchor=(0.5, 1.05), fontsize=9)
    common.save_figure(fig, OUTPUT)
    plt.close(fig)


if __name__ == "__main__":
    main()
