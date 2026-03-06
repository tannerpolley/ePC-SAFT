from __future__ import annotations

from pathlib import Path
import sys

import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import _common as common

OUTPUT = Path(__file__).with_name("figure_5.png")
DATA_ROOT = common.REPO_ROOT / "data" / "MIAC" / "water"

SERIES = {
    "Li": {"color": common.ORGANIC_COLOR, "marker": "^"},
    "Na": {"color": common.GREEN_COLOR, "marker": "s"},
    "K": {"color": common.GRAY_COLOR, "marker": "o"},
}
PANELS = [
    ("a)", ["LiCl", "NaCl", "KCl"], "Cl$^-$ salts"),
    ("b)", ["LiBr", "NaBr", "KBr"], "Br$^-$ salts"),
]


def main() -> None:
    common.configure_style()
    fig, axes = plt.subplots(1, 2, figsize=(7.4, 3.35), sharey=True)
    for ax, (label, salts, title) in zip(axes, PANELS):
        for salt in salts:
            cation = salt.split("C")[0] if "Cl" in salt else salt.split("B")[0]
            style = SERIES[cation]
            data = common.read_miac_dataset(DATA_ROOT / f"water-{salt}.csv", "water")
            x_data = [row["molality"] for row in data]
            y_data = [row["miac_m"] for row in data]
            m_grid, y_model = common.mean_ionic_activity_curve("figiel_2025", salt, "water", {"water": 1.0}, 6.0, points=600)
            ax.plot(m_grid, y_model, color=style["color"], linewidth=1.8)
            ax.scatter(x_data, y_data, marker=style["marker"], s=26, facecolor="none", edgecolor=style["color"], linewidth=1.0, label=salt)
        common.panel_label(ax, label)
        ax.set_title(title, fontsize=10)
        ax.set_xlim(0.0, 6.0)
        ax.set_ylim(0.4, 4.6)
        ax.set_xlabel(r"$\bar{m}_{salt}$ / mol kg$^{-1}$")
    axes[0].set_ylabel(r"$\gamma_{\pm}^{m,*}$ / -")
    axes[0].legend(loc="upper left", fontsize=8)
    axes[1].legend(loc="upper left", fontsize=8)
    common.save_figure(fig, OUTPUT)
    plt.close(fig)


if __name__ == "__main__":
    main()
