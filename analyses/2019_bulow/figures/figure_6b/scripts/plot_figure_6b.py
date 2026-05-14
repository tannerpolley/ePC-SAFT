from __future__ import annotations

from pathlib import Path
import sys

import matplotlib
import numpy as np

ANALYSIS_SCRIPTS = Path(__file__).resolve().parents[3] / "scripts"
if str(ANALYSIS_SCRIPTS) not in sys.path:
    sys.path.insert(0, str(ANALYSIS_SCRIPTS))

from _common import ANION_SERIES, scan_temperature_branch
from scripts.plot_outputs import paper_validation_path, save_plot_figure

matplotlib.use("Agg")
import matplotlib.pyplot as plt

COLORS = ["0.5", "tab:orange", "black", "green"]


def main() -> None:
    t = np.linspace(278.15, 364.55, 8)
    fig, ax = plt.subplots(figsize=(6.8, 5.1))
    for anion, color in zip(ANION_SERIES, COLORS):
        scan = scan_temperature_branch("C8mim+", anion, t, use_kij=True, model_mode="epc")
        if scan["T"].size == 0:
            continue
        label = anion.replace("-", "")
        ax.plot(scan["x_il_il_rich"], scan["T"], color=color, linewidth=1.9, label=label)
    ax.set_xlabel(r"IL mole fraction, $x_{IL}$")
    ax.set_ylabel("temperature / K")
    ax.set_title("Bulow 2019 Figure 6b style")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8, title="[C8mim][X]")
    out = paper_validation_path(__file__, "figure_6b.png")
    fig.tight_layout()
    save_plot_figure(fig, out, dpi=220, bbox_inches=None)
    plt.close(fig)
    print(f"Wrote: {out}")


if __name__ == "__main__":
    main()

