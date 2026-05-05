from __future__ import annotations

from pathlib import Path
import sys

import matplotlib
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from _common import IMIDAZOLIUM_NTF2, scan_temperature_branch
from scripts.plot_outputs import paper_validation_path, save_plot_figure

matplotlib.use("Agg")
import matplotlib.pyplot as plt

COLORS = ["black", "0.5", "tab:orange", "green"]


def main() -> None:
    t = np.linspace(288.15, 333.15, 8)
    fig, ax = plt.subplots(figsize=(6.8, 5.1))
    for cation, color in zip(IMIDAZOLIUM_NTF2, COLORS):
        scan = scan_temperature_branch(cation, "NTf2-", t, use_kij=True, model_mode="epc")
        if scan["T"].size == 0:
            continue
        label = cation.replace("+", "")
        ax.plot(scan["x_il_il_rich"], scan["T"], color=color, linewidth=1.9, label=label)
        if cation == "C8mim+":
            ax.plot(scan["x_il_water_rich"], scan["T"], color=color, linewidth=1.2, linestyle="--")
    ax.set_xlabel(r"IL mole fraction, $x_{IL}$")
    ax.set_ylabel("temperature / K")
    ax.set_title("Bulow 2019 Figure 6a style")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8, title="[Cnmim][NTf2]")
    out = paper_validation_path(__file__, "figure_6a.png")
    fig.tight_layout()
    save_plot_figure(fig, out, dpi=220, bbox_inches=None)
    plt.close(fig)
    print(f"Wrote: {out}")


if __name__ == "__main__":
    main()
