from __future__ import annotations

from pathlib import Path
import sys

import matplotlib
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from _common import IMIDAZOLIUM_NTF2, water_solubility_in_il
from scripts.plot_outputs import paper_validation_path, save_plot_figure

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def main() -> None:
    labels = [cation.replace("+", "") for cation in IMIDAZOLIUM_NTF2]
    epc = [water_solubility_in_il(cation, "NTf2-", use_kij=False, model_mode="epc") for cation in IMIDAZOLIUM_NTF2]
    orig = [water_solubility_in_il(cation, "NTf2-", use_kij=False, model_mode="orig_water") for cation in IMIDAZOLIUM_NTF2]
    xpos = np.arange(len(labels))

    fig, ax = plt.subplots(figsize=(7.0, 4.8))
    ax.bar(xpos, epc, width=0.55, color="green", alpha=0.85, label="ePC-SAFT prediction")
    ax.plot(xpos, orig, color="tab:blue", linewidth=2.0, marker="o", label=r"original ePC-SAFT ($\varepsilon_r=80$)")
    ax.set_xticks(xpos, labels)
    ax.set_ylabel(r"water mole fraction in IL-rich phase, $x_w$")
    ax.set_title("Bulow 2019 Figure 4 style")
    ax.grid(True, axis="y", alpha=0.3)
    ax.legend(fontsize=8)
    out = paper_validation_path(__file__, "figure_4.png")
    fig.tight_layout()
    save_plot_figure(fig, out, dpi=220, bbox_inches=None)
    plt.close(fig)
    print(f"Wrote: {out}")


if __name__ == "__main__":
    main()
