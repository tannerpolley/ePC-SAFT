from __future__ import annotations

from pathlib import Path
import sys

import matplotlib
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from _common import IL_DIELC, WATER_DIELC

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def main() -> None:
    x_il = np.linspace(0.0, 1.0, 200)
    eps = WATER_DIELC * (1.0 - x_il) + IL_DIELC * x_il

    fig, ax = plt.subplots(figsize=(6.4, 4.8))
    ax.plot(x_il, eps, color="green", linewidth=2.2, label="linear ePC-SAFT mixing rule")
    ax.scatter([0.0, 1.0], [WATER_DIELC, IL_DIELC], color="black", s=28, zorder=4, label="pure-component anchors")
    ax.set_xlabel(r"IL mole fraction, $x_{IL}$")
    ax.set_ylabel(r"relative dielectric constant, $\varepsilon_r$")
    ax.set_title("Bulow 2019 Figure 1 style")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8)
    out = Path(__file__).resolve().parent / "figure_1.png"
    fig.tight_layout()
    fig.savefig(out, dpi=220)
    plt.close(fig)
    print(f"Wrote: {out}")


if __name__ == "__main__":
    main()
