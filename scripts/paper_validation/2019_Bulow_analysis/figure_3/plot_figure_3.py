from __future__ import annotations

from pathlib import Path
import sys

import matplotlib
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from _common import scan_temperature_branch
from scripts.plot_outputs import paper_validation_path, save_plot_figure

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def _plot_branch(ax, scan: dict[str, np.ndarray], *, color: str, label: str, linestyle: str = "-") -> None:
    if scan["T"].size == 0:
        return
    ax.plot(scan["x_il_water_rich"], scan["T"], color=color, linestyle=linestyle, linewidth=1.8)
    ax.plot(scan["x_il_il_rich"], scan["T"], color=color, linestyle=linestyle, linewidth=1.8, label=label)


def main() -> None:
    t = np.linspace(288.15, 360.05, 10)
    epc = scan_temperature_branch("C4mim+", "NTf2-", t, use_kij=False, model_mode="epc")
    water80 = scan_temperature_branch("C4mim+", "NTf2-", t, use_kij=False, model_mode="orig_water")
    il11 = scan_temperature_branch("C4mim+", "NTf2-", t, use_kij=False, model_mode="orig_il")

    fig, ax = plt.subplots(figsize=(6.7, 5.0))
    _plot_branch(ax, water80, color="tab:blue", label=r"original ePC-SAFT ($\varepsilon_r=80$)")
    _plot_branch(ax, il11, color="tab:orange", label=r"original ePC-SAFT ($\varepsilon_r=11$)")
    if water80["T"].size and il11["T"].size:
        ax.plot(water80["x_il_water_rich"], water80["T"], color="0.5", linewidth=1.8)
        ax.plot(il11["x_il_il_rich"], il11["T"], color="0.5", linewidth=1.8, label=r"phase-dependent $\varepsilon$ approx.")
    _plot_branch(ax, epc, color="green", label="ePC-SAFT")
    ax.set_xlabel(r"IL mole fraction, $x_{IL}$")
    ax.set_ylabel("temperature / K")
    ax.set_title("Bulow 2019 Figure 3 style: water + [C4mim][NTf2]")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8)
    out = paper_validation_path(__file__, "figure_3.png")
    fig.tight_layout()
    save_plot_figure(fig, out, dpi=220, bbox_inches=None)
    plt.close(fig)
    print(f"Wrote: {out}")


if __name__ == "__main__":
    main()
