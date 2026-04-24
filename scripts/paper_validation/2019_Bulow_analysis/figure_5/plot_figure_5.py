from __future__ import annotations

from pathlib import Path
import sys

import matplotlib
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from _common import scan_temperature_branch
from scripts.plot_outputs import paper_validation_path

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def _plot(ax, scan: dict[str, np.ndarray], color: str, linestyle: str, label: str) -> None:
    if scan["T"].size == 0:
        return
    ax.plot(scan["x_il_water_rich"], scan["T"], color=color, linestyle=linestyle, linewidth=1.8)
    ax.plot(scan["x_il_il_rich"], scan["T"], color=color, linestyle=linestyle, linewidth=1.8, label=label)


def main() -> None:
    t = np.linspace(288.15, 360.05, 10)
    orig_0 = scan_temperature_branch("C4mim+", "NTf2-", t, use_kij=False, model_mode="orig_water")
    orig_k = scan_temperature_branch("C4mim+", "NTf2-", t, use_kij=True, model_mode="orig_water")
    epc_0 = scan_temperature_branch("C4mim+", "NTf2-", t, use_kij=False, model_mode="epc")
    epc_k = scan_temperature_branch("C4mim+", "NTf2-", t, use_kij=True, model_mode="epc")

    fig, ax = plt.subplots(figsize=(6.7, 5.0))
    _plot(ax, orig_0, "tab:blue", "-", r"original ePC-SAFT, $k_{ij}=0$")
    _plot(ax, orig_k, "tab:blue", "--", r"original ePC-SAFT, fitted $k_{ij}$")
    _plot(ax, epc_0, "green", "-", r"ePC-SAFT, $k_{ij}=0$")
    _plot(ax, epc_k, "green", "--", r"ePC-SAFT, fitted $k_{ij}$")
    ax.set_xlabel(r"IL mole fraction, $x_{IL}$")
    ax.set_ylabel("temperature / K")
    ax.set_title("Bulow 2019 Figure 5 style: water + [C4mim][NTf2]")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8)
    out = paper_validation_path(__file__, "figure_5.png")
    fig.tight_layout()
    fig.savefig(out, dpi=220)
    plt.close(fig)
    print(f"Wrote: {out}")


if __name__ == "__main__":
    main()
