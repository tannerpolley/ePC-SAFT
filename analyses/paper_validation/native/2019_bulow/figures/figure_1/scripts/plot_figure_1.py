from __future__ import annotations

from pathlib import Path
import sys as _bootstrap_sys
from pathlib import Path as _BootstrapPath

for _candidate in _BootstrapPath(__file__).resolve().parents:
    if (_candidate / "scripts" / "plot_outputs.py").is_file():
        if str(_candidate) not in _bootstrap_sys.path:
            _bootstrap_sys.path.insert(0, str(_candidate))
        break
else:
    raise ModuleNotFoundError("Could not locate repo root containing scripts/plot_outputs.py")
from scripts.plot_outputs import analysis_scripts_dir
import sys



import matplotlib
import numpy as np

ANALYSIS_SCRIPTS = analysis_scripts_dir(__file__)
if str(ANALYSIS_SCRIPTS) not in sys.path:
    sys.path.insert(0, str(ANALYSIS_SCRIPTS))

from _common import IL_DIELC, WATER_DIELC
from scripts.plot_outputs import paper_validation_path, save_plot_figure

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
    out = paper_validation_path(__file__, "figure_1.png")
    fig.tight_layout()
    save_plot_figure(fig, out, dpi=220, bbox_inches=None)
    plt.close(fig)
    print(f"Wrote: {out}")


if __name__ == "__main__":
    main()

