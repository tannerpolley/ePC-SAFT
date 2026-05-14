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

from _common import FIG7_ILS, il_label, water_solubility_in_il
from scripts.plot_outputs import paper_validation_path, save_plot_figure

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def main() -> None:
    labels = [il_label(c, a) for c, a in FIG7_ILS]
    vals = [water_solubility_in_il(c, a, use_kij=True, model_mode="epc") for c, a in FIG7_ILS]
    xpos = np.arange(len(labels))

    fig, ax = plt.subplots(figsize=(12.5, 5.0))
    ax.bar(xpos, vals, color="green", alpha=0.85)
    ax.set_xticks(xpos, labels, rotation=60, ha="right")
    ax.set_ylabel(r"water mole fraction in IL-rich phase, $x_w$")
    ax.set_title("Bulow 2019 Figure 7 style")
    ax.grid(True, axis="y", alpha=0.3)
    out = paper_validation_path(__file__, "figure_7.png")
    fig.tight_layout()
    save_plot_figure(fig, out, dpi=220, bbox_inches=None)
    plt.close(fig)
    print(f"Wrote: {out}")


if __name__ == "__main__":
    main()

