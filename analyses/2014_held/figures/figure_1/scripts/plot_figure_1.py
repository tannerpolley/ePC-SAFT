from __future__ import annotations

from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np

ANALYSIS_SCRIPTS = Path(__file__).resolve().parents[3] / "scripts"
if str(ANALYSIS_SCRIPTS) not in sys.path:
    sys.path.insert(0, str(ANALYSIS_SCRIPTS))

import _common as common

OUTPUT = Path(__file__).with_name("figure_1.png")
INPUT_DIR = Path(__file__).resolve().parents[1] / "input"


def main() -> None:
    common.configure_style()

    m_plot = np.linspace(0.01, 4.0, 101)

    data = {}
    for salt in ("NaCl", "KBr"):
        m_exp, phi_exp = common.load_osmotic_data(INPUT_DIR / f"{salt}.csv")
        data[salt] = {
            "m_exp": m_exp,
            "phi_exp": phi_exp,
            "phi_2008": common.calc_osmotic_curve(salt, m_plot, "2008"),
            "phi_2014": common.calc_osmotic_curve(salt, m_plot, "2014"),
        }

    fig, ax = plt.subplots(figsize=(7.8, 4.9))
    marker_map = {"NaCl": "o", "KBr": "s"}
    line_style = {"NaCl": "-", "KBr": "--"}

    for salt in ("NaCl", "KBr"):
        ax.scatter(
            data[salt]["m_exp"],
            data[salt]["phi_exp"],
            marker=marker_map[salt],
            s=34,
            facecolors="none",
            edgecolors="black",
            linewidths=1.0,
            label=f"{salt} data",
            zorder=5,
        )
        ax.plot(
            m_plot,
            data[salt]["phi_2008"],
            color="0.5",
            linewidth=1.7,
            linestyle=line_style[salt],
            label=f"{salt} strategy 1 (2008)",
            zorder=2,
        )
        ax.plot(
            m_plot,
            data[salt]["phi_2014"],
            color="black",
            linewidth=2.0,
            linestyle=line_style[salt],
            label=f"{salt} strategy 2 (2014)",
            zorder=3,
        )

    ax.set_xlim(0.0, 4.0)
    ax.set_ylim(0.8, 1.2)
    ax.set_xlabel(r"molality, $m$ / mol kg$^{-1}$")
    ax.set_ylabel(r"molal osmotic coefficient, $\phi_m$")
    ax.set_title("Held 2014 Fig. 1 reproduction (NaCl and KBr, 298.15 K)")
    ax.grid(True, alpha=0.25)
    ax.legend(fontsize=8, ncol=2)

    common.save_figure(fig, OUTPUT)
    plt.close(fig)


if __name__ == "__main__":
    main()


