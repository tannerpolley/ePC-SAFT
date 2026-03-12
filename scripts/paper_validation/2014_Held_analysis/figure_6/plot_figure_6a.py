from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

import _shared

ROOT = Path(__file__).resolve().parents[1]

import _common as common

OUTPUT = Path(__file__).with_name("figure_6a.png")


def main() -> None:
    common.configure_style()
    exp_rows = _shared.load_experimental_rows()
    model_rows = _shared.solve_model_rows()

    fig, ax = plt.subplots(figsize=(7.4, 4.9))

    x_exp = np.asarray([row["w_nh4cl_aq"] * 100.0 for row in exp_rows], dtype=float)
    y_org_exp = np.asarray([row["w_buoh_org"] * 100.0 for row in exp_rows], dtype=float)
    y_aq_exp = np.asarray([row["w_buoh_aq"] * 100.0 for row in exp_rows], dtype=float)
    ax.scatter(x_exp, y_org_exp, marker="o", s=42, color="black", label="1-butanol in organic phase (data)", zorder=7)
    ax.scatter(
        x_exp,
        y_aq_exp,
        marker="s",
        s=38,
        facecolors="0.65",
        edgecolors="0.25",
        linewidths=0.8,
        label="1-butanol in aqueous phase (data)",
        zorder=7,
    )

    if model_rows:
        x_model = np.asarray([row["w_nh4cl_aq"] * 100.0 for row in model_rows], dtype=float)
        y_org_model = np.asarray([row["w_buoh_org"] * 100.0 for row in model_rows], dtype=float)
        y_aq_model = np.asarray([row["w_buoh_aq"] * 100.0 for row in model_rows], dtype=float)
        order = np.argsort(x_model)
        x_model = x_model[order]
        y_org_model = y_org_model[order]
        y_aq_model = y_aq_model[order]
        ax.plot(x_model, y_org_model, color="black", linewidth=2.0, label="Organic phase model")
        ax.plot(x_model, y_aq_model, color="0.45", linewidth=1.8, label="Aqueous phase model")

    ax.set_xlim(0.0, 25.0)
    ax.set_ylim(0.0, 100.0)
    ax.set_xticks(np.arange(0.0, 26.0, 5.0))
    ax.set_xlabel(r"$w_{NH_4Cl}^{aqueous}$ / wt.%")
    ax.set_ylabel(r"$w_{1-butanol}$ / wt.%")
    ax.set_title("Held 2014 Fig. 6a reproduction")
    ax.grid(True, alpha=0.24)
    ax.legend(fontsize=8, ncol=1, loc="center right")

    common.save_figure(fig, OUTPUT)
    plt.close(fig)


if __name__ == "__main__":
    main()
