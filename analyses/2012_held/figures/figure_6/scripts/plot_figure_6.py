from __future__ import annotations

from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np

ANALYSIS_SCRIPTS = Path(__file__).resolve().parents[3] / "scripts"
if str(ANALYSIS_SCRIPTS) not in sys.path:
    sys.path.insert(0, str(ANALYSIS_SCRIPTS))

import _common as common

OUTPUT = Path(__file__).with_name("figure_6.png")
DATA_PATH = common.analysis_data_path(__file__, "water-methanol-NaBr.csv", kind="input")
SERIES = [
    {
        "w_methanol": 0.6,
        "marker": "o",
        "s": 46,
        "facecolors": "white",
        "edgecolors": "#4c4c4c",
        "linecolor": "#2b2b2b",
        "label": r"$x_{\mathrm{H_2O}}^0 = 0.542,\ x_{\mathrm{MeOH}}^0 = 0.458$",
    },
    {
        "w_methanol": 0.8,
        "marker": "s",
        "s": 34,
        "facecolors": "#cfcfcf",
        "edgecolors": "#686868",
        "linecolor": "#5d5d5d",
        "label": r"$x_{\mathrm{H_2O}}^0 = 0.308,\ x_{\mathrm{MeOH}}^0 = 0.692$",
    },
]


def _weight_to_mole_comp(weights: dict[str, float]) -> dict[str, float]:
    moles = {}
    for solvent, weight in weights.items():
        moles[solvent] = float(weight) / common.SOLVENT_MW[solvent]
    total = sum(moles.values())
    if total <= 0.0:
        raise ValueError(f"Invalid solvent weights: {weights}")
    return {solvent: moles[solvent] / total for solvent in moles}


def _load_grouped_data() -> dict[float, list[dict[str, float]]]:
    _, rows = common.read_csv_rows(DATA_PATH)
    grouped: dict[float, list[dict[str, float]]] = {}
    for row in rows:
        w_methanol = common.parse_float(row.get("x_Methanol")) or common.parse_float(row.get("w_Methanol"))
        molality = common.parse_float(row.get("molality"))
        gamma = (
            common.parse_float(row.get("miac_m"))
            or common.parse_float(row.get("gamma"))
            or common.parse_float(row.get("miac"))
        )
        if w_methanol is None or molality is None or gamma is None:
            continue
        key = round(float(w_methanol), 6)
        grouped.setdefault(key, []).append({"molality": float(molality), "miac_m": float(gamma)})
    for values in grouped.values():
        values.sort(key=lambda item: item["molality"])
    return grouped


def _plot_binary_reference(ax, solvent: str, color: str, label: str, m_max: float = 3.0) -> None:
    grid, y_model = common.mean_ionic_activity_curve(
        dataset="2012_Held",
        salt="NaBr",
        solvent_system=solvent,
        comp={solvent: 1.0},
        m_max=m_max,
        points=500,
    )
    ax.plot(grid, y_model, color=color, linewidth=1.5, zorder=1, label=label)


def main() -> None:
    common.configure_style()

    # Ye 1994 Table I reports the solvent composition as methanol weight fraction.
    grouped = _load_grouped_data()

    fig, ax = plt.subplots(figsize=(7.6, 4.8))
    _plot_binary_reference(ax, "water", color="#d2d2d2", label="binary water / NaBr")
    _plot_binary_reference(ax, "methanol", color="#b5b5b5", label="binary methanol / NaBr", m_max=1.6)

    for series in SERIES:
        rows = grouped[round(series["w_methanol"], 6)]
        comp = _weight_to_mole_comp({"water": 1.0 - series["w_methanol"], "methanol": series["w_methanol"]})

        m_data = np.asarray([row["molality"] for row in rows], dtype=float)
        y_data = np.asarray([row["miac_m"] for row in rows], dtype=float)
        mask = m_data <= 3.0 + 1e-12

        grid, y_model = common.mean_ionic_activity_curve(
            dataset="2012_Held",
            salt="NaBr",
            solvent_system="water-methanol",
            comp=comp,
            m_max=3.0,
            points=600,
        )

        ax.plot(grid, y_model, color=series["linecolor"], linewidth=2.0, zorder=3)
        ax.scatter(
            m_data[mask],
            y_data[mask],
            s=series["s"],
            marker=series["marker"],
            facecolors=series["facecolors"],
            edgecolors=series["edgecolors"],
            linewidths=1.0,
            zorder=5,
            label=series["label"],
        )

    ax.set_xlim(0.0, 3.0)
    ax.set_ylim(0.3, 1.0)
    ax.set_xlabel(r"molality, $m_{\mathrm{NaBr}}$ / mol kg$^{-1}$")
    ax.set_ylabel(r"$\gamma_{\pm}^{m,*}$ / -")
    ax.set_title("2012 Fig. 6 style: NaBr in water/methanol mixtures (298.15 K)")
    ax.legend(loc="upper right", fontsize=8.6, handlelength=2.4)
    ax.grid(False)

    common.save_figure(fig, OUTPUT)
    plt.close(fig)


if __name__ == "__main__":
    main()

