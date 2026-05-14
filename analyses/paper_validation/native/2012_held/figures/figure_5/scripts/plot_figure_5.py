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



import matplotlib.pyplot as plt
import numpy as np

ANALYSIS_SCRIPTS = analysis_scripts_dir(__file__)
if str(ANALYSIS_SCRIPTS) not in sys.path:
    sys.path.insert(0, str(ANALYSIS_SCRIPTS))

import _common as common

OUTPUT = Path(__file__).with_name("figure_5.png")
DATA_PATH = common.analysis_data_path(__file__, "esteso_1989_table1_selected_series.csv", kind="input")
SERIES = [
    {
        "xw": 0.9109,
        "marker": "o",
        "s": 46,
        "facecolors": "#d0d0d0",
        "edgecolors": "#8f8f8f",
        "label_pos": (2.12, 0.615),
    },
    {
        "xw": 0.7932,
        "marker": "s",
        "s": 32,
        "facecolors": "#a8a8a8",
        "edgecolors": "#747474",
        "label_pos": (1.72, 0.470),
    },
    {
        "xw": 0.6303,
        "marker": "*",
        "s": 80,
        "facecolors": "#ececec",
        "edgecolors": "#8a8a8a",
        "label_pos": (1.05, 0.345),
    },
    {
        "xw": 0.3899,
        "marker": "D",
        "s": 34,
        "facecolors": "none",
        "edgecolors": "#5f5f5f",
        "label_pos": (0.12, 0.335),
    },
]


def _find_signature(grouped: dict, xw_target: float):
    for signature in grouped:
        water_x = round(float(dict(signature).get("water", 0.0)), 4)
        if abs(water_x - xw_target) <= 5e-4:
            return signature
    raise KeyError(f"Could not resolve x_water={xw_target:.4f} in {DATA_PATH}.")


def _plot_pure_reference(ax, solvent_system: str, color: str, lw: float) -> None:
    grid, y_model = common.mean_ionic_activity_curve(
        dataset="2012_Held",
        salt="NaCl",
        solvent_system=solvent_system,
        comp={solvent_system: 1.0},
        m_max=3.0,
        points=500,
    )
    ax.plot(grid, y_model, color=color, linewidth=lw, zorder=1)


def main() -> None:
    common.configure_style()

    entries = common.read_miac_dataset(DATA_PATH, "water-ethanol")
    grouped = common.group_by_signature(entries)

    fig, ax = plt.subplots(figsize=(7.6, 4.8))

    # Binary water/NaCl and ethanol/NaCl references span the full panel width.
    _plot_pure_reference(ax, "water", color="#c2c2c2", lw=1.7)
    _plot_pure_reference(ax, "ethanol", color="#a4a4a4", lw=1.5)

    for series in SERIES:
        signature = _find_signature(grouped, series["xw"])
        rows = grouped[signature]
        comp = {key: float(val) for key, val in signature}

        m_data = np.asarray([float(row["molality"]) for row in rows], dtype=float)
        y_data = np.asarray([float(row["miac_m"]) for row in rows], dtype=float)

        grid, y_model = common.mean_ionic_activity_curve(
            dataset="2012_Held",
            salt="NaCl",
            solvent_system="water-ethanol",
            comp=comp,
            m_max=float(np.max(m_data)),
            points=360,
        )

        ax.plot(grid, y_model, color="black", linewidth=2.0, zorder=3)
        scatter_kwargs = {
            "s": series["s"],
            "marker": series["marker"],
            "linewidths": 1.0,
            "zorder": 5,
        }
        scatter_kwargs["facecolors"] = series["facecolors"]
        scatter_kwargs["edgecolors"] = series["edgecolors"]
        ax.scatter(m_data, y_data, **scatter_kwargs)
        ax.text(
            series["label_pos"][0],
            series["label_pos"][1],
            rf"$x_w = {series['xw']:.4f}$",
            fontsize=10,
            color="#3d3d3d",
        )

    ax.set_xlim(0.0, 3.0)
    ax.set_ylim(0.3, 1.0)
    ax.set_xlabel(r"molality, $m_{\mathrm{NaCl}}$ / mol kg$^{-1}$")
    ax.set_ylabel(r"$\gamma_{\pm}^{m,*}$ / -")
    ax.set_title("2012 Fig. 5 style: NaCl in water/ethanol mixtures (298.15 K)")
    ax.grid(False)

    common.save_figure(fig, OUTPUT)
    plt.close(fig)


if __name__ == "__main__":
    main()

