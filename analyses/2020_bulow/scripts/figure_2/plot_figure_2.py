from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

SCRIPT_DIR = Path(__file__).resolve().parent
ANALYSIS_ROOT = SCRIPT_DIR.parent
REPO_ROOT = ANALYSIS_ROOT.parents[2]

if str(ANALYSIS_ROOT) not in sys.path:
    sys.path.insert(0, str(ANALYSIS_ROOT))
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts._env import require_epcsaft_install

require_epcsaft_install()

import _plot_common as common
from epcsaft.parameters import get_prop_dict
from scripts._epcsaft_oop import epcsaft_density, epcsaft_solvation_free_energy

DATA_PATH = SCRIPT_DIR / "data" / "water_comparisons.csv"
OUTPUT_PATH = SCRIPT_DIR / "figure_2.png"
T_REF = 298.15
P_REF = 1.0e5

SERIES = [
    ("data", "Literature data", "#bdbdbd", None),
    ("SAFT-VR", "SAFT-VR", "#8e63c7", None),
    ("advanced", "ePC-SAFT advanced (paper)", "#2ca02c", None),
    ("advanced_calc", "ePC-SAFT advanced (epcsaft)", "#2ca02c", "////"),
    ("revised", "ePC-SAFT revised (paper)", "#e67e22", None),
    ("revised_calc", "ePC-SAFT revised (epcsaft)", "#e67e22", "////"),
]

ADVANCED_CASES = {
    "Li+": ("2020_Bulow", ["Li+", "Cl-", "Water"]),
    "Na+": ("2020_Bulow", ["Na+", "Cl-", "Water"]),
    "K+": ("2020_Bulow", ["K+", "Cl-", "Water"]),
    "F-": ("2020_Bulow", ["Na+", "F-", "Water"]),
    "Cl-": ("2020_Bulow", ["Na+", "Cl-", "Water"]),
    "Br-": ("2020_Bulow", ["Na+", "Br-", "Water"]),
    "I-": ("2020_Bulow", ["Na+", "I-", "Water"]),
}

REVISED_CASES = {
    "Li+": ("2014_Held", ["Li+", "Cl-", "Water"]),
    "Na+": ("2014_Held", ["Na+", "Cl-", "Water"]),
    "K+": ("2014_Held", ["K+", "Cl-", "Water"]),
    "F-": ("2014_Held", ["Na+", "F-", "Water"]),
    "Cl-": ("2014_Held", ["Na+", "Cl-", "Water"]),
    "Br-": ("2014_Held", ["Na+", "Br-", "Water"]),
    "I-": ("2014_Held", ["Na+", "I-", "Water"]),
}


def _compute_gsolv(dataset_name: str, species: list[str], ion: str) -> float:
    x = np.asarray([1.0e-8, 1.0e-8, 1.0 - 2.0e-8], dtype=float)
    params = get_prop_dict(dataset_name, species, x, T_REF, user_options={})
    rho = epcsaft_density(T_REF, P_REF, x, params, phase="liq")
    values = epcsaft_solvation_free_energy(T_REF, rho, x, params, species=species)
    return float(values[ion]) / 1000.0


def _computed_row(ions: list[str], case_map: dict[str, tuple[str, list[str]]]) -> np.ndarray:
    out = np.full(len(ions), np.nan, dtype=float)
    for idx, ion in enumerate(ions):
        case = case_map.get(ion)
        if case is None:
            continue
        dataset_name, species = case
        out[idx] = _compute_gsolv(dataset_name, species, ion)
    return out


def main() -> None:
    common.configure_style()
    frame = common.load_indexed_csv(DATA_PATH)
    ions = list(frame.columns)

    row_map = {
        "data": frame.values("data", ions),
        "SAFT-VR": frame.values("SAFT-VR", ions),
        "advanced": frame.values("advanced", ions),
        "advanced_calc": _computed_row(ions, ADVANCED_CASES),
        "revised": frame.values("revised", ions),
        "revised_calc": _computed_row(ions, REVISED_CASES),
    }

    x = np.arange(len(ions), dtype=float)
    width = 0.13
    offsets = np.linspace(-2.5 * width, 2.5 * width, len(SERIES))

    fig, ax = plt.subplots(figsize=(12.8, 6.4))
    for offset, (row_key, label, color, hatch) in zip(offsets, SERIES):
        values = row_map[row_key]
        valid = np.isfinite(values)
        if not np.any(valid):
            continue
        bars = ax.bar(
            x[valid] + offset,
            values[valid],
            width=width,
            label=label,
            color=color,
            edgecolor="black",
            linewidth=0.45,
            hatch=hatch,
            alpha=0.9 if hatch is None else 0.75,
        )
        common.annotate_bar_values(ax, bars, fontsize=6)

    ax.axhline(0.0, color="black", linewidth=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(ions)
    ax.set_ylabel(r"$\Delta G_{\mathrm{hyd},i}^{\infty}$ / kJ mol$^{-1}$")
    ax.set_title("Bulow 2020 Part I Figure 2: Gibbs Energy of Hydration with epcsaft Overlay")
    ax.grid(axis="y", alpha=0.22)

    values = np.vstack(list(row_map.values()))
    y_min = float(np.nanmin(values))
    y_max = float(np.nanmax(values))
    ax.set_ylim(y_min - 50.0, max(40.0, y_max + 55.0))

    common.annotate_percent_deltas(
        ax,
        x + offsets[3],
        row_map["advanced"],
        row_map["advanced_calc"],
        xs_ref=x + offsets[2],
    )
    common.annotate_percent_deltas(
        ax,
        x + offsets[5],
        row_map["revised"],
        row_map["revised_calc"],
        xs_ref=x + offsets[4],
    )
    common.add_percent_note(ax)

    ax.legend(ncol=3, frameon=True)
    common.save_figure(fig, OUTPUT_PATH)
    plt.close(fig)


if __name__ == "__main__":
    main()
