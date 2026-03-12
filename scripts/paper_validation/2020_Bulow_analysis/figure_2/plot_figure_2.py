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

import _plot_common as common
from data.epcsaft_properties import get_prop_dict
from pcsaft import pcsaft_den, pcsaft_gsolv


DATA_PATH = SCRIPT_DIR / "data" / "water_comparisons.csv"
OUTPUT_PATH = SCRIPT_DIR / "figure_2.png"
T_REF = 298.15
P_REF = 1.0e5

SERIES = [
    ("data", "Literature data", "#bdbdbd", None),
    ("SAFT-VR", "SAFT-VR", "#8e63c7", None),
    ("advanced", "ePC-SAFT advanced (paper)", "#2ca02c", None),
    ("advanced_calc", "ePC-SAFT advanced (pcsaft)", "#2ca02c", "////"),
    ("revised", "ePC-SAFT revised (paper)", "#e67e22", None),
    ("revised_calc", "ePC-SAFT revised (pcsaft)", "#e67e22", "////"),
]

ADVANCED_CASES = {
    "Li+": ("bulow_2020", ["Li+", "Cl-", "Water"]),
    "Na+": ("bulow_2020", ["Na+", "Cl-", "Water"]),
    "K+": ("bulow_2020", ["K+", "Cl-", "Water"]),
    "F-": ("bulow_2020", ["Na+", "F-", "Water"]),
    "Cl-": ("bulow_2020", ["Na+", "Cl-", "Water"]),
    "Br-": ("bulow_2020", ["Na+", "Br-", "Water"]),
    "I-": ("bulow_2020", ["Na+", "I-", "Water"]),
}

REVISED_CASES = {
    "Li+": ("held_2014", ["Li+", "Cl-", "Water"]),
    "Na+": ("held_2014", ["Na+", "Cl-", "Water"]),
    "K+": ("held_2014", ["K+", "Cl-", "Water"]),
    "F-": ("held_2014", ["Na+", "F-", "Water"]),
    "Cl-": ("held_2014", ["Na+", "Cl-", "Water"]),
    "Br-": ("held_2014", ["Na+", "Br-", "Water"]),
    "I-": ("held_2014", ["Na+", "I-", "Water"]),
}


def _compute_gsolv(dataset_name: str, species: list[str], ion: str) -> float:
    x = np.asarray([1.0e-8, 1.0e-8, 1.0 - 2.0e-8], dtype=float)
    params = get_prop_dict(dataset_name, species, x, T_REF, user_options={})
    rho = pcsaft_den(T_REF, P_REF, x, params, phase="liq")
    values = pcsaft_gsolv(T_REF, rho, x, params, species=species)
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
        ax.bar(
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

    ax.axhline(0.0, color="black", linewidth=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(ions)
    ax.set_ylabel(r"$\Delta G_{\mathrm{hyd},i}^{\infty}$ / kJ mol$^{-1}$")
    ax.set_title("Bulow 2020 Part I Figure 2: Gibbs Energy of Hydration with pcsaft Overlay")
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
        fontsize=7,
    )
    common.annotate_percent_deltas(
        ax,
        x + offsets[5],
        row_map["revised"],
        row_map["revised_calc"],
        fontsize=7,
    )
    common.add_percent_note(ax)

    ax.legend(ncol=3, frameon=True)
    common.save_figure(fig, OUTPUT_PATH)
    plt.close(fig)


if __name__ == "__main__":
    main()
