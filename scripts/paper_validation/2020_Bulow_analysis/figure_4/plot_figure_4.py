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

DATASETS = [
    ("a", "Water to methanol", "methanol", SCRIPT_DIR / "data" / "water-methanol-comparison.csv"),
    ("b", "Water to ethanol", "ethanol", SCRIPT_DIR / "data" / "water-ethanol-comparison.csv"),
]
SERIES = [
    ("data median", "Literature median", "#bdbdbd", None),
    ("advanced", "ePC-SAFT advanced (paper)", "#2ca02c", None),
    ("advanced_calc", "ePC-SAFT advanced (epcsaft)", "#2ca02c", "////"),
    ("revised", "ePC-SAFT revised (paper)", "#e67e22", None),
    ("revised_calc", "ePC-SAFT revised (epcsaft)", "#e67e22", "////"),
]
T_REF = 298.15
P_REF = 1.0e5
EPS = 1.0e-8
VARIANT_DATASET = {"advanced": "2020_Bulow", "revised": "2014_Held"}
SOLVENT_SPECIES = {"water": "Water", "methanol": "Methanol", "ethanol": "Ethanol"}


def _species_for_ion(ion: str, solvent: str) -> list[str]:
    solvent_species = SOLVENT_SPECIES[solvent]
    if ion in {"Li+", "Na+", "K+"}:
        return [ion, "Cl-", solvent_species]
    if ion == "F-":
        return ["Na+", "F-", solvent_species]
    return ["Na+", ion, solvent_species]


def _gsolv_ion(variant: str, ion: str, solvent: str) -> float:
    species = _species_for_ion(ion, solvent)
    dataset_name = VARIANT_DATASET[variant]
    x = np.asarray([EPS, EPS, 1.0 - 2.0 * EPS], dtype=float)
    params = get_prop_dict(dataset_name, species, x, T_REF, user_options={})
    rho = epcsaft_density(T_REF, P_REF, x, params, phase="liq")
    values = epcsaft_solvation_free_energy(T_REF, rho, x, params, species=species)
    return float(values[ion]) / 1000.0


def _transfer_total(variant: str, ion: str, solvent: str) -> float:
    return _gsolv_ion(variant, ion, solvent) - _gsolv_ion(variant, ion, "water")


def _computed_values(ions: list[str], solvent: str, variant: str) -> np.ndarray:
    out = np.empty(len(ions), dtype=float)
    for idx, ion in enumerate(ions):
        out[idx] = _transfer_total(variant, ion, solvent)
    return out


def _plot_one(panel_tag: str, title: str, solvent: str, data_path: Path) -> None:
    print(f"Computing Figure 4{panel_tag} ({solvent})", flush=True)
    frame = common.load_indexed_csv(data_path)
    ions = list(frame.columns)
    row_map = {
        "data low": frame.values("data low", ions),
        "data median": frame.values("data median", ions),
        "data high": frame.values("data high", ions),
        "advanced": frame.values("advanced", ions),
        "advanced_calc": _computed_values(ions, solvent, "advanced"),
        "revised": frame.values("revised", ions),
        "revised_calc": _computed_values(ions, solvent, "revised"),
    }

    x = np.arange(len(ions), dtype=float)
    width = 0.14
    offsets = np.linspace(-2.0 * width, 2.0 * width, len(SERIES))

    fig, ax = plt.subplots(figsize=(11.8, 5.9))
    for offset, (row_key, label, color, hatch) in zip(offsets, SERIES):
        values = row_map[row_key]
        bars = ax.bar(
            x + offset,
            values,
            width=width,
            color=color,
            edgecolor="black",
            linewidth=0.45,
            hatch=hatch,
            alpha=0.9 if hatch is None else 0.75,
            label=label,
        )
        common.annotate_bar_values(ax, bars, fontsize=6)

    ax.errorbar(
        x + offsets[0],
        row_map["data median"],
        yerr=np.vstack(
            [
                row_map["data median"] - row_map["data low"],
                row_map["data high"] - row_map["data median"],
            ]
        ),
        fmt="none",
        ecolor="0.35",
        elinewidth=1.1,
        capsize=4.0,
        zorder=5,
    )

    values = np.vstack([row_map[key] for key in row_map])
    y_min = float(np.nanmin(values))
    y_max = float(np.nanmax(values))

    ax.axhline(0.0, color="black", linewidth=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(ions)
    ax.set_ylabel(r"$\Delta G_{\mathrm{tr},i}^{\infty}$ / kJ mol$^{-1}$")
    ax.set_title(f"Bulow 2020 Part I Figure 4{panel_tag}: {title}")
    ax.set_ylim(y_min - 5.0, y_max + 8.0)
    ax.grid(axis="y", alpha=0.22)
    common.annotate_percent_deltas(
        ax,
        x + offsets[2],
        row_map["advanced"],
        row_map["advanced_calc"],
        xs_ref=x + offsets[1],
    )
    common.annotate_percent_deltas(
        ax,
        x + offsets[4],
        row_map["revised"],
        row_map["revised_calc"],
        xs_ref=x + offsets[3],
    )
    common.add_percent_note(ax)
    ax.legend(ncol=3, frameon=True)

    output_path = SCRIPT_DIR / f"figure_4{panel_tag}.png"
    common.save_figure(fig, output_path)
    plt.close(fig)
    print(f"Wrote {output_path}", flush=True)


def main() -> None:
    common.configure_style()
    for panel_tag, title, solvent, data_path in DATASETS:
        _plot_one(panel_tag, title, solvent, data_path)


if __name__ == "__main__":
    main()
