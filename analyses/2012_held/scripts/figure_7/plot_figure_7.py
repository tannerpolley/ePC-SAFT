from __future__ import annotations

from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import _common as common

OUTPUT = Path(__file__).with_name("figure_7.png")
DATA_PATH = Path(__file__).with_name("hernandez_2007_table6.csv")
DIELC_PATH = Path(__file__).with_name("water-methanol-ethanol.csv")
SERIES = [
    {
        "weights": {"water": 0.85, "methanol": 0.10, "ethanol": 0.05},
        "marker": "o",
        "s": 42,
        "facecolors": "white",
        "edgecolors": "#3f3f3f",
        "linecolor": "#222222",
        "label": r"$x_{\mathrm{H_2O}}^0 = 0.918,\ x_{\mathrm{MeOH}}^0 = 0.061$",
    },
    {
        "weights": {"water": 0.70, "methanol": 0.15, "ethanol": 0.15},
        "marker": "s",
        "s": 34,
        "facecolors": "#c7c7c7",
        "edgecolors": "#606060",
        "linecolor": "#505050",
        "label": r"$x_{\mathrm{H_2O}}^0 = 0.830,\ x_{\mathrm{MeOH}}^0 = 0.100$",
    },
    {
        "weights": {"water": 0.60, "methanol": 0.10, "ethanol": 0.30},
        "marker": "*",
        "s": 84,
        "facecolors": "#ededed",
        "edgecolors": "#7a7a7a",
        "linecolor": "#7a7a7a",
        "label": r"$x_{\mathrm{H_2O}}^0 = 0.776,\ x_{\mathrm{MeOH}}^0 = 0.073$",
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


def _weight_signature(weights: dict[str, float]) -> tuple[tuple[str, float], ...]:
    return tuple((solvent, round(float(weights[solvent]), 6)) for solvent in ("water", "methanol", "ethanol"))


def _load_dielc_lookup() -> dict[tuple[tuple[str, float], ...], float]:
    _, rows = common.read_csv_rows(DIELC_PATH)
    lookup: dict[tuple[tuple[str, float], ...], float] = {}
    for row in rows:
        weights = {
            "water": common.parse_float(row.get("w_H2O")),
            "methanol": common.parse_float(row.get("w_Methanol")),
            "ethanol": common.parse_float(row.get("w_Ethanol")),
        }
        dielc = common.parse_float(row.get("dielc"))
        if dielc is None or any(value is None for value in weights.values()):
            continue
        lookup[_weight_signature(weights)] = float(dielc)
    return lookup


def _load_grouped_data() -> dict[tuple[tuple[str, float], ...], list[dict[str, float]]]:
    _, rows = common.read_csv_rows(DATA_PATH)
    grouped: dict[tuple[tuple[str, float], ...], list[dict[str, float]]] = {}
    for row in rows:
        weights = {
            "water": common.parse_float(row.get("w_H2O")),
            "methanol": common.parse_float(row.get("w_Methanol")),
            "ethanol": common.parse_float(row.get("w_Ethanol")),
        }
        molality = common.parse_float(row.get("molality"))
        gamma = common.parse_float(row.get("miac_m"))
        if molality is None or gamma is None or any(value is None for value in weights.values()):
            continue
        key = _weight_signature(weights)
        grouped.setdefault(key, []).append({"molality": float(molality), "miac_m": float(gamma)})
    for values in grouped.values():
        values.sort(key=lambda item: item["molality"])
    return grouped


def main() -> None:
    common.configure_style()

    # Hernandez-Hernandez 2007 Table 6 reports salt-free solvent composition in mass fractions,
    # and the tabulated gamma values are on the molality basis. Table 4 gives
    # exact dielectric constants for the same mass-fraction solvent mixtures.
    grouped = _load_grouped_data()
    dielc_lookup = _load_dielc_lookup()

    fig, ax = plt.subplots(figsize=(7.6, 4.8))

    for series in SERIES:
        key = _weight_signature(series["weights"])
        rows = grouped[key]
        comp = _weight_to_mole_comp(series["weights"])
        expected_dielc = dielc_lookup[key]

        params_check = common.build_params(
            dataset="2012_Held",
            salt="NaCl",
            solvent_system="water-methanol-ethanol",
            comp=comp,
        )
        applied_dielc = common.parse_float(params_check.get("mixed_solvent_rel_perm"))
        if applied_dielc is None or abs(applied_dielc - expected_dielc) > 1.0e-6:
            raise ValueError(
                f"Exact mixed dielectric mismatch for {series['weights']}: "
                f"expected {expected_dielc:.6f}, got {applied_dielc!r}."
            )

        m_data = np.asarray([row["molality"] for row in rows], dtype=float)
        y_data = np.asarray([row["miac_m"] for row in rows], dtype=float)
        mask = m_data <= 1.0 + 1e-12

        grid, y_model = common.mean_ionic_activity_curve(
            dataset="2012_Held",
            salt="NaCl",
            solvent_system="water-methanol-ethanol",
            comp=comp,
            m_max=1.0,
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

    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.4, 1.0)
    ax.set_xlabel(r"molality, $m_{\mathrm{NaCl}}$ / mol kg$^{-1}$")
    ax.set_ylabel(r"$\gamma_{\pm}^{m,*}$ / -")
    ax.set_title("2012 Fig. 7 style: NaCl in water/methanol/ethanol mixtures (298.15 K)")
    ax.legend(loc="upper right", fontsize=8.8, handlelength=2.3)
    ax.grid(False)

    common.save_figure(fig, OUTPUT)
    plt.close(fig)


if __name__ == "__main__":
    main()
