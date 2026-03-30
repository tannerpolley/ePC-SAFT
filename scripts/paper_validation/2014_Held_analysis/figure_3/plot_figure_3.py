from __future__ import annotations

import copy
import csv
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
REPO_ROOT = Path(__file__).resolve().parents[4]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

import _common as common
from scripts._env import require_pcsaft_install

require_pcsaft_install()

from pcsaft.parameters import get_prop_dict
from pcsaft import pcsaft_den

T_REF = 298.15
P_REF = 1.0e5
OUTPUT = Path(__file__).with_name("figure_3.png")
DIGITIZED = Path(__file__).with_name("figure_3_digitized.csv")
RAW_DATA_SOURCES = (
    common.REPO_ROOT / "data" / "osmotic" / "water" / "LiAc-NaAc-KAc.csv",
    Path(__file__).with_name("LiAc-NaAc-KAc.csv"),
)

SALT_SPECS = {
    "LiAc": ("Li+", "Ac-"),
    "NaAc": ("Na+", "Ac-"),
    "KAc": ("K+", "Ac-"),
}

def _load_digitized() -> dict[str, tuple[np.ndarray, np.ndarray]]:
    raw_path = next((path for path in RAW_DATA_SOURCES if path.exists()), None)
    if raw_path is not None:
        grouped: dict[str, list[tuple[float, float]]] = {k: [] for k in SALT_SPECS}
        with raw_path.open("r", newline="", encoding="utf-8-sig") as handle:
            reader = csv.reader(handle)
            next(reader, None)
            for row in reader:
                if len(row) >= 2 and row[0] and row[1]:
                    grouped["KAc"].append((float(row[0]), float(row[1])))
                if len(row) >= 4 and row[2] and row[3]:
                    grouped["NaAc"].append((float(row[2]), float(row[3])))
                if len(row) >= 6 and row[4] and row[5]:
                    grouped["LiAc"].append((float(row[4]), float(row[5])))

        out: dict[str, tuple[np.ndarray, np.ndarray]] = {}
        for salt, pts in grouped.items():
            pts_sorted = sorted(pts, key=lambda p: p[0])
            out[salt] = (
                np.asarray([p[0] for p in pts_sorted], dtype=float),
                np.asarray([p[1] for p in pts_sorted], dtype=float),
            )
        return out

    grouped: dict[str, list[tuple[float, float]]] = {k: [] for k in SALT_SPECS}
    with DIGITIZED.open("r", newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            salt = str(row.get("salt", "")).strip()
            if salt not in grouped:
                continue
            try:
                m = float(row["molality"])
                phi = float(row["osmotic"])
            except (TypeError, ValueError, KeyError):
                continue
            grouped[salt].append((m, phi))

    out: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    for salt, pts in grouped.items():
        pts_sorted = sorted(pts, key=lambda p: p[0])
        out[salt] = (
            np.asarray([p[0] for p in pts_sorted], dtype=float),
            np.asarray([p[1] for p in pts_sorted], dtype=float),
        )
    return out


def _build_strategy2_params(salt: str) -> dict:
    cation, anion = SALT_SPECS[salt]
    species = [cation, anion, "H2O"]
    x_ref = common.mole_fraction_from_molality_11(1e-8)
    return get_prop_dict("2014_Held", species, x_ref, T_REF, user_options={})


def _build_strategy1_like_params(salt: str) -> dict:
    cation, anion = SALT_SPECS[salt]
    species = [cation, anion, "H2O"]
    x_ref = common.mole_fraction_from_molality_11(1e-8)
    params = get_prop_dict("2009_Held", species, x_ref, T_REF, user_options={})
    return {
        key: (np.asarray(val, dtype=float).copy() if isinstance(val, np.ndarray) else copy.deepcopy(val))
        for key, val in params.items()
    }


def _calc_curve(params: dict, m_grid: np.ndarray) -> np.ndarray:
    out = np.empty_like(m_grid, dtype=float)
    for i, m in enumerate(m_grid):
        m_eval = max(float(m), 1e-12)
        x = common.mole_fraction_from_molality_11(m_eval)
        rho = pcsaft_den(T_REF, P_REF, x, params, phase="liq")
        out[i] = common.osmotic_molality_from_fugacity(T_REF, rho, x, params)
    return out


def main() -> None:
    common.configure_style()
    data = _load_digitized()

    m_grid = np.linspace(0.0, 4.0, 401)
    curves_s1 = {salt: _calc_curve(_build_strategy1_like_params(salt), m_grid) for salt in SALT_SPECS}
    curves_s2 = {salt: _calc_curve(_build_strategy2_params(salt), m_grid) for salt in SALT_SPECS}

    fig, axes = plt.subplots(1, 2, figsize=(10.0, 4.8), sharey=True)

    marker_map = {"LiAc": "^", "NaAc": "o", "KAc": "*"}
    ls_map = {"LiAc": "-", "NaAc": "-", "KAc": "--"}
    color_map = {"LiAc": "0.70", "NaAc": "black", "KAc": "black"}

    for ax, curves, panel_title in (
        (axes[0], curves_s1, "(a) Classical ePC-SAFT (Held 2009)"),
        (axes[1], curves_s2, "(b) Strategy 2 (Held 2014)"),
    ):
        for salt in ("KAc", "NaAc", "LiAc"):
            m_exp, phi_exp = data[salt]
            ax.scatter(
                m_exp,
                phi_exp,
                marker=marker_map[salt],
                s=38,
                edgecolors="black",
                facecolors="0.85",
                linewidths=0.9,
                zorder=5,
                label=f"{salt} data",
            )
            ax.plot(
                m_grid,
                curves[salt],
                color=color_map[salt],
                linewidth=2.0 if salt == "NaAc" else 1.6,
                linestyle=ls_map[salt],
                zorder=3,
                label=f"{salt} model",
            )

        ax.set_xlim(0.0, 4.0)
        ax.set_ylim(0.90, 1.30)
        ax.set_xlabel(r"molality, $m_{alkali\ acetate}$ / mol kg$^{-1}$")
        ax.set_title(panel_title)
        ax.grid(True, alpha=0.24)

    axes[0].set_ylabel(r"molal osmotic coefficient, $\phi_m$")
    axes[1].legend(fontsize=7, ncol=2, loc="upper left")

    fig.suptitle("Held 2014 Fig. 3 reproduction (LiAc, NaAc, KAc)", y=1.02)
    common.save_figure(fig, OUTPUT)
    plt.close(fig)


if __name__ == "__main__":
    main()

