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
from scripts._env import require_epcsaft_install

require_epcsaft_install()

from epcsaft.parameters import get_prop_dict
from scripts._epcsaft_oop import epcsaft_density, epcsaft_fugacity_coefficient, epcsaft_pressure

T_REF = 298.15
P_REF = 1.0e5
MW_WATER = 18.01528e-3
OUTPUT = Path(__file__).with_name("figure_4a.png")
DIGITIZED = common.analysis_data_path(__file__, "figure_4a_digitized.csv", kind="input")

ION_STRATEGY1 = {
    "K+": (2.9698, 271.05),
    "Cl-": (3.0575, 47.29),
}


def _load_digitized() -> dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]]:
    grouped: dict[str, list[tuple[float, float, float]]] = {"alanine": [], "glycine": []}
    with DIGITIZED.open("r", newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            solute = str(row.get("solute", "")).strip().lower()
            if solute not in grouped:
                continue
            try:
                m_kcl = float(row["m_kcl"])
                m_aa = float(row["m_amino_acid"])
                phi = float(row["osmotic"])
            except (TypeError, ValueError, KeyError):
                continue
            grouped[solute].append((m_kcl, m_aa, phi))

    out: dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]] = {}
    for solute, pts in grouped.items():
        pts_sorted = sorted(pts, key=lambda p: p[1])
        out[solute] = (
            np.asarray([p[0] for p in pts_sorted], dtype=float),
            np.asarray([p[1] for p in pts_sorted], dtype=float),
            np.asarray([p[2] for p in pts_sorted], dtype=float),
        )
    return out


def _mole_fraction_from_molalities(m_kcl: float, m_aa: float) -> np.ndarray:
    n_w = 1.0 / MW_WATER
    n_k = float(m_kcl)
    n_cl = float(m_kcl)
    n_aa = float(m_aa)
    n_total = n_w + n_k + n_cl + n_aa
    if n_total <= 0.0:
        raise ValueError("Non-positive total moles in molality conversion.")
    return np.asarray([n_k / n_total, n_cl / n_total, n_aa / n_total, n_w / n_total], dtype=float)


def _build_strategy2_params(solute_species: str) -> dict:
    species = ["K+", "Cl-", solute_species, "H2O"]
    x_ref = _mole_fraction_from_molalities(1.0, 1e-8)
    return get_prop_dict("2014_Held", species, x_ref, T_REF, user_options={})


def _build_strategy1_like_params(solute_species: str) -> dict:
    params = _build_strategy2_params(solute_species)
    p = {
        key: (np.asarray(val, dtype=float).copy() if isinstance(val, np.ndarray) else copy.deepcopy(val))
        for key, val in params.items()
    }

    s_arr = np.asarray(p["s"], dtype=float).copy()
    e_arr = np.asarray(p["e"], dtype=float).copy()
    s_arr[0], e_arr[0] = ION_STRATEGY1["K+"]
    s_arr[1], e_arr[1] = ION_STRATEGY1["Cl-"]
    p["s"] = s_arr
    p["e"] = e_arr

    k_ij = np.asarray(p["k_ij"], dtype=float).copy()
    idx_k, idx_cl, idx_aa, idx_w = 0, 1, 2, 3
    k_ij[idx_k, idx_w] = 0.0
    k_ij[idx_w, idx_k] = 0.0
    k_ij[idx_cl, idx_w] = 0.0
    k_ij[idx_w, idx_cl] = 0.0
    k_ij[idx_k, idx_cl] = 1.0
    k_ij[idx_cl, idx_k] = 1.0
    k_ij[idx_k, idx_aa] = 0.0
    k_ij[idx_aa, idx_k] = 0.0
    k_ij[idx_cl, idx_aa] = 0.0
    k_ij[idx_aa, idx_cl] = 0.0
    p["k_ij"] = k_ij
    return p


def _calc_osmotic_curve(params: dict, m_kcl: float, m_aa_grid: np.ndarray) -> np.ndarray:
    out = np.empty_like(m_aa_grid, dtype=float)
    for i, m_aa in enumerate(m_aa_grid):
        x = _mole_fraction_from_molalities(m_kcl, float(m_aa))
        rho = epcsaft_density(T_REF, P_REF, x, params, phase="liq")
        out[i] = _osmotic_molality_from_fugacity(T_REF, rho, x, params)
    return out


def _osmotic_molality_from_fugacity(T: float, rho: float, x: np.ndarray, params: dict) -> float:
    idx_water = 3

    molality = np.asarray(
        [
            x[0] / (x[idx_water] * MW_WATER),
            x[1] / (x[idx_water] * MW_WATER),
            x[2] / (x[idx_water] * MW_WATER),
            0.0,
        ],
        dtype=float,
    )

    x0 = np.zeros_like(x)
    x0[idx_water] = 1.0

    fugcoef = np.asarray(epcsaft_fugacity_coefficient(T, rho, x, params), dtype=float).reshape(-1)
    p_mix = epcsaft_pressure(T, rho, x, params)
    rho0 = epcsaft_density(T, p_mix, x0, params, phase="liq")
    fugcoef0 = np.asarray(epcsaft_fugacity_coefficient(T, rho0, x0, params), dtype=float).reshape(-1)
    gamma_w = fugcoef[idx_water] / fugcoef0[idx_water]

    return float(-1000.0 * np.log(x[idx_water] * gamma_w) / 18.0153 / np.sum(molality))


def main() -> None:
    common.configure_style()
    data = _load_digitized()

    m_grid = np.linspace(0.0, 3.0, 160)

    params_s2_ala = _build_strategy2_params("Alanine")
    params_s2_gly = _build_strategy2_params("Glycine")
    params_s1_ala = _build_strategy1_like_params("Alanine")
    params_s1_gly = _build_strategy1_like_params("Glycine")

    m_kcl_ala = float(np.median(data["alanine"][0]))
    m_kcl_gly = float(np.median(data["glycine"][0]))

    phi_s2_ala = _calc_osmotic_curve(params_s2_ala, m_kcl_ala, m_grid)
    phi_s2_gly = _calc_osmotic_curve(params_s2_gly, m_kcl_gly, m_grid)
    phi_s1_ala = _calc_osmotic_curve(params_s1_ala, m_kcl_ala, m_grid)
    phi_s1_gly = _calc_osmotic_curve(params_s1_gly, m_kcl_gly, m_grid)

    fig, ax = plt.subplots(figsize=(7.2, 5.0))

    _, m_ala, phi_ala = data["alanine"]
    _, m_gly, phi_gly = data["glycine"]

    ax.scatter(m_gly, phi_gly, marker="s", s=35, color="black", label="Glycine data", zorder=6)
    ax.scatter(
        m_ala,
        phi_ala,
        marker="o",
        s=54,
        facecolors="0.85",
        edgecolors="black",
        linewidths=0.9,
        label="Alanine data",
        zorder=6,
    )

    ax.plot(m_grid, phi_s2_gly, color="black", linewidth=2.2, label="Glycine model (strategy 2)")
    ax.plot(m_grid, phi_s2_ala, color="0.45", linewidth=2.2, label="Alanine model (strategy 2)")

    ax.plot(m_grid, phi_s1_gly, color="black", linewidth=1.1, linestyle="--", label="Glycine model (strategy 1-like)")
    ax.plot(m_grid, phi_s1_ala, color="0.55", linewidth=1.1, linestyle="--", label="Alanine model (strategy 1-like)")

    ax.set_xlim(0.0, 3.0)
    ax.set_ylim(0.85, 1.05)
    ax.set_xlabel(r"amino-acid molality, $m_{amino\ acid}$ / mol kg$^{-1}$")
    ax.set_ylabel(r"molal osmotic coefficient, $\phi_m$")
    ax.set_title("Held 2014 Fig. 4a reproduction (1 m KCl)")
    ax.grid(True, alpha=0.24)
    ax.legend(fontsize=8, ncol=2)

    common.save_figure(fig, OUTPUT)
    plt.close(fig)


if __name__ == "__main__":
    main()
