from __future__ import annotations

import csv
import os
import platform
import sys
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO_ROOT = Path(__file__).resolve().parents[3]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts._env import require_pcsaft_install

require_pcsaft_install()

# Avoid WMI stalls from platform.machine() during scipy import on some Windows sessions.
def _fast_machine() -> str:
    return os.environ.get("PROCESSOR_ARCHITECTURE", "AMD64")

platform.machine = _fast_machine

from pcsaft import dielc_water, pcsaft_den, pcsaft_fugcoef, pcsaft_p

T_REF = 298.15
P_REF = 1.0e5
MW_WATER = 18.0153e-3
def configure_style() -> None:
    plt.rcParams.update(
        {
            "font.size": 10,
            "font.family": "DejaVu Serif",
            "axes.linewidth": 1.0,
            "axes.grid": False,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
            "legend.frameon": False,
            "mathtext.default": "regular",
        }
    )


def save_figure(fig, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=300, bbox_inches="tight")


def water_sigma(T: float) -> float:
    return float(2.7927 + 10.11 * np.exp(-0.01775 * T) - 1.417 * np.exp(-0.01146 * T))


def water_dielc(T: float) -> float:
    # Held 2014 Section 3.1.1 uses the temperature-dependent water dielectric
    # relation from Held 2008; only the composition dependence is neglected.
    return float(dielc_water(T))


def mole_fraction_from_molality_11(molality: float) -> np.ndarray:
    n_water = 1.0 / MW_WATER
    n_cation = float(molality)
    n_anion = float(molality)
    n_total = n_water + n_cation + n_anion
    return np.asarray([n_cation, n_anion, n_water], dtype=float) / n_total


def _elec_model_no_born_constant() -> dict:
    return {
        "rel_perm": {"rule": 0, "differential_mode": "analytical"},
        "DH_model": {"d_ion_mode": 1, "bjeruum_treatment": False},
        "include_born_model": False,
        "born_model": {
            "d_Born_mode": 0,
            "solvation_shell_model": False,
            "dielectric_saturation": False,
            "bulk_mode": "mix",
            "mu_born_model": {
                "differential_mode": "analytical",
                "comp_dep_rel_perm": True,
                "include_sum_term": True,
                "comp_dep_delta_d": False,
            },
        },
    }


def _base_params(m: np.ndarray, s: np.ndarray, e: np.ndarray, k_ij: np.ndarray) -> dict:
    return {
        "m": m,
        "s": s,
        "e": e,
        "e_assoc": np.asarray([0.0, 0.0, 2425.7], dtype=float),
        "vol_a": np.asarray([0.0, 0.0, 0.0451], dtype=float),
        "z": np.asarray([1.0, -1.0, 0.0], dtype=float),
        "k_ij": np.asarray(k_ij, dtype=float),
        "dielc": np.asarray([8.0, 8.0, np.nan], dtype=float),
        "elec_model": _elec_model_no_born_constant(),
        "debug": False,
    }


def build_params(salt: str, strategy: str, T: float = T_REF) -> dict:
    sigma_w = water_sigma(T)
    dielc_w = water_dielc(T)
    strategy_key = str(strategy).strip()

    if strategy_key == "2008":
        if salt == "NaCl":
            m = np.asarray([1.0, 1.0, 1.2047], dtype=float)
            s = np.asarray([2.4122, 3.0575, sigma_w], dtype=float)
            e = np.asarray([646.05, 47.29, 353.9449], dtype=float)
        elif salt == "KBr":
            m = np.asarray([1.0, 1.0, 1.2047], dtype=float)
            s = np.asarray([2.9698, 3.4573, sigma_w], dtype=float)
            e = np.asarray([271.05, 60.22, 353.9449], dtype=float)
        elif salt == "KCl":
            m = np.asarray([1.0, 1.0, 1.2047], dtype=float)
            s = np.asarray([2.9698, 3.0575, sigma_w], dtype=float)
            e = np.asarray([271.05, 47.29, 353.9449], dtype=float)
        else:
            raise ValueError(f"Unsupported salt '{salt}' for strategy 2008.")

        k_ij = np.zeros((3, 3), dtype=float)
        k_ij[0, 1] = 1.0
        k_ij[1, 0] = 1.0
        params = _base_params(m, s, e, k_ij)
        params["dielc"][2] = dielc_w
        return params

    if strategy_key == "2014":
        if salt == "NaCl":
            m = np.asarray([1.0, 1.0, 1.2047], dtype=float)
            s = np.asarray([2.8232, 2.7560, sigma_w], dtype=float)
            e = np.asarray([230.00, 170.00, 353.9449], dtype=float)
            k_ij = np.zeros((3, 3), dtype=float)
            k_na_w = -0.007981 * T + 2.37999
            k_ij[0, 2] = k_na_w
            k_ij[2, 0] = k_na_w
            k_ij[1, 2] = -0.25
            k_ij[2, 1] = -0.25
            k_ij[0, 1] = 0.317
            k_ij[1, 0] = 0.317
            params = _base_params(m, s, e, k_ij)
            params["dielc"][2] = dielc_w
            return params

        if salt == "KBr":
            m = np.asarray([1.0, 1.0, 1.2047], dtype=float)
            s = np.asarray([3.3417, 3.0707, sigma_w], dtype=float)
            e = np.asarray([200.00, 190.00, 353.9449], dtype=float)
            k_ij = np.zeros((3, 3), dtype=float)
            k_kw = -0.004012 * T + 1.3959
            k_ij[0, 2] = k_kw
            k_ij[2, 0] = k_kw
            k_ij[1, 2] = -0.25
            k_ij[2, 1] = -0.25
            k_ij[0, 1] = -0.102
            k_ij[1, 0] = -0.102
            params = _base_params(m, s, e, k_ij)
            params["dielc"][2] = dielc_w
            return params

        if salt == "KCl":
            m = np.asarray([1.0, 1.0, 1.2047], dtype=float)
            s = np.asarray([3.3417, 2.7560, sigma_w], dtype=float)
            e = np.asarray([200.00, 170.00, 353.9449], dtype=float)
            k_ij = np.zeros((3, 3), dtype=float)
            k_kw = -0.004012 * T + 1.3959
            k_ij[0, 2] = k_kw
            k_ij[2, 0] = k_kw
            k_ij[1, 2] = -0.25
            k_ij[2, 1] = -0.25
            k_ij[0, 1] = 0.064
            k_ij[1, 0] = 0.064
            params = _base_params(m, s, e, k_ij)
            params["dielc"][2] = dielc_w
            return params

        raise ValueError(f"Unsupported salt '{salt}' for strategy 2014.")

    raise ValueError(f"Unsupported strategy '{strategy}'.")


def osmotic_molality_from_fugacity(T: float, rho: float, x: np.ndarray, params: dict) -> float:
    z = np.asarray(params["z"], dtype=float)
    idx_water = np.where(np.abs(z) <= 1e-12)[0]
    if idx_water.size != 1:
        raise ValueError(f"Expected exactly one neutral solvent, found {idx_water.size}.")
    iw = int(idx_water[0])

    molality = x / (x[iw] * MW_WATER)
    molality[iw] = 0.0

    x0 = np.zeros_like(x)
    x0[iw] = 1.0

    fugcoef = np.asarray(pcsaft_fugcoef(T, rho, x, params), dtype=float).reshape(-1)
    p_mix = pcsaft_p(T, rho, x, params)
    phase0 = "vap" if rho < 900.0 else "liq"
    rho0 = pcsaft_den(T, p_mix, x0, params, phase=phase0)
    fugcoef0 = np.asarray(pcsaft_fugcoef(T, rho0, x0, params), dtype=float).reshape(-1)
    gamma_w = fugcoef[iw] / fugcoef0[iw]

    return float(-1000.0 * np.log(x[iw] * gamma_w) / 18.0153 / np.sum(molality))


def calc_osmotic_curve(salt: str, m_values: np.ndarray, strategy: str, T: float = T_REF) -> np.ndarray:
    params = build_params(salt, strategy, T=T)
    out = np.empty_like(m_values, dtype=float)
    for i, m_salt in enumerate(np.asarray(m_values, dtype=float)):
        m_eval = max(float(m_salt), 1e-12)
        x = mole_fraction_from_molality_11(m_eval)
        rho = pcsaft_den(T, P_REF, x, params, phase="liq")
        out[i] = osmotic_molality_from_fugacity(T, rho, x, params)
    return out


def load_osmotic_data(salt: str, m_min: float = 0.0, m_max: float = 4.0) -> tuple[np.ndarray, np.ndarray]:
    candidates = [
        REPO_ROOT / "data" / "osmotic" / "water" / f"{salt}.csv",
        REPO_ROOT / "data" / "osmotic" / f"{salt}.csv",
    ]
    data_path = next((p for p in candidates if p.exists()), None)
    if data_path is None:
        raise FileNotFoundError(f"No osmotic CSV found for {salt}. Checked: {candidates}")

    molality = []
    osmotic = []
    with data_path.open("r", newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        fields = reader.fieldnames or []
        m_key = next((k for k in fields if k and "m" in k.lower()), None)
        phi_key = next((k for k in fields if k and "osmotic" in k.lower()), None)
        if m_key is None or phi_key is None:
            raise KeyError(f"Missing expected columns in {data_path}. Columns={fields}")
        for row in reader:
            try:
                m_val = float(row[m_key])
                phi_val = float(row[phi_key])
            except (TypeError, ValueError):
                continue
            if m_min <= m_val <= m_max:
                molality.append(m_val)
                osmotic.append(phi_val)

    m_arr = np.asarray(molality, dtype=float)
    phi_arr = np.asarray(osmotic, dtype=float)
    order = np.argsort(m_arr)
    return m_arr[order], phi_arr[order]


