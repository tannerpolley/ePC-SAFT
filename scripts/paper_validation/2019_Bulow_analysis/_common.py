from __future__ import annotations

import os
import platform
import sys
from pathlib import Path

import numpy as np

# Avoid WMI stalls from platform.machine() during scipy / extension imports.
platform.machine = lambda: os.environ.get("PROCESSOR_ARCHITECTURE", "AMD64")

ROOT = Path(__file__).resolve().parent
REPO_ROOT = ROOT.parents[2]

if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts._env import require_epcsaft_install

require_epcsaft_install()

from epcsaft.parameters import get_prop_dict


P_REF = 1.0e5
DATASET = "2019_Bulow"
WATER_DIELC = 78.09
IL_DIELC = 11.0

IMIDAZOLIUM_NTF2 = ["C2mim+", "C4mim+", "C6mim+", "C8mim+"]
ANION_SERIES = ["BF4-", "TFO-", "NTf2-", "PF6-"]
FIG7_ILS = [
    ("C2mim+", "NTf2-"),
    ("C3mim+", "NTf2-"),
    ("C4mim+", "NTf2-"),
    ("C5mim+", "NTf2-"),
    ("C6mim+", "NTf2-"),
    ("C7mim+", "NTf2-"),
    ("C8mim+", "NTf2-"),
    ("C3mpyr+", "NTf2-"),
    ("C4mpyr+", "NTf2-"),
    ("C6mpyr+", "NTf2-"),
    ("C4m4py+", "NTf2-"),
    ("C6mim+", "BF4-"),
    ("C8mim+", "BF4-"),
    ("C8mim+", "TFO-"),
    ("C4py+", "NTf2-"),
    ("C6py+", "NTf2-"),
    ("C8py+", "NTf2-"),
]


def il_label(cation: str, anion: str) -> str:
    return f"[{cation.replace('+', '')}][{anion.replace('-', '')}]"


def build_params(
    species: list[str],
    t: float,
    *,
    use_kij: bool = True,
    model_mode: str = "epc",
) -> dict:
    z_feed = np.asarray([0.5, 0.25, 0.25], dtype=float)
    params = get_prop_dict(DATASET, species, z_feed, t, user_options={})
    params["debug"] = False

    if not use_kij:
        k_ij = np.asarray(params["k_ij"], dtype=float).copy()
        k_ij[0, 1:] = 0.0
        k_ij[1:, 0] = 0.0
        params["k_ij"] = k_ij

    if model_mode == "epc":
        return params

    params = dict(params)
    dielc = np.asarray(params["dielc"], dtype=float).copy()
    dielc[:] = WATER_DIELC if model_mode == "orig_water" else IL_DIELC
    params["dielc"] = dielc

    elec = dict(params["elec_model"])
    elec["include_born_model"] = False
    rel_perm = dict(elec["rel_perm"])
    rel_perm["rule"] = 0
    rel_perm["differential_mode"] = 0
    elec["rel_perm"] = rel_perm
    dh_model = dict(elec["DH_model"])
    mu_dh = dict(dh_model["mu_DH_model"])
    mu_dh["comp_dep_rel_perm"] = False
    dh_model["mu_DH_model"] = mu_dh
    elec["DH_model"] = dh_model
    params["elec_model"] = elec
    return params


def solve_binary_lle(
    cation: str,
    anion: str,
    t: float,
    *,
    use_kij: bool = True,
    model_mode: str = "epc",
    z_feed: np.ndarray | None = None,
) -> dict:
    raise NotImplementedError("The legacy multiphase LLE workflow has been removed and will be rewritten later.")


def split_branches(out: dict) -> tuple[np.ndarray, np.ndarray]:
    phases = out.get("phases", [])
    if len(phases) != 2:
        raise RuntimeError(f"Expected 2 phases, got {len(phases)}")
    x0 = np.asarray(phases[0]["x"], dtype=float)
    x1 = np.asarray(phases[1]["x"], dtype=float)
    if x0[0] <= x1[0]:
        il_rich, water_rich = x0, x1
    else:
        il_rich, water_rich = x1, x0
    return il_rich, water_rich


def il_mole_fraction(x: np.ndarray) -> float:
    return float(x[1] + x[2])


def water_mole_fraction(x: np.ndarray) -> float:
    return float(x[0])


def scan_temperature_branch(
    cation: str,
    anion: str,
    temperatures: np.ndarray,
    *,
    use_kij: bool = True,
    model_mode: str = "epc",
) -> dict[str, np.ndarray]:
    t_keep: list[float] = []
    il_rich_xil: list[float] = []
    water_rich_xil: list[float] = []

    for t in np.asarray(temperatures, dtype=float):
        try:
            out = solve_binary_lle(cation, anion, float(t), use_kij=use_kij, model_mode=model_mode)
            il_rich, water_rich = split_branches(out)
        except Exception:
            continue
        t_keep.append(float(t))
        il_rich_xil.append(il_mole_fraction(il_rich))
        water_rich_xil.append(il_mole_fraction(water_rich))

    return {
        "T": np.asarray(t_keep, dtype=float),
        "x_il_il_rich": np.asarray(il_rich_xil, dtype=float),
        "x_il_water_rich": np.asarray(water_rich_xil, dtype=float),
    }


def water_solubility_in_il(cation: str, anion: str, *, t: float = 298.15, use_kij: bool = True, model_mode: str = "epc") -> float:
    try:
        out = solve_binary_lle(cation, anion, t, use_kij=use_kij, model_mode=model_mode)
        il_rich, _ = split_branches(out)
    except Exception:
        return float("nan")
    return water_mole_fraction(il_rich)


