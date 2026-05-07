from __future__ import annotations

import csv
import math
import os
import platform
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO_ROOT = Path(__file__).resolve().parents[3]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts._env import require_epcsaft_install
from scripts.plot_outputs import analysis_data_path, paper_validation_output_path, save_plot_figure

require_epcsaft_install()


# Avoid WMI stalls from platform.machine() during native extension imports on some Windows sessions.
def _fast_machine() -> str:
    return os.environ.get("PROCESSOR_ARCHITECTURE", "AMD64")


platform.machine = _fast_machine

from epcsaft.parameters import get_prop_dict
from scripts._epcsaft_oop import epcsaft_activity_coefficient, epcsaft_density, epcsaft_fugacity_coefficient

T_REF = 298.15
P_REF = 1.0e5

SOLVENT_MW = {
    "water": 18.01528e-3,
    "methanol": 32.04e-3,
    "ethanol": 46.068e-3,
}

SALT_SPECS = {
    "LiCl": {"cation": "Li+", "anion": "Cl-", "z_cation": 1, "z_anion": -1},
    "LiBr": {"cation": "Li+", "anion": "Br-", "z_cation": 1, "z_anion": -1},
    "LiI": {"cation": "Li+", "anion": "I-", "z_cation": 1, "z_anion": -1},
    "NaCl": {"cation": "Na+", "anion": "Cl-", "z_cation": 1, "z_anion": -1},
    "NaBr": {"cation": "Na+", "anion": "Br-", "z_cation": 1, "z_anion": -1},
    "NaI": {"cation": "Na+", "anion": "I-", "z_cation": 1, "z_anion": -1},
    "KCl": {"cation": "K+", "anion": "Cl-", "z_cation": 1, "z_anion": -1},
    "KBr": {"cation": "K+", "anion": "Br-", "z_cation": 1, "z_anion": -1},
    "KI": {"cation": "K+", "anion": "I-", "z_cation": 1, "z_anion": -1},
    "NH4Cl": {"cation": "NH4+", "anion": "Cl-", "z_cation": 1, "z_anion": -1},
    "NH4Br": {"cation": "NH4+", "anion": "Br-", "z_cation": 1, "z_anion": -1},
    "NH4I": {"cation": "NH4+", "anion": "I-", "z_cation": 1, "z_anion": -1},
}


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
    path = paper_validation_output_path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    try:
        save_plot_figure(fig, path, dpi=300)
    except MemoryError:
        # Fallback for intermittent Agg allocator spikes on high-DPI tight bboxes.
        save_plot_figure(fig, path, dpi=180, bbox_inches=None)
    except RuntimeError as exc:
        if "bad allocation" not in str(exc).lower():
            raise
        # Fallback for intermittent Agg allocator spikes on high-DPI tight bboxes.
        save_plot_figure(fig, path, dpi=180, bbox_inches=None)


def parse_float(value) -> float | None:
    if value is None:
        return None
    text = str(value).strip()
    if not text:
        return None
    try:
        val = float(text)
    except ValueError:
        return None
    if not math.isfinite(val):
        return None
    return float(val)


def read_csv_rows(path: Path) -> Tuple[List[str], List[Dict[str, str]]]:
    with path.open("r", newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        fields = [h.strip() for h in (reader.fieldnames or []) if h and h.strip()]
        rows: List[Dict[str, str]] = []
        for row in reader:
            clean: Dict[str, str] = {}
            for key, value in row.items():
                if not key:
                    continue
                ks = key.strip()
                if not ks:
                    continue
                clean[ks] = value.strip() if isinstance(value, str) else value
            rows.append(clean)
    return fields, rows


def stoich_for_salt(salt: str) -> Tuple[int, int]:
    spec = SALT_SPECS[salt]
    zc = abs(int(round(spec["z_cation"])))
    za = abs(int(round(spec["z_anion"])))
    g = math.gcd(zc, za)
    return za // g, zc // g


def species_for_combo(salt: str, solvent_system: str) -> List[str]:
    spec = SALT_SPECS[salt]
    solvents = [s for s in solvent_system.split("-") if s]
    solvent_species = []
    for s in solvents:
        if s == "water":
            solvent_species.append("H2O")
        elif s == "methanol":
            solvent_species.append("Methanol")
        elif s == "ethanol":
            solvent_species.append("Ethanol")
        else:
            raise ValueError(f"Unsupported solvent '{s}'.")
    return [spec["cation"], spec["anion"], *solvent_species]


def normalized_comp(solvent_system: str, comp: Dict[str, float]) -> Dict[str, float]:
    solvents = [s for s in solvent_system.split("-") if s]
    if len(solvents) == 1:
        return {solvents[0]: 1.0}
    vals = {s: float(comp.get(s, 0.0)) for s in solvents}
    total = sum(vals.values())
    if total <= 0.0:
        return {s: 1.0 / len(solvents) for s in solvents}
    return {s: vals[s] / total for s in solvents}


def _weight_to_mole_fractions(weights: Dict[str, float]) -> Dict[str, float]:
    moles = {}
    for s, w in weights.items():
        mw = SOLVENT_MW[s]
        moles[s] = w / mw if mw > 0.0 else 0.0
    total = sum(moles.values())
    if total <= 0.0:
        return {s: 0.0 for s in weights}
    return {s: moles[s] / total for s in weights}


def extract_comp(row: Dict[str, str], solvent_system: str) -> Dict[str, float]:
    solvents = [s for s in solvent_system.split("-") if s]
    if len(solvents) == 1:
        return {solvents[0]: 1.0}

    xmap = {
        "x_h2o": "water",
        "x_water": "water",
        "x_methanol": "methanol",
        "x_meoh": "methanol",
        "x_ethanol": "ethanol",
        "x_etoh": "ethanol",
    }
    wmap = {
        "w_h2o": "water",
        "w_water": "water",
        "w_h2o_salt_free": "water",
        "w_water_salt_free": "water",
        "w_methanol": "methanol",
        "w_meoh": "methanol",
        "w_methanol_salt_free": "methanol",
        "w_meoh_salt_free": "methanol",
        "w_ethanol": "ethanol",
        "w_etoh": "ethanol",
        "w_ethanol_salt_free": "ethanol",
        "w_etoh_salt_free": "ethanol",
    }

    xs: Dict[str, float] = {}
    ws: Dict[str, float] = {}
    for key, value in row.items():
        lk = key.strip().lower()
        if lk in xmap and xmap[lk] in solvents:
            val = parse_float(value)
            if val is not None:
                xs[xmap[lk]] = val
        if lk in wmap and wmap[lk] in solvents:
            val = parse_float(value)
            if val is not None:
                ws[wmap[lk]] = val

    if xs:
        if len(solvents) == 2 and len(xs) == 1:
            known = next(iter(xs))
            other = [s for s in solvents if s != known][0]
            xs[other] = 1.0 - xs[known]
        return normalized_comp(solvent_system, xs)

    if ws:
        if len(solvents) == 2 and len(ws) == 1:
            known = next(iter(ws))
            other = [s for s in solvents if s != known][0]
            ws[other] = 1.0 - ws[known]
        mole_frac = _weight_to_mole_fractions(ws)
        return normalized_comp(solvent_system, mole_frac)

    return normalized_comp(solvent_system, {})


def comp_signature(solvent_system: str, comp: Dict[str, float]) -> Tuple[Tuple[str, float], ...]:
    solvents = [s for s in solvent_system.split("-") if s]
    if len(solvents) <= 1:
        return tuple()
    return tuple((s, round(float(comp.get(s, 0.0)), 6)) for s in solvents)


def signature_label(signature: Tuple[Tuple[str, float], ...]) -> str:
    if not signature:
        return "pure solvent"
    key_map = {"water": "H2O", "methanol": "MeOH", "ethanol": "EtOH"}
    chunks = [f"x_{key_map.get(sol, sol)}^0={val:.2f}" for sol, val in signature]
    return ", ".join(chunks)


def mw_mix(solvent_system: str, comp: Dict[str, float]) -> float:
    frac = normalized_comp(solvent_system, comp)
    return sum(frac[s] * SOLVENT_MW[s] for s in frac)


def molality_to_species_molefraction(
    molality: float, salt: str, solvent_system: str, comp: Dict[str, float]
) -> np.ndarray:
    species = species_for_combo(salt, solvent_system)
    solvents = [s for s in solvent_system.split("-") if s]
    solvent_species = species[2:]
    frac = normalized_comp(solvent_system, comp)

    n_solv_total = 1.0 / mw_mix(solvent_system, frac)
    nu_cat, nu_an = stoich_for_salt(salt)

    n_totals: Dict[str, float] = {sp: 0.0 for sp in species}
    for s_key, sp in zip(solvents, solvent_species):
        n_totals[sp] += frac[s_key] * n_solv_total
    n_totals[species[0]] += nu_cat * float(molality)
    n_totals[species[1]] += nu_an * float(molality)

    total = sum(n_totals.values())
    if total <= 0.0:
        raise ValueError("Non-positive total moles while converting molality to mole fractions.")

    return np.asarray([n_totals[sp] / total for sp in species], dtype=float)


def pair_key_for_salt(salt: str) -> str:
    spec = SALT_SPECS[salt]
    return f"{spec['cation']}{spec['anion']}"


def _resolve_pair_key(result: Dict[str, float], salt: str) -> str:
    target = pair_key_for_salt(salt)
    if target in result:
        return target
    cat = SALT_SPECS[salt]["cation"]
    an = SALT_SPECS[salt]["anion"]
    for key in result.keys():
        if cat in key and an in key:
            return key
    raise KeyError(f"Could not resolve mean-ionic key for salt {salt}. Keys={list(result.keys())}")


def build_params(
    dataset: str, salt: str, solvent_system: str, comp: Dict[str, float], user_options: dict | None = None
) -> dict:
    x_ref = molality_to_species_molefraction(1e-8, salt, solvent_system, comp)
    return get_prop_dict(
        dataset, species_for_combo(salt, solvent_system), x_ref, T_REF, user_options=user_options or {}
    )


def mean_ionic_activity_curve(
    dataset: str,
    salt: str,
    solvent_system: str,
    comp: Dict[str, float],
    m_max: float,
    points: int = 600,
    user_options: dict | None = None,
) -> Tuple[np.ndarray, np.ndarray]:
    grid = np.linspace(0.0, float(m_max), int(points))
    params = build_params(dataset, salt, solvent_system, comp, user_options=user_options)
    species = species_for_combo(salt, solvent_system)

    gamma = np.empty_like(grid)
    for idx, m in enumerate(grid):
        m_eval = max(float(m), 1e-12)
        x = molality_to_species_molefraction(m_eval, salt, solvent_system, comp)
        rho = epcsaft_density(T_REF, P_REF, x, params, phase="liq")
        vals = epcsaft_activity_coefficient(
            T_REF, rho, x, params, species=species, mean_ionic_form=True, basis="molality"
        )
        key = _resolve_pair_key(vals, salt)
        gamma[idx] = float(vals[key])

    return grid, gamma


def mean_ionic_activity_curve_x(
    dataset: str,
    salt: str,
    solvent_system: str,
    comp: Dict[str, float],
    m_max: float,
    points: int = 600,
    user_options: dict | None = None,
) -> Tuple[np.ndarray, np.ndarray]:
    grid = np.linspace(0.0, float(m_max), int(points))
    params = build_params(dataset, salt, solvent_system, comp, user_options=user_options)
    species = species_for_combo(salt, solvent_system)

    gamma = np.empty_like(grid)
    for idx, m in enumerate(grid):
        m_eval = max(float(m), 1e-12)
        x = molality_to_species_molefraction(m_eval, salt, solvent_system, comp)
        rho = epcsaft_density(T_REF, P_REF, x, params, phase="liq")
        vals = epcsaft_activity_coefficient(T_REF, rho, x, params, species=species, mean_ionic_form=True, basis="mole")
        key = _resolve_pair_key(vals, salt)
        gamma[idx] = float(vals[key])

    return grid, gamma


def solvent_activity_curve(
    dataset: str,
    salt: str,
    solvent_system: str,
    m_max: float,
    points: int = 600,
    user_options: dict | None = None,
) -> Tuple[np.ndarray, np.ndarray]:
    solvents = [s for s in solvent_system.split("-") if s]
    if len(solvents) != 1:
        raise ValueError("solvent_activity_curve currently supports only single-solvent systems.")

    comp = {solvents[0]: 1.0}
    species = species_for_combo(salt, solvent_system)
    params = build_params(dataset, salt, solvent_system, comp, user_options=user_options)
    solvent_idx = 2

    x_pure = np.zeros(len(species), dtype=float)
    x_pure[solvent_idx] = 1.0
    rho_pure = epcsaft_density(T_REF, P_REF, x_pure, params, phase="liq")
    fugcoef_pure = np.asarray(epcsaft_fugacity_coefficient(T_REF, rho_pure, x_pure, params), dtype=float)
    solvent_ref = float(fugcoef_pure[solvent_idx])
    if solvent_ref <= 0.0 or not math.isfinite(solvent_ref):
        raise ValueError("Pure-solvent fugacity coefficient is not finite and positive.")

    grid = np.linspace(0.0, float(m_max), int(points))
    gamma = np.empty_like(grid)
    for idx, m in enumerate(grid):
        m_eval = max(float(m), 1e-12)
        x = molality_to_species_molefraction(m_eval, salt, solvent_system, comp)
        rho = epcsaft_density(T_REF, P_REF, x, params, phase="liq")
        fugcoef = np.asarray(epcsaft_fugacity_coefficient(T_REF, rho, x, params), dtype=float)
        gamma[idx] = float(fugcoef[solvent_idx] / solvent_ref)

    return grid, gamma


def read_miac_dataset(path: Path, solvent_system: str) -> List[Dict[str, object]]:
    fields, rows = read_csv_rows(path)
    lookup = {field.lower(): field for field in fields}

    molality_key = next((lookup[c] for c in ("molality", "m") if c in lookup), None)
    miac_m_key = next((lookup[c] for c in ("miac_m", "gamma") if c in lookup), None)
    miac_key = next((lookup[c] for c in ("miac", "y") if c in lookup), None)

    if molality_key is None or miac_m_key is None:
        raise ValueError(f"Missing required columns in {path}.")

    data: List[Dict[str, object]] = []
    for row in rows:
        m = parse_float(row.get(molality_key))
        gm = parse_float(row.get(miac_m_key))
        if m is None or gm is None:
            continue
        comp = extract_comp(row, solvent_system)
        data.append(
            {
                "molality": m,
                "miac_m": gm,
                "miac": parse_float(row.get(miac_key)) if miac_key else None,
                "comp": comp,
                "signature": comp_signature(solvent_system, comp),
            }
        )

    data.sort(key=lambda x: float(x["molality"]))
    return data


def group_by_signature(
    entries: List[Dict[str, object]],
) -> Dict[Tuple[Tuple[str, float], ...], List[Dict[str, object]]]:
    grouped: Dict[Tuple[Tuple[str, float], ...], List[Dict[str, object]]] = defaultdict(list)
    for entry in entries:
        grouped[entry["signature"]].append(entry)
    for group in grouped.values():
        group.sort(key=lambda x: float(x["molality"]))
    return dict(grouped)
