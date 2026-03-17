from __future__ import annotations

import csv
import math
import sys
from pathlib import Path

import matplotlib
import numpy as np

ANALYSIS_ROOT = Path(__file__).resolve().parents[1]
REPO_ROOT = Path(__file__).resolve().parents[4]
if str(ANALYSIS_ROOT) not in sys.path:
    sys.path.insert(0, str(ANALYSIS_ROOT))
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from data.epcsaft_properties import get_prop_dict
from pcsaft import pcsaft_den, pcsaft_miac
from _plot_common import configure_style, save_figure

matplotlib.use("Agg")
import matplotlib.pyplot as plt

T_REF = 298.15
P_REF = 1.0e5
DATA_DIR = Path(__file__).resolve().parent / "data"

SOLVENT_MW = {
    "methanol": 32.04e-3,
    "ethanol": 46.068e-3,
}

SALT_SPECS = {
    "LiBr": {"cation": "Li+", "anion": "Br-"},
    "LiCl": {"cation": "Li+", "anion": "Cl-"},
    "NaBr": {"cation": "Na+", "anion": "Br-"},
}

PANELS = [
    {"id": "7a", "salt": "LiBr", "solvent": "ethanol", "data": DATA_DIR / "ethanol-LiBr.csv", "xlim": (0.0, 0.2), "ylim": (0.0, 4.0)},
    {"id": "7b", "salt": "NaBr", "solvent": "ethanol", "data": DATA_DIR / "ethanol-NaBr.csv", "xlim": (0.0, 0.01), "ylim": (0.0, 1.0)},
    {"id": "7c", "salt": "LiCl", "solvent": "methanol", "data": DATA_DIR / "methanol-LiCl.csv", "xlim": (0.0, 0.15), "ylim": (0.0, 2.0)},
    {"id": "7d", "salt": "NaBr", "solvent": "methanol", "data": DATA_DIR / "methanol-NaBr.csv", "xlim": (0.0, 0.05), "ylim": (0.0, 1.0)},
]


def _species_for_combo(salt: str, solvent: str) -> list[str]:
    spec = SALT_SPECS[salt]
    solvent_name = "Methanol" if solvent == "methanol" else "Ethanol"
    return [spec["cation"], spec["anion"], solvent_name]


def _resolve_pair_key(result: dict[str, float], salt: str) -> str:
    spec = SALT_SPECS[salt]
    target = f"{spec['cation']}{spec['anion']}"
    if target in result:
        return target
    for key in result:
        if spec["cation"] in key and spec["anion"] in key:
            return key
    raise KeyError(f"Could not resolve MIAC key for {salt}. Keys={list(result.keys())}")


def _molality_for_salt_mole_fraction(x_salt: np.ndarray | float, solvent: str) -> np.ndarray:
    x = np.asarray(x_salt, dtype=float)
    n_solv = 1.0 / SOLVENT_MW[solvent]
    return np.where(x <= 0.0, 0.0, (x * n_solv) / np.maximum(1.0 - x, 1.0e-12))


def _molality_to_species_molefraction(molality: float, salt: str, solvent: str) -> np.ndarray:
    m = max(float(molality), 0.0)
    n_solv = 1.0 / SOLVENT_MW[solvent]
    n_cat = m
    n_an = m
    n_total = n_solv + n_cat + n_an
    return np.asarray([n_cat / n_total, n_an / n_total, n_solv / n_total], dtype=float)


def _build_params(dataset: str, salt: str, solvent: str) -> dict:
    x_ref = _molality_to_species_molefraction(1e-8, salt, solvent)
    return get_prop_dict(dataset, _species_for_combo(salt, solvent), x_ref, T_REF)


def _miac_curve(dataset: str, salt: str, solvent: str, x_grid: np.ndarray) -> np.ndarray:
    params = _build_params(dataset, salt, solvent)
    species = _species_for_combo(salt, solvent)
    pair_key = None
    y = np.empty_like(x_grid, dtype=float)
    m_grid = _molality_for_salt_mole_fraction(x_grid, solvent)
    for idx, m in enumerate(m_grid):
        x = _molality_to_species_molefraction(float(max(m, 1.0e-12)), salt, solvent)
        rho = pcsaft_den(T_REF, P_REF, x, params, phase="liq")
        vals = pcsaft_miac(T_REF, rho, x, params, species=species)
        if pair_key is None:
            pair_key = _resolve_pair_key(vals, salt)
        y[idx] = float(vals[pair_key])
    return y


def _read_dataset(path: Path, solvent: str) -> tuple[np.ndarray, np.ndarray]:
    xs: list[float] = []
    ys: list[float] = []
    with path.open("r", newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            try:
                x_val = float(str(row.get("mole_fraction", "")).strip())
            except ValueError:
                x_val = math.nan
            try:
                y_val = float(str(row.get("miac", "")).strip())
            except ValueError:
                y_val = math.nan
            if not math.isfinite(y_val):
                try:
                    m_val = float(str(row.get("molality", "")).strip())
                    y_m = float(str(row.get("miac_m", row.get("gamma", ""))).strip())
                    y_val = y_m * (1.0 + SOLVENT_MW[solvent] * m_val * 2.0)
                except ValueError:
                    y_val = math.nan
            if not math.isfinite(x_val):
                try:
                    m_val = float(str(row.get("molality", "")).strip())
                    x_val = float(_molality_for_salt_mole_fraction(np.asarray([m_val], dtype=float), solvent)[0])
                except ValueError:
                    x_val = math.nan
            if math.isfinite(x_val) and math.isfinite(y_val):
                xs.append(float(x_val))
                ys.append(float(y_val))
    x_arr = np.asarray(xs, dtype=float)
    y_arr = np.asarray(ys, dtype=float)
    order = np.argsort(x_arr)
    return x_arr[order], y_arr[order]


def _plot_panel(panel: dict[str, object]) -> Path:
    salt = str(panel["salt"])
    solvent = str(panel["solvent"])
    x_data, y_data = _read_dataset(Path(panel["data"]), solvent)
    x_min, x_max = panel["xlim"]
    y_min, y_max = panel["ylim"]

    x_grid = np.linspace(float(x_min), float(x_max), 1201)
    y_rev = _miac_curve("held_2014", salt, solvent, x_grid)
    y_adv = _miac_curve("bulow_2020", salt, solvent, x_grid)

    fig, ax = plt.subplots(figsize=(5.2, 4.0))
    ax.scatter(x_data, y_data, s=24, facecolor="white", edgecolor="black", linewidth=0.8, label="Literature", zorder=5)
    ax.plot(x_grid, y_rev, color="black", linewidth=1.8, label="ePC-SAFT revised", zorder=3)
    ax.plot(x_grid, y_adv, color="#228b22", linewidth=1.9, label="ePC-SAFT advanced", zorder=4)
    ax.set_xlim(float(x_min), float(x_max))
    ax.set_ylim(float(y_min), float(y_max))
    ax.set_xlabel(r"salt mole fraction, $x_{salt}$")
    ax.set_ylabel(r"$\gamma_{\pm}^{*}$ / -")
    ax.set_title(f"Bulow 2020 Figure {panel['id']}: {salt} in {solvent} at 298.15 K and 1 bar")
    ax.grid(True, alpha=0.25)
    ax.legend(loc="best")
    out = Path(__file__).resolve().parent / f"figure_{panel['id']}.png"
    save_figure(fig, out)
    plt.close(fig)
    print(f"Wrote {out}")
    return out


def main() -> None:
    configure_style()
    for panel in PANELS:
        _plot_panel(panel)


if __name__ == "__main__":
    main()
