"""Dataset-driven MIAC fit validation (v2).

This script validates MIAC_m datasets using parameter sets from
`data/pcsaft_parameters/<dataset>/` through `data.epcsaft_properties`.

It writes canonical fit plots to:
  data/MIAC_m/<solvent>/plot_fits/maic_m_<solvent>_<Salt>_fit.png
"""

from __future__ import annotations

import csv
import math
import os
import sys
from collections import OrderedDict
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import matplotlib
import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from data.epcsaft_properties import get_prop_dict, molality_to_molefraction
from pcsaft import pcsaft_den, pcsaft_miac_m

matplotlib.use("Agg")
import matplotlib.pyplot as plt


T_REF = 298.15
P_REF = 1.0e5
AXIS_LABEL_SIZE = 13
AXIS_TICK_SIZE = 11

SALT_SPECS = {
    "LiBr": {"cation": "Li+", "anion": "Br-"},
    "LiCl": {"cation": "Li+", "anion": "Cl-"},
    "LiI": {"cation": "Li+", "anion": "I-"},
    "NaBr": {"cation": "Na+", "anion": "Br-"},
    "NaCl": {"cation": "Na+", "anion": "Cl-"},
    "NaI": {"cation": "Na+", "anion": "I-"},
    "KCl": {"cation": "K+", "anion": "Cl-"},
    "KBr": {"cation": "K+", "anion": "Br-"},
    "KI": {"cation": "K+", "anion": "I-"},
}

DATASET_VARIANTS: "OrderedDict[str, Dict[str, object]]" = OrderedDict(
    [
        ("cameretti_2005", {"label": "2005", "color": "black", "linestyle": "--", "lw": 1.6, "user_options": {}}),
        ("held_2008", {"label": "2008", "color": "dimgray", "linestyle": "--", "lw": 1.6, "user_options": {}}),
        ("held_2014", {"label": "2014", "color": "silver", "linestyle": "--", "lw": 1.6, "user_options": {}}),
        ("bulow_2020", {"label": "2020", "color": "tab:purple", "linestyle": "--", "lw": 1.8, "user_options": {}}),
        ("figiel_2025", {"label": "2025", "color": "tab:red", "linestyle": "--", "lw": 1.8, "user_options": {}}),
    ]
)


def _requested_scope() -> Tuple[str, str]:
    solvent = os.getenv("MIAC_SOLVENT", "all").strip().lower()
    if solvent not in {"all", "water", "methanol", "ethanol"}:
        raise ValueError("MIAC_SOLVENT must be one of: all, water, methanol, ethanol.")
    salt = os.getenv("MIAC_SALT", "").strip()
    return solvent, salt


def _canonical_salt_token(salt: str) -> str:
    for item in SALT_SPECS:
        if item.lower() == salt.lower():
            return item
    return salt


def discover_combos(solvent_scope: str | None = None, salt_scope: str | None = None) -> List[Dict[str, object]]:
    selected_solvent, selected_salt = _requested_scope()
    if solvent_scope is not None:
        selected_solvent = solvent_scope.strip().lower()
    if salt_scope is not None:
        selected_salt = salt_scope.strip()

    combos: List[Dict[str, object]] = []
    for solvent in ("Water", "Methanol", "Ethanol"):
        solvent_lower = solvent.lower()
        if selected_solvent != "all" and selected_solvent != solvent_lower:
            continue

        data_dir = REPO_ROOT / "data" / "MIAC_m" / solvent_lower
        if not data_dir.exists():
            continue

        if solvent == "Water":
            candidates = sorted(data_dir.glob("*.csv"))
        else:
            prefix = f"{solvent_lower}-"
            candidates = sorted(data_dir.glob(f"{prefix}*.csv"))

        for path in candidates:
            if solvent == "Water":
                salt = path.stem
            else:
                salt = path.stem.replace(prefix, "", 1)
            if salt not in SALT_SPECS:
                continue
            if selected_salt and salt.lower() != selected_salt.lower():
                continue
            salt_token = _canonical_salt_token(salt)
            output_dir = path.parent / "plot_fits"
            combos.append(
                {
                    "salt": salt,
                    "solvent": solvent,
                    "data_path": path,
                    "output": output_dir / f"maic_m_{solvent_lower}_{salt_token}_fit.png",
                }
            )

    if not combos:
        scope = f"solvent={selected_solvent}, salt={selected_salt or '<all>'}"
        raise FileNotFoundError(f"No MIAC CSV files found for {scope}.")
    return combos


def _species_for_combo(salt: str, solvent: str) -> List[str]:
    salt_spec = SALT_SPECS[salt]
    if solvent == "Water":
        water_species = "H2O-2B-Li" if salt.startswith("Li") else "H2O-2B-NaCl"
        return [salt_spec["cation"], salt_spec["anion"], water_species]
    return [salt_spec["cation"], salt_spec["anion"], solvent]


def _pair_key(salt: str) -> str:
    salt_spec = SALT_SPECS[salt]
    return f"{salt_spec['cation']}{salt_spec['anion']}"


def load_exp_data(combo: Dict[str, object]) -> Tuple[np.ndarray, np.ndarray]:
    path = Path(combo["data_path"])
    salt = str(combo["salt"])
    molal: List[float] = []
    miac: List[float] = []

    with path.open("r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        fields = set(reader.fieldnames or [])

        m_key = next((candidate for candidate in ("molality (kg/mol)", f"m_{salt}", "molality", "m") if candidate in fields), None)
        if m_key is None:
            raise ValueError(
                f"Missing molality column in {path}. Tried: 'molality (kg/mol)', 'm_{salt}', 'molality', 'm'."
            )

        y_key = next((candidate for candidate in ("miac_m", "gamma") if candidate in fields), None)
        if y_key is None:
            raise ValueError(f"Missing MIAC column in {path}. Tried: 'miac_m', 'gamma'.")

        for row in reader:
            molal.append(float(row[m_key]))
            miac.append(float(row[y_key]))

    molal_arr = np.asarray(molal, dtype=float)
    miac_arr = np.asarray(miac, dtype=float)
    if molal_arr.size == 0:
        raise ValueError(f"No data rows in {path}.")
    if not np.all(np.isfinite(molal_arr)) or not np.all(np.isfinite(miac_arr)):
        raise ValueError(f"Non-finite experimental values in {path}.")

    order = np.argsort(molal_arr)
    return molal_arr[order], miac_arr[order]


def build_params_for_variant(dataset_name: str, species: List[str], user_options: dict | None = None) -> Dict[str, object]:
    x_ref = molality_to_molefraction(1e-8, species=species)
    return get_prop_dict(dataset_name, species, x_ref, T_REF, user_options=user_options)


def calc_curve(combo: Dict[str, object], dataset_name: str, molal_grid: np.ndarray) -> np.ndarray:
    salt = str(combo["salt"])
    solvent = str(combo["solvent"])
    species = _species_for_combo(salt, solvent)
    user_options = dict(DATASET_VARIANTS[dataset_name].get("user_options", {}))
    params = build_params_for_variant(dataset_name, species, user_options=user_options)
    pair_key = _pair_key(salt)

    gamma = np.empty_like(molal_grid, dtype=float)
    for idx, m in enumerate(molal_grid):
        m_eval = float(m) if m > 0.0 else 1e-12
        x = molality_to_molefraction(m_eval, species=species)
        rho = pcsaft_den(T_REF, P_REF, x, params, phase="liq")
        gamma[idx] = pcsaft_miac_m(T_REF, rho, x, params, species=species)[pair_key]

    if not np.all(np.isfinite(gamma)):
        raise ValueError(f"Non-finite MIAC values for {salt}/{solvent} in dataset {dataset_name}.")
    return gamma


def plot_combo(combo: Dict[str, object], output_path: Path | None = None, save: bool = True, close: bool = True, ax=None):
    salt = str(combo["salt"])
    solvent = str(combo["solvent"])
    molal_exp, miac_exp = load_exp_data(combo)

    max_molality = float(np.max(molal_exp))
    if max_molality < 1.0:
        xmax = round(max(0.1, math.ceil(max_molality * 10.0) / 10.0), 1)
    else:
        xmax = float(min(10, int(math.ceil(max_molality))))
    visible_mask = molal_exp <= (xmax + 1e-12)
    if not np.any(visible_mask):
        visible_mask = np.ones_like(molal_exp, dtype=bool)
    ymax = float(max(1, int(math.ceil(float(np.max(miac_exp[visible_mask]))))))

    molal_grid = np.linspace(0.0, xmax, 701)

    curves = {}
    for dataset_name in DATASET_VARIANTS:
        curves[dataset_name] = calc_curve(combo, dataset_name, molal_grid)

    created_fig = False
    if ax is None:
        fig, ax = plt.subplots(figsize=(8.0, 5.4))
        created_fig = True
    else:
        fig = ax.figure

    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")

    ax.scatter(
        molal_exp,
        miac_exp,
        color="black",
        marker="o",
        s=34,
        facecolors="none",
        label=f"{salt} data ({solvent})",
    )

    for dataset_name, cfg in DATASET_VARIANTS.items():
        ax.plot(
            molal_grid,
            curves[dataset_name],
            color=str(cfg["color"]),
            linestyle=str(cfg["linestyle"]),
            linewidth=float(cfg["lw"]),
            label=f"{salt} {cfg['label']}",
        )

    ax.set_xlim(0.0, xmax)
    ax.set_ylim(0.0, ymax)
    ax.set_xlabel(r"molality, $m$ / mol kg$^{-1}$", fontsize=AXIS_LABEL_SIZE)
    ax.set_ylabel(r"mean ionic activity coefficient, $\gamma_{\pm}^{m}$", fontsize=AXIS_LABEL_SIZE)
    ax.xaxis.label.set_color("black")
    ax.yaxis.label.set_color("black")
    ax.title.set_color("black")
    ax.tick_params(axis="both", labelsize=AXIS_TICK_SIZE)
    ax.tick_params(colors="black")
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color("black")
        spine.set_linewidth(1.0)
    ax.set_title(f"{salt} in {solvent.lower()} at 298.15 K")
    ax.grid(True, alpha=0.3, color="0.7")
    legend = ax.legend(fontsize=8)
    frame = legend.get_frame()
    frame.set_facecolor("white")
    frame.set_edgecolor("black")
    frame.set_alpha(1.0)
    for text in legend.get_texts():
        text.set_color("black")

    if save:
        if output_path is None:
            output_path = Path(combo["output"])
        output_path.parent.mkdir(parents=True, exist_ok=True)
        fig.tight_layout()
        fig.savefig(output_path, dpi=220)
        if not output_path.exists():
            raise FileNotFoundError(f"Expected plot was not written: {output_path}")

    if created_fig and close:
        plt.close(fig)

    return {
        "figure": fig,
        "axis": ax,
        "output_path": output_path,
        "molality_exp": molal_exp,
        "miac_exp": miac_exp,
        "molality_grid": molal_grid,
        "curves": curves,
    }


def run_validate_miac_fits_v2() -> List[Path]:
    combos = discover_combos()
    generated: List[Path] = []

    for combo in combos:
        result = plot_combo(combo, save=True, close=True)
        generated.append(Path(result["output_path"]))

    print("Dataset variants:")
    for dataset_name, cfg in DATASET_VARIANTS.items():
        print(f"- {dataset_name} -> {cfg['label']}")
    print("Generated validation plots:")
    for path in generated:
        print(f"- {path}")
    return generated


def test_validate_miac_fits_v2() -> None:
    run_validate_miac_fits_v2()


if __name__ == "__main__":
    run_validate_miac_fits_v2()