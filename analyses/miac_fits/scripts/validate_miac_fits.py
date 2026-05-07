"""Dataset-driven MIAC fit validation.

This script validates MIAC datasets using parameter sets from the packaged
`epcsaft.parameters` datasets.

Experimental data source is canonical `data/reference/MIAC/**` with `miac` and `miac_m` values.

It writes fit plots to:
  analyses/miac_fits/results/miac/<solvent_system>/miac_m/<plot_set>/<plot_set>.png
  analyses/miac_fits/results/miac/<solvent_system>/miac/<plot_set>/<plot_set>.png
"""

from __future__ import annotations

import csv
import math
import os
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from collections import OrderedDict
from pathlib import Path
from typing import Dict, List, Literal, Tuple

import matplotlib
import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[3]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts._env import require_epcsaft_install

require_epcsaft_install()

from epcsaft.parameters import get_prop_dict
from scripts._epcsaft_oop import as_mixture
from scripts.plot_outputs import fits_plot_path, save_plot_figure

matplotlib.use("Agg")
import matplotlib.pyplot as plt

T_REF = 298.15
P_REF = 1.0e5
AXIS_LABEL_SIZE = 13
AXIS_TICK_SIZE = 11
SOURCE_COLOR_CYCLE = [
    "tab:blue",
    "tab:orange",
    "tab:green",
    "tab:red",
    "tab:purple",
    "tab:brown",
    "tab:pink",
    "tab:olive",
    "tab:cyan",
]

SALT_SPECS = {
    "LiBr": {"cation": "Li+", "anion": "Br-", "z_cation": 1, "z_anion": -1},
    "LiCl": {"cation": "Li+", "anion": "Cl-", "z_cation": 1, "z_anion": -1},
    "LiI": {"cation": "Li+", "anion": "I-", "z_cation": 1, "z_anion": -1},
    "NaBr": {"cation": "Na+", "anion": "Br-", "z_cation": 1, "z_anion": -1},
    "NaCl": {"cation": "Na+", "anion": "Cl-", "z_cation": 1, "z_anion": -1},
    "NaI": {"cation": "Na+", "anion": "I-", "z_cation": 1, "z_anion": -1},
    "KCl": {"cation": "K+", "anion": "Cl-", "z_cation": 1, "z_anion": -1},
    "KBr": {"cation": "K+", "anion": "Br-", "z_cation": 1, "z_anion": -1},
    "KI": {"cation": "K+", "anion": "I-", "z_cation": 1, "z_anion": -1},
}

CATION_ORDER = {"Li+": 1, "Na+": 2, "K+": 3}
ANION_ORDER = {"Cl-": 1, "Br-": 2, "I-": 3}
GRID_SALTS = [
    ["LiCl", "LiBr", "LiI"],
    ["NaCl", "NaBr", "NaI"],
    ["KCl", "KBr", "KI"],
]

SOLVENT_MW = {
    "water": 18.01528e-3,
    "methanol": 32.04e-3,
    "ethanol": 46.068e-3,
}

SOLVENT_SHORT = {
    "water": "h2o",
    "methanol": "meoh",
    "ethanol": "etoh",
}

SALT_CANONICAL = {k.lower(): k for k in SALT_SPECS}
SALT_CANONICAL["li"] = "LiI"


DATASET_VARIANTS: "OrderedDict[str, Dict[str, object]]" = OrderedDict(
    [
        ("2005_Cameretti", {"label": "2005", "color": "black", "linestyle": "--", "lw": 1.6, "user_options": {}}),
        ("2008_Held", {"label": "2008", "color": "dimgray", "linestyle": "--", "lw": 1.6, "user_options": {}}),
        ("2012_Held", {"label": "2012", "color": "tab:blue", "linestyle": "--", "lw": 1.7, "user_options": {}}),
        ("2014_Held", {"label": "2014", "color": "silver", "linestyle": "--", "lw": 1.6, "user_options": {}}),
        ("2020_Bulow", {"label": "2020", "color": "tab:purple", "linestyle": "--", "lw": 1.8, "user_options": {}}),
        ("2025_Figiel", {"label": "2025", "color": "tab:red", "linestyle": "--", "lw": 1.8, "user_options": {}}),
    ]
)

AQUEOUS_ONLY_DATASETS = {"2005_Cameretti", "2008_Held"}


def _variant_names_for_solvent_system(solvent_system: str) -> List[str]:
    if solvent_system == "water":
        # Water plots now compare modern sets only.
        return ["2012_Held", "2014_Held", "2020_Bulow", "2025_Figiel"]
    return [name for name in DATASET_VARIANTS if name not in AQUEOUS_ONLY_DATASETS]


def _available_solvent_systems() -> List[str]:
    data_root = REPO_ROOT / "data" / "reference" / "MIAC"
    if not data_root.exists():
        return []
    return sorted(p.name.lower().replace("_", "-") for p in data_root.iterdir() if p.is_dir())


def _requested_scope() -> Tuple[str, str]:
    solvent = os.getenv("MIAC_SOLVENT", "all").strip().lower().replace("_", "-")
    available = set(_available_solvent_systems())
    if solvent != "all" and solvent not in available:
        raise ValueError(f"MIAC_SOLVENT must be one of: all, {', '.join(sorted(available))}.")
    salt = os.getenv("MIAC_SALT", "").strip()
    return solvent, salt


def _requested_quantities() -> List[Literal["miac_m", "miac"]]:
    token = os.getenv("MIAC_QUANTITY", "all").strip().lower()
    if token in {"all", "both", ""}:
        return ["miac_m", "miac"]
    if token in {"miac_m", "miac"}:
        return [token]
    raise ValueError("MIAC_QUANTITY must be one of: all, miac_m, miac.")


def _requested_grid_points() -> int:
    token = os.getenv("MIAC_GRID_POINTS", "701").strip()
    try:
        value = int(token)
    except ValueError as exc:
        raise ValueError("MIAC_GRID_POINTS must be an integer.") from exc
    if value < 11:
        raise ValueError("MIAC_GRID_POINTS must be at least 11.")
    return value


def _requested_workers(total_jobs: int) -> int:
    token = os.getenv("MIAC_WORKERS", "").strip()
    if not token:
        return max(1, min(total_jobs, os.cpu_count() or 1))
    try:
        value = int(token)
    except ValueError as exc:
        raise ValueError("MIAC_WORKERS must be an integer.") from exc
    if value < 1:
        raise ValueError("MIAC_WORKERS must be at least 1.")
    return min(total_jobs, value)


def _canonical_salt_token(salt: str) -> str:
    return SALT_CANONICAL.get(salt.lower(), salt)


def _parse_salt_from_stem(stem: str) -> str:
    token = stem.split("-")[-1] if "-" in stem else stem
    return _canonical_salt_token(token)


def _solvent_key_from_col(col: str) -> str | None:
    c = col.strip().lower()
    if c in {"x_h2o", "x_water"}:
        return "water"
    if c in {"x_methanol", "x_meoh"}:
        return "methanol"
    if c in {"x_ethanol", "x_etoh"}:
        return "ethanol"
    if c.startswith("x_"):
        tail = c[2:]
        if tail in SOLVENT_MW:
            return tail
    return None


def _read_csv_rows(path: Path) -> Tuple[List[str], List[Dict[str, str]]]:
    with path.open("r", newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        fields = [h.strip() for h in (reader.fieldnames or []) if h and h.strip()]
        rows = []
        for row in reader:
            clean = {}
            for k, v in row.items():
                if not k:
                    continue
                ks = k.strip()
                if not ks:
                    continue
                clean[ks] = v.strip() if isinstance(v, str) else v
            rows.append(clean)
    return fields, rows


def _extract_comp(row: Dict[str, str], solvent_system: str) -> Dict[str, float]:
    solvents = [s for s in solvent_system.split("-") if s]
    if len(solvents) <= 1:
        return {solvents[0]: 1.0} if solvents else {}

    comp: Dict[str, float] = {}
    for col, val in row.items():
        key = _solvent_key_from_col(col)
        if key is None:
            continue
        if key not in solvents:
            continue
        try:
            f = float(val)
        except (TypeError, ValueError):
            continue
        if math.isfinite(f):
            comp[key] = f

    if len(solvents) == 2 and len(comp) == 1:
        s_known = next(iter(comp.keys()))
        s_other = [s for s in solvents if s != s_known][0]
        comp[s_other] = 1.0 - comp[s_known]

    if not comp:
        return {s: 1.0 / len(solvents) for s in solvents}

    for s in solvents:
        comp.setdefault(s, 0.0)
    denom = sum(comp.values())
    if abs(denom) < 1e-12:
        return {s: 1.0 / len(solvents) for s in solvents}
    return {s: comp[s] / denom for s in solvents}


def _comp_signature(comp: Dict[str, float], solvent_system: str) -> Tuple[Tuple[str, float], ...]:
    solvents = [s for s in solvent_system.split("-") if s]
    if len(solvents) <= 1:
        return tuple()
    return tuple((s, round(float(comp.get(s, 0.0)), 6)) for s in solvents)


def _comp_suffix(comp: Dict[str, float], solvent_system: str) -> str:
    sig = _comp_signature(comp, solvent_system)
    if not sig:
        return ""
    parts = []
    for s, x in sig:
        short = SOLVENT_SHORT.get(s, s[:4])
        parts.append(f"{short}{int(round(100 * x)):02d}")
    return "_".join(parts)


def _comp_label(comp: Dict[str, float], solvent_system: str) -> str:
    sig = _comp_signature(comp, solvent_system)
    if not sig:
        return ""
    pieces = []
    for s, x in sig:
        if s == "water":
            name = "H2O"
        else:
            name = s.capitalize()
        pieces.append(f"x_{name}={x:.2f}")
    return ", ".join(pieces)


def _salt_rank(salt: str) -> int:
    spec = SALT_SPECS[salt]
    c_rank = CATION_ORDER.get(spec["cation"])
    a_rank = ANION_ORDER.get(spec["anion"])
    if c_rank is None or a_rank is None:
        raise KeyError(f"Unsupported salt for ranking: {salt}")
    return (c_rank - 1) * 3 + a_rank


def _stoich_for_salt(salt: str) -> Tuple[int, int]:
    spec = SALT_SPECS[salt]
    zc = int(round(abs(spec["z_cation"])))
    za = int(round(abs(spec["z_anion"])))
    g = math.gcd(zc, za)
    nu_cat = za // g
    nu_an = zc // g
    return nu_cat, nu_an


def _sum_nu_for_salt(salt: str) -> int:
    nu_cat, nu_an = _stoich_for_salt(salt)
    return int(nu_cat + nu_an)


def _normalized_comp(solvent_system: str, comp: Dict[str, float]) -> Dict[str, float]:
    solvents = [s for s in solvent_system.split("-") if s]
    if not solvents:
        raise ValueError("Empty solvent system while normalizing composition.")
    if len(solvents) == 1:
        return {solvents[0]: 1.0}

    frac = {s: float(comp.get(s, 0.0)) for s in solvents}
    denom = sum(frac.values())
    if abs(denom) < 1e-12:
        return {s: 1.0 / len(solvents) for s in solvents}
    return {s: frac[s] / denom for s in solvents}


def _mw_mix(solvent_system: str, comp: Dict[str, float]) -> float:
    frac = _normalized_comp(solvent_system, comp)
    mw_mix = 0.0
    for solvent, x_solvent in frac.items():
        if solvent not in SOLVENT_MW:
            raise ValueError(f"Missing solvent MW for '{solvent}'.")
        mw_mix += x_solvent * SOLVENT_MW[solvent]
    if mw_mix <= 0.0:
        raise ValueError("Invalid mixed-solvent molar mass while computing conversion factor.")
    return float(mw_mix)


def _miac_conversion_factor(molality: np.ndarray, salt: str, solvent_system: str, comp: Dict[str, float]) -> np.ndarray:
    mw_mix = _mw_mix(solvent_system, comp)
    sum_nu = _sum_nu_for_salt(salt)
    m = np.asarray(molality, dtype=float)
    return 1.0 + mw_mix * m * float(sum_nu)


def _salt_mole_fraction_from_molality(molality: np.ndarray, solvent_system: str, comp: Dict[str, float]) -> np.ndarray:
    mw_mix = _mw_mix(solvent_system, comp)
    m = np.asarray(molality, dtype=float)
    n_solv = 1.0 / mw_mix
    denom = m + n_solv
    x_salt = np.where(denom > 0.0, m / denom, 0.0)
    return x_salt.astype(float)


def _molality_axis_upper(max_molality: float) -> float:
    if max_molality < 1.0:
        return round(max(0.1, math.ceil(max_molality * 10.0) / 10.0), 1)
    return float(min(10, int(math.ceil(max_molality))))


def _mole_fraction_axis_upper(max_x: float) -> float:
    if max_x <= 0.0:
        return 0.01
    if max_x < 0.02:
        return round(math.ceil(max_x * 1000.0) / 1000.0, 3)
    if max_x < 0.2:
        return round(math.ceil(max_x * 100.0) / 100.0, 2)
    return round(math.ceil(max_x * 10.0) / 10.0, 1)


def _source_label(raw: str) -> str:
    s = str(raw or "").strip()
    return s if s else "Unspecified source"


def _color_for_source(source_label: str) -> str:
    if source_label == "Unspecified source":
        return "black"
    idx = sum(ord(ch) for ch in source_label) % len(SOURCE_COLOR_CYCLE)
    return SOURCE_COLOR_CYCLE[idx]


def _high_outlier_mask(values: np.ndarray) -> np.ndarray:
    """Return mask that removes only a single extreme high spike, if clearly separated."""
    y = np.asarray(values, dtype=float)
    finite = np.isfinite(y)
    positive = y > 0.0
    good_idx = np.where(finite & positive)[0]
    keep = finite.copy()

    # Conservative rule: only drop the single max point when it is both absolutely
    # large and strongly separated from the second-highest value.
    if good_idx.size >= 6:
        vals = y[good_idx]
        order = np.argsort(vals)
        v_max = float(vals[order[-1]])
        v_second = float(vals[order[-2]]) if vals.size >= 2 else 0.0
        if v_max > 5.0 and v_max > 3.0 * max(v_second, 1e-12):
            drop_idx = int(good_idx[order[-1]])
            keep[drop_idx] = False

    return keep


def discover_combos(solvent_scope: str | None = None, salt_scope: str | None = None) -> List[Dict[str, object]]:
    selected_solvent, selected_salt = _requested_scope()
    if solvent_scope is not None:
        selected_solvent = solvent_scope.strip().lower().replace("_", "-")
    if salt_scope is not None:
        selected_salt = salt_scope.strip()

    combos: List[Dict[str, object]] = []
    data_root = REPO_ROOT / "data" / "reference" / "MIAC"

    for data_dir in sorted([p for p in data_root.iterdir() if p.is_dir()]):
        solvent_system = data_dir.name.lower().replace("_", "-")
        if selected_solvent != "all" and selected_solvent != solvent_system:
            continue

        for path in sorted(data_dir.glob("*.csv")):
            salt = _parse_salt_from_stem(path.stem)
            if salt not in SALT_SPECS:
                continue
            if selected_salt and salt.lower() != selected_salt.lower():
                continue

            fields, rows = _read_csv_rows(path)
            comp_cols = [c for c in fields if _solvent_key_from_col(c) is not None]

            comp_groups: Dict[Tuple[Tuple[str, float], ...], Dict[str, float]] = {}
            if comp_cols:
                for row in rows:
                    comp = _extract_comp(row, solvent_system)
                    comp_groups[_comp_signature(comp, solvent_system)] = comp
            else:
                comp = _extract_comp({}, solvent_system)
                comp_groups[_comp_signature(comp, solvent_system)] = comp

            output_dir_miac_m = fits_plot_path("miac", solvent_system, "miac_m", "_placeholder").parent
            output_dir_miac = fits_plot_path("miac", solvent_system, "miac", "_placeholder").parent
            for sig, comp in sorted(comp_groups.items()):
                suffix = _comp_suffix(comp, solvent_system)
                stem = f"{solvent_system}_{_salt_rank(salt)}_{salt}"
                if suffix:
                    stem += f"_{suffix}"
                combos.append(
                    {
                        "salt": salt,
                        "solvent_system": solvent_system,
                        "data_path": path,
                        "comp_signature": sig,
                        "comp": comp,
                        "comp_label": _comp_label(comp, solvent_system),
                        "output_miac_m": output_dir_miac_m / f"miac_m_{stem}.png",
                        "output_miac": output_dir_miac / f"miac_{stem}.png",
                    }
                )

    if not combos:
        scope = f"solvent={selected_solvent}, salt={selected_salt or '<all>'}"
        raise FileNotFoundError(f"No MIAC CSV files found for {scope}.")
    return combos


def _species_for_combo(salt: str, solvent_system: str) -> List[str]:
    salt_spec = SALT_SPECS[salt]
    solvents = [s for s in solvent_system.split("-") if s]

    solvent_species: List[str] = []
    for s in solvents:
        if s == "water":
            solvent_species.append("H2O-2B-Li" if salt.startswith("Li") else "H2O-2B-NaCl")
        elif s == "methanol":
            solvent_species.append("Methanol")
        elif s == "ethanol":
            solvent_species.append("Ethanol")
        else:
            raise ValueError(f"Unsupported solvent key '{s}' for combo species mapping.")

    return [salt_spec["cation"], salt_spec["anion"], *solvent_species]


def _pair_key(salt: str) -> str:
    salt_spec = SALT_SPECS[salt]
    return f"{salt_spec['cation']}{salt_spec['anion']}"


def _molality_to_molefraction_combo(
    molality: float, salt: str, solvent_system: str, comp: Dict[str, float]
) -> np.ndarray:
    species = _species_for_combo(salt, solvent_system)
    solvents = [s for s in solvent_system.split("-") if s]
    solvent_species = species[2:]

    if len(solvents) != len(solvent_species):
        raise ValueError("Solvent mapping mismatch while converting molality to mole fractions.")

    nu_cat, nu_an = _stoich_for_salt(salt)

    frac = _normalized_comp(solvent_system, comp)

    m_mix = _mw_mix(solvent_system, frac)
    n_solv_total = 1.0 / m_mix

    n_totals: Dict[str, float] = {sp: 0.0 for sp in species}
    for s_key, sp_name in zip(solvents, solvent_species):
        n_totals[sp_name] += frac[s_key] * n_solv_total

    n_totals[species[0]] += nu_cat * molality
    n_totals[species[1]] += nu_an * molality

    n_total = sum(n_totals.values())
    if n_total <= 0.0:
        raise ValueError("Computed non-positive total moles.")

    return np.asarray([n_totals[sp] / n_total for sp in species], dtype=float)


def load_exp_data(
    combo: Dict[str, object],
    quantity: Literal["miac_m", "miac"] = "miac_m",
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    path = Path(combo["data_path"])
    target_sig = tuple(combo.get("comp_signature", tuple()))
    solvent_system = str(combo["solvent_system"])
    comp_ref = dict(combo.get("comp", {}))

    fields, rows = _read_csv_rows(path)

    lookup = {f.lower(): f for f in fields}
    m_key = None
    for candidate in ("molality", "molality (kg/mol)", "m", "m (mol/kg)", "x"):
        if candidate in lookup:
            m_key = lookup[candidate]
            break
    if m_key is None:
        raise ValueError(f"Missing molality column in {path}.")

    y_key = None
    if quantity == "miac_m":
        candidates = ("miac_m", "gamma")
    else:
        candidates = ("miac", "y")
    for candidate in candidates:
        if candidate in lookup:
            y_key = lookup[candidate]
            break
    if y_key is None:
        raise ValueError(f"Missing {quantity}-compatible column in {path}. Tried: {candidates}.")

    x_key = lookup.get("mole_fraction") if quantity == "miac" else None
    source_key = lookup.get("source")

    molal: List[float] = []
    values: List[float] = []
    x_plot: List[float] = []
    src_labels: List[str] = []

    for row in rows:
        comp = _extract_comp(row, solvent_system)
        sig = _comp_signature(comp, solvent_system)
        if target_sig and sig != target_sig:
            continue

        m_val = row.get(m_key)
        y_val = row.get(y_key)
        try:
            m = float(m_val)
            y = float(y_val)
        except (TypeError, ValueError):
            continue
        if not (math.isfinite(m) and math.isfinite(y)):
            continue

        if quantity == "miac":
            x_val = None
            if x_key is not None:
                try:
                    x_candidate = float(row.get(x_key, ""))
                    if math.isfinite(x_candidate):
                        x_val = x_candidate
                except (TypeError, ValueError):
                    x_val = None
            if x_val is None:
                x_val = float(
                    _salt_mole_fraction_from_molality(np.asarray([m], dtype=float), solvent_system, comp_ref)[0]
                )
            x_plot.append(x_val)
        else:
            x_plot.append(m)

        src_raw = row.get(source_key, "") if source_key else ""
        src_labels.append(_source_label(src_raw))
        molal.append(m)
        values.append(y)

    molal_arr = np.asarray(molal, dtype=float)
    values_arr = np.asarray(values, dtype=float)
    x_plot_arr = np.asarray(x_plot, dtype=float)
    src_arr = np.asarray(src_labels, dtype=object)
    if molal_arr.size == 0:
        raise ValueError(f"No usable data rows in {path} for composition {combo.get('comp_label', '')}.")

    order = np.argsort(molal_arr)
    return molal_arr[order], values_arr[order], x_plot_arr[order], src_arr[order]


def _filter_exp_series(
    molal_exp_raw: np.ndarray,
    values_exp_raw: np.ndarray,
    x_exp_raw: np.ndarray,
    source_exp_raw: np.ndarray,
) -> Dict[str, object]:
    keep = _high_outlier_mask(values_exp_raw)
    removed = int(np.count_nonzero(~keep))
    if np.any(keep):
        return {
            "molality_exp": molal_exp_raw[keep],
            "values_exp": values_exp_raw[keep],
            "x_exp": x_exp_raw[keep],
            "source_exp": source_exp_raw[keep],
            "removed": removed,
        }
    return {
        "molality_exp": molal_exp_raw,
        "values_exp": values_exp_raw,
        "x_exp": x_exp_raw,
        "source_exp": source_exp_raw,
        "removed": 0,
    }


def build_params_for_variant(
    dataset_name: str, combo: Dict[str, object], user_options: dict | None = None
) -> Dict[str, object]:
    x_ref = _molality_to_molefraction_combo(
        1e-8, str(combo["salt"]), str(combo["solvent_system"]), dict(combo.get("comp", {}))
    )
    species = _species_for_combo(str(combo["salt"]), str(combo["solvent_system"]))
    return get_prop_dict(dataset_name, species, x_ref, T_REF, user_options=user_options)


def calc_gamma_m_curve(
    combo: Dict[str, object],
    dataset_name: str,
    molal_grid: np.ndarray,
) -> np.ndarray:
    salt = str(combo["salt"])
    solvent_system = str(combo["solvent_system"])
    comp = dict(combo.get("comp", {}))
    species = _species_for_combo(salt, solvent_system)
    user_options = dict(DATASET_VARIANTS[dataset_name].get("user_options", {}))
    params = build_params_for_variant(dataset_name, combo, user_options=user_options)
    mixture = as_mixture(params, species=species)
    pair_key = _pair_key(salt)

    gamma_m = np.empty_like(molal_grid, dtype=float)
    for idx, m in enumerate(molal_grid):
        m_eval = float(m) if m > 0.0 else 1e-12
        x = _molality_to_molefraction_combo(m_eval, salt, solvent_system, comp)
        state = mixture.state(T=T_REF, x=x, P=P_REF, phase="liq")
        gamma_m[idx] = float(
            state.activity_coefficient(species=species, mean_ionic_form=True, basis="molality")[pair_key]
        )

    if not np.all(np.isfinite(gamma_m)):
        raise ValueError(f"Non-finite MIAC_m values for {salt}/{solvent_system} in dataset {dataset_name}.")
    return gamma_m


def calc_curve(
    combo: Dict[str, object],
    dataset_name: str,
    molal_grid: np.ndarray,
    quantity: Literal["miac_m", "miac"] = "miac_m",
) -> np.ndarray:
    gamma_m = calc_gamma_m_curve(combo, dataset_name, molal_grid)
    if quantity == "miac_m":
        return gamma_m

    salt = str(combo["salt"])
    solvent_system = str(combo["solvent_system"])
    comp = dict(combo.get("comp", {}))
    factor = _miac_conversion_factor(molal_grid, salt, solvent_system, comp)
    gamma = gamma_m * factor
    if not np.all(np.isfinite(gamma)):
        raise ValueError(f"Non-finite MIAC values for {salt}/{solvent_system} in dataset {dataset_name}.")
    return gamma


def prepare_combo_payload(combo: Dict[str, object], grid_points: int = 701) -> Dict[str, object]:
    salt = str(combo["salt"])
    solvent_system = str(combo["solvent_system"])
    comp = dict(combo.get("comp", {}))

    exp_miac_m = _filter_exp_series(*load_exp_data(combo, quantity="miac_m"))
    exp_miac = _filter_exp_series(*load_exp_data(combo, quantity="miac"))

    max_molality = max(
        float(np.max(np.asarray(exp_miac_m["molality_exp"], dtype=float))),
        float(np.max(np.asarray(exp_miac["molality_exp"], dtype=float))),
    )
    m_upper = _molality_axis_upper(max_molality)
    molal_grid = np.linspace(0.0, m_upper, grid_points)
    x_grid_miac = _salt_mole_fraction_from_molality(molal_grid, solvent_system, comp)

    curves_gamma_m: Dict[str, np.ndarray] = {}
    curves_miac: Dict[str, np.ndarray] = {}
    active_variants: List[str] = []
    for dataset_name in _variant_names_for_solvent_system(solvent_system):
        try:
            gamma_m = calc_gamma_m_curve(combo, dataset_name, molal_grid)
        except Exception:
            continue
        curves_gamma_m[dataset_name] = gamma_m
        curves_miac[dataset_name] = gamma_m * _miac_conversion_factor(molal_grid, salt, solvent_system, comp)
        active_variants.append(dataset_name)

    return {
        "combo": combo,
        "molality_grid": molal_grid,
        "x_grids": {
            "miac_m": molal_grid,
            "miac": x_grid_miac,
        },
        "exp": {
            "miac_m": exp_miac_m,
            "miac": exp_miac,
        },
        "curves": {
            "miac_m": curves_gamma_m,
            "miac": curves_miac,
        },
        "active_variants": active_variants,
    }


def plot_combo(
    combo: Dict[str, object],
    output_path: Path | None = None,
    save: bool = True,
    close: bool = True,
    ax=None,
    show_legend: bool = True,
    quantity: Literal["miac_m", "miac"] = "miac_m",
    payload: Dict[str, object] | None = None,
):
    salt = str(combo["salt"])
    solvent_system = str(combo["solvent_system"])

    if payload is None:
        payload = prepare_combo_payload(combo)
    exp_payload = dict(payload["exp"][quantity])
    molal_exp = np.asarray(exp_payload["molality_exp"], dtype=float)
    values_exp = np.asarray(exp_payload["values_exp"], dtype=float)
    x_exp = np.asarray(exp_payload["x_exp"], dtype=float)
    source_exp = np.asarray(exp_payload["source_exp"], dtype=object)
    removed = int(exp_payload["removed"])

    if removed > 0:
        print(
            f"[outlier-filter] {salt}/{solvent_system} [{quantity}] removed {removed} high outlier experimental point(s)."
        )

    molal_grid = np.asarray(payload["molality_grid"], dtype=float)
    if quantity == "miac_m":
        x_grid = molal_grid
        x_upper = float(molal_grid[-1])
        x_label = r"molality, $m$ / mol kg$^{-1}$"
    else:
        x_grid = np.asarray(payload["x_grids"]["miac"], dtype=float)
        x_upper = _mole_fraction_axis_upper(float(np.max(x_exp)))
        x_label = r"salt mole fraction, $x_{salt}$"

    visible_mask = x_exp <= (x_upper + 1e-12)
    if not np.any(visible_mask):
        visible_mask = np.ones_like(x_exp, dtype=bool)
    ymax = float(max(1, int(math.ceil(float(np.max(values_exp[visible_mask]))))))

    curves = dict(payload["curves"][quantity])
    active_variants = list(payload["active_variants"])

    created_fig = False
    if ax is None:
        fig, ax = plt.subplots(figsize=(8.0, 5.4))
        created_fig = True
    else:
        fig = ax.figure

    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")

    unique_sources = []
    for src in source_exp.tolist():
        if src not in unique_sources:
            unique_sources.append(src)

    if len(unique_sources) <= 1 and unique_sources[0] == "Unspecified source":
        ax.scatter(
            x_exp,
            values_exp,
            color="black",
            marker="o",
            s=34,
            facecolors="none",
            label=f"{salt} data ({solvent_system})",
        )
    else:
        for src in unique_sources:
            mask = source_exp == src
            color = _color_for_source(str(src))
            src_label = f"{salt} data - {src}" if src != "Unspecified source" else f"{salt} data - Unspecified source"
            ax.scatter(
                x_exp[mask],
                values_exp[mask],
                color=color,
                marker="o",
                s=36,
                facecolors="none",
                linewidths=1.1,
                label=src_label,
            )

    for dataset_name in active_variants:
        cfg = DATASET_VARIANTS[dataset_name]
        ax.plot(
            x_grid,
            curves[dataset_name],
            color=str(cfg["color"]),
            linestyle=str(cfg["linestyle"]),
            linewidth=float(cfg["lw"]),
            label=f"{salt} {cfg['label']}",
        )

    ax.set_xlim(0.0, x_upper)
    ax.set_ylim(0.0, ymax)
    ax.set_xlabel(x_label, fontsize=AXIS_LABEL_SIZE)
    if quantity == "miac_m":
        ax.set_ylabel(r"mean ionic activity coefficient, $\gamma_{\pm}^{m,*}$", fontsize=AXIS_LABEL_SIZE)
    else:
        ax.set_ylabel(r"mean ionic activity coefficient, $\gamma_{\pm}^{*}$", fontsize=AXIS_LABEL_SIZE)
    ax.xaxis.label.set_color("black")
    ax.yaxis.label.set_color("black")
    ax.title.set_color("black")
    ax.tick_params(axis="both", labelsize=AXIS_TICK_SIZE)
    ax.tick_params(colors="black")
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color("black")
        spine.set_linewidth(1.0)

    comp_label = str(combo.get("comp_label", "")).strip()
    q_label = "MIAC_m" if quantity == "miac_m" else "MIAC"
    if comp_label:
        ax.set_title(f"{salt} in {solvent_system} at 298.15 K ({comp_label}) [{q_label}]")
    else:
        ax.set_title(f"{salt} in {solvent_system} at 298.15 K [{q_label}]")

    ax.grid(True, alpha=0.3, color="0.7")
    if show_legend:
        legend = ax.legend(fontsize=8)
        frame = legend.get_frame()
        frame.set_facecolor("white")
        frame.set_edgecolor("black")
        frame.set_alpha(1.0)
        for text in legend.get_texts():
            text.set_color("black")

    if save:
        if output_path is None:
            output_key = "output_miac_m" if quantity == "miac_m" else "output_miac"
            output_path = Path(combo[output_key])
        output_path.parent.mkdir(parents=True, exist_ok=True)
        fig.tight_layout()
        save_plot_figure(fig, output_path, dpi=220, bbox_inches=None)
        if not output_path.exists():
            raise FileNotFoundError(f"Expected plot was not written: {output_path}")

    if created_fig and close:
        plt.close(fig)

    return {
        "figure": fig,
        "axis": ax,
        "output_path": output_path,
        "quantity": quantity,
        "molality_exp": molal_exp,
        "x_exp": x_exp,
        "source_exp": source_exp,
        "values_exp": values_exp,
        "molality_grid": molal_grid,
        "x_grid": x_grid,
        "curves": curves,
        "active_variants": active_variants,
    }


def _grid_output_path(solvent_system: str) -> Path:
    return fits_plot_path("miac", solvent_system, "miac_m", f"miac_m_{solvent_system}_grid_3x3.png")


def plot_single_solvent_grid(
    solvent_system: str,
    combos: List[Dict[str, object]],
    close: bool = True,
    payload_map: Dict[Tuple[str, str], Dict[str, object]] | None = None,
) -> Path:
    combo_map = {
        str(combo["salt"]): combo
        for combo in combos
        if str(combo["solvent_system"]) == solvent_system and not combo.get("comp_signature")
    }

    grid_salts = GRID_SALTS[:2] if solvent_system == "ethanol" else GRID_SALTS
    nrows = len(grid_salts)
    fig_height = 7.8 if solvent_system == "ethanol" else 11.2
    fig, axes = plt.subplots(nrows, 3, figsize=(12.8, fig_height), sharex=False, sharey=False)
    fig.patch.set_facecolor("white")
    axes = np.atleast_2d(axes)

    legend_handles = None
    legend_labels = None
    for r, salt_row in enumerate(grid_salts):
        for c, salt in enumerate(salt_row):
            ax = axes[r, c]
            ax.set_facecolor("white")
            combo = combo_map.get(salt)
            if combo is None:
                ax.axis("off")
                continue

            payload = None
            if payload_map is not None:
                payload = payload_map.get((solvent_system, salt))
            plot_combo(
                combo,
                save=False,
                close=False,
                ax=ax,
                show_legend=False,
                quantity="miac_m",
                payload=payload,
            )
            if legend_handles is None:
                legend_handles, legend_labels = ax.get_legend_handles_labels()

    if legend_handles and legend_labels:
        legend = fig.legend(
            legend_handles,
            legend_labels,
            loc="lower center",
            ncol=min(5, len(legend_labels)),
            fontsize=8,
            frameon=True,
            bbox_to_anchor=(0.5, 0.01),
        )
        frame = legend.get_frame()
        frame.set_facecolor("white")
        frame.set_edgecolor("black")
        frame.set_alpha(1.0)
        for text in legend.get_texts():
            text.set_color("black")

    fig.suptitle(f"{solvent_system.capitalize()} MIAC_m fits at 298.15 K", color="black", fontsize=14)
    fig.tight_layout(rect=(0.0, 0.06, 1.0, 0.96))

    out = _grid_output_path(solvent_system)
    out.parent.mkdir(parents=True, exist_ok=True)
    save_plot_figure(fig, out, dpi=220, bbox_inches=None)
    if close:
        plt.close(fig)
    return out


def _render_combo_outputs(
    combo: Dict[str, object],
    quantities: List[Literal["miac_m", "miac"]],
    grid_points: int,
) -> Dict[str, object]:
    payload = prepare_combo_payload(combo, grid_points=grid_points)
    generated = []
    for quantity in quantities:
        result = plot_combo(combo, save=True, close=True, quantity=quantity, payload=payload)
        generated.append(Path(result["output_path"]))
    grid_payload = None
    if not combo.get("comp_signature") and str(combo["solvent_system"]) in {"water", "methanol", "ethanol"}:
        grid_payload = {
            "key": (str(combo["solvent_system"]), str(combo["salt"])),
            "payload": payload,
        }
    return {"generated": generated, "grid_payload": grid_payload}


def run_validate_miac_fits_v2() -> List[Path]:
    combos = discover_combos()
    quantities = _requested_quantities()
    grid_points = _requested_grid_points()
    workers = _requested_workers(len(combos))
    generated: List[Path] = []
    grid_payload_map: Dict[Tuple[str, str], Dict[str, object]] = {}

    if workers == 1:
        for combo in combos:
            result = _render_combo_outputs(combo, quantities, grid_points)
            generated.extend(Path(p) for p in result["generated"])
            grid_payload = result["grid_payload"]
            if grid_payload is not None:
                grid_payload_map[tuple(grid_payload["key"])] = dict(grid_payload["payload"])
    else:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            futures = [pool.submit(_render_combo_outputs, combo, quantities, grid_points) for combo in combos]
            for future in as_completed(futures):
                result = future.result()
                generated.extend(Path(p) for p in result["generated"])
                grid_payload = result["grid_payload"]
                if grid_payload is not None:
                    grid_payload_map[tuple(grid_payload["key"])] = dict(grid_payload["payload"])

    if "miac_m" in quantities:
        for solvent_system in ("water", "methanol", "ethanol"):
            if any(
                str(combo["solvent_system"]) == solvent_system and not combo.get("comp_signature") for combo in combos
            ):
                generated.append(
                    plot_single_solvent_grid(solvent_system, combos, close=True, payload_map=grid_payload_map)
                )

    generated = sorted(generated, key=lambda path: str(path))

    print("Dataset variants:")
    for dataset_name, cfg in DATASET_VARIANTS.items():
        print(f"- {dataset_name} -> {cfg['label']}")
    print(f"Workers: {workers}")
    print(f"Grid points: {grid_points}")
    print("Generated validation plots:")
    for path in generated:
        print(f"- {path}")
    return generated


def test_validate_miac_fits_v2() -> None:
    run_validate_miac_fits_v2()


if __name__ == "__main__":
    run_validate_miac_fits_v2()
