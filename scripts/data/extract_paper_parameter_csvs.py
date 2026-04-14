"""Extract per-paper binary interaction matrices and pure-component parameter CSVs.

Outputs:
- data/epcsaft_parameters/<dataset_key>/pure/<name>.csv
- data/epcsaft_parameters/<dataset_key>/mixed/binary_interaction/k_ij.csv
- optional:
  - data/epcsaft_parameters/<dataset_key>/mixed/binary_interaction/l_ij.csv
  - data/epcsaft_parameters/<dataset_key>/mixed/binary_interaction/khb_ij.csv

Target papers: 2005, 2008, 2014, 2020, 2021, 2025
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import shutil
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple


REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts._env import require_epcsaft_install

require_epcsaft_install()

PAPER_DIR = REPO_ROOT / "docs" / "papers" / "md"
OUT_BASE_DIR = REPO_ROOT / "data" / "epcsaft_parameters"
OUT_BINARY_DIR = OUT_BASE_DIR / "binary_interaction_parameters"
OUT_PURE_DIR = OUT_BASE_DIR / "pure_component_parameters"
LEGACY_FLAT_DIRS = [OUT_BINARY_DIR, OUT_PURE_DIR]
OPTIONAL_INTERACTION_TYPES = ("l_ij", "khb_ij")
PURE_FILENAME_BY_DATASET = {
    "2008_Held": "water.csv",
    "2012_Held": None,
}
COMPONENT_ORDER = [
    "water",
    "methanol",
    "ethanol",
    "Li+",
    "Na+",
    "K+",
    "Cl-",
    "Br-",
    "I-",
]
TARGET_COMPONENTS = set(COMPONENT_ORDER)

DISPLAY_NAME_MAP = {
    "water": "H2O",
    "methanol": "Methanol",
    "ethanol": "Ethanol",
    "Li+": "Li+",
    "Na+": "Na+",
    "K+": "K+",
    "Cl-": "Cl-",
    "Br-": "Br-",
    "I-": "I-",
}
DISPLAY_COMPONENT_ORDER = [DISPLAY_NAME_MAP[c] for c in COMPONENT_ORDER]
PURE_WIDE_COLUMNS = [
    "component",
    "m",
    "s",
    "e",
    "e_assoc",
    "vol_a",
    "assoc_scheme",
    "dipm",
    "dip_num",
    "z",
    "dielc",
    "d_born",
    "f_solv",
]
PURE_PARAM_MAP = {
    "m_seg": "m",
    "sigma": "s",
    "sigma_expr": "s",
    "u_over_kb": "e",
    "association_energy": "e_assoc",
    "association_volume": "vol_a",
    "d_born": "d_born",
    "f_k": "f_solv",
}

PAPER_FILES = {
    "2005_Cameretti": "Cameretti, Sadowski, Mollerup - 2005 - Modeling of Aqueous Electrolyte Solutions with Perturbed-Chai.md",
    "2008_Held": "Held, Cameretti, Sadowski - 2008 - Modeling aqueous electrolyte solutions. Part 1. Fully dissociated.md",
    "2014_Held": "Held et al. - 2014 - ePC-SAFT Revised.md",
    "2020_Bulow": "Bülow, Ascani, Held - 2020 - ePC-SAFT advanced - Part I Physical meaning of including a concentratio.md",
    "2021_Bulow": "Bülow, Ascani, Held - 2021 - ePC-SAFT advanced – Part II Application to Salt Solubility in Ionic and.md",
    "2025_Figiel": "Figiel, Yu, Held - 2025 - Predicting Thermodynamic Properties of Ions in Single Solvents and in Mixe.md",
}

EXPECTED_BASENAMES = list(PAPER_FILES.keys())
PAPER_KEY_BY_DATASET = {
    "2005_Cameretti": "2005",
    "2008_Held": "2008",
    "2014_Held": "2014",
    "2020_Bulow": "2020",
    "2021_Bulow": "2021",
    "2025_Figiel": "2025",
}
CANONICAL_USER_OPTIONS_SOURCE = "scripts/fits/validate_miac_fits.py"
RESOLVER_SOURCE = "epcsaft.parameters::_resolve_runtime_options"
RUNTIME_REQUIRED_KEYS = {
    "dielc_rule",
    "dielc_diff_mode",
    "d_ion_mode",
    "mu_DH_diff_mode",
    "mu_DH_comp_dep_rel_perm",
    "mu_DH_include_sum_term",
    "include_born_model",
    "d_Born_mode",
    "born_solvation_shell_model",
    "born_dielectric_saturation",
    "born_bulk_mode",
    "mu_born_diff_mode",
    "mu_born_comp_dep_rel_perm",
    "mu_born_include_sum_term",
    "mu_born_comp_dep_delta_d",
    "bjeruum_treatment",
    "debug",
}
RUNTIME_SENTINELS: Dict[str, Dict[str, Any]] = {
    "2005_Cameretti": {"include_born_model": False, "dielc_rule": 0, "d_ion_mode": 1},
    "2008_Held": {"include_born_model": False, "dielc_rule": 0, "d_ion_mode": 1},
    "2014_Held": {"include_born_model": False, "dielc_rule": 0, "d_ion_mode": 1},
    "2020_Bulow": {"include_born_model": True, "dielc_rule": 1, "d_Born_mode": 1, "d_ion_mode": 1},
    "2021_Bulow": {"include_born_model": True, "dielc_rule": 3, "d_Born_mode": 0, "d_ion_mode": 1},
    "2025_Figiel": {
        "include_born_model": True,
        "dielc_rule": 4,
        "dielc_diff_mode": 1,
        "d_Born_mode": 3,
        "born_solvation_shell_model": True,
        "born_dielectric_saturation": True,
        "mu_born_diff_mode": 1,
    },
}
DATASET_USER_OPTIONS: Dict[str, Dict[str, Any]] = {
    "2005_Cameretti": {
        "elec_model": {
            "rel_perm": {"rule": "constant", "differential_mode": "analytical"},
            "DH_model": {"d_ion_mode": 1, "bjeruum_treatment": False, "mu_DH_model": {"differential_mode": "analytical", "comp_dep_rel_perm": True, "include_sum_term": True}},
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
        },
        "debug": False,
    },
    "2008_Held": {
        "elec_model": {
            "rel_perm": {"rule": "constant", "differential_mode": "analytical"},
            "DH_model": {"d_ion_mode": 1, "bjeruum_treatment": False, "mu_DH_model": {"differential_mode": "analytical", "comp_dep_rel_perm": True, "include_sum_term": True}},
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
        },
        "debug": False,
    },
    "2014_Held": {
        "elec_model": {
            "rel_perm": {"rule": "constant", "differential_mode": "analytical"},
            "DH_model": {"d_ion_mode": 1, "bjeruum_treatment": False, "mu_DH_model": {"differential_mode": "analytical", "comp_dep_rel_perm": True, "include_sum_term": True}},
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
        },
        "debug": False,
    },
    "2020_Bulow": {
        "elec_model": {
            "rel_perm": {"rule": 1, "differential_mode": "analytical"},
            "DH_model": {"d_ion_mode": 1, "bjeruum_treatment": False, "mu_DH_model": {"differential_mode": "analytical", "comp_dep_rel_perm": True, "include_sum_term": True}},
            "include_born_model": True,
            "born_model": {
                "d_Born_mode": 1,
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
        },
        "debug": False,
    },
    "2021_Bulow": {
        "elec_model": {
            "rel_perm": {"rule": 3, "differential_mode": "analytical"},
            "DH_model": {"d_ion_mode": 1, "bjeruum_treatment": False, "mu_DH_model": {"differential_mode": "analytical", "comp_dep_rel_perm": True, "include_sum_term": True}},
            "include_born_model": True,
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
        },
        "debug": False,
    },
    "2025_Figiel": {
        "elec_model": {
            "rel_perm": {"rule": "empirical", "differential_mode": "numerical"},
            "DH_model": {"d_ion_mode": 1, "bjeruum_treatment": False, "mu_DH_model": {"differential_mode": "analytical", "comp_dep_rel_perm": True, "include_sum_term": True}},
            "include_born_model": True,
            "born_model": {
                "d_Born_mode": 3,
                "solvation_shell_model": True,
                "dielectric_saturation": True,
                "bulk_mode": "mix",
                "mu_born_model": {
                    "differential_mode": "numerical",
                    "comp_dep_rel_perm": True,
                    "include_sum_term": True,
                    "comp_dep_delta_d": True,
                },
            },
        },
        "debug": False,
    },
}
BINARY_TABLE_PRECEDENCE = {
    "2005_Cameretti": [],
    "2008_Held": [],
    "2014_Held": ["table2_ion_water", "table3_cation_anion"],
    "2020_Bulow": ["table2_ion_water", "table3_cation_anion"],
    "2021_Bulow": ["table3_ion_water", "table4_cation_anion", "table8_ion_organic"],
    "2025_Figiel": ["table4_cation_anion", "table5_ion_solvent"],
}


@dataclass
class PureRow:
    component: str
    parameter: str
    value: str
    unit: str
    source_table: str
    paper: str
    notes: str


def _clean_cell(cell: str) -> str:
    s = cell.strip()
    s = s.replace("\u2212", "-")
    s = s.replace("−", "-")
    s = s.replace("•", "*")
    s = s.replace("\\bullet", "*")
    s = s.replace("\\cdot", "*")
    s = s.replace("\\times", "*")
    s = re.sub(r"\$\{?\s*\^\{[^}]*\}\s*", "", s)
    s = re.sub(r"\$\{?\s*\\text\s*\{[^}]*\}\s*", "", s)
    s = re.sub(r"\$\{?\s*\\mathrm\s*\{", "", s)
    s = s.replace("\\mathrm", "")
    s = s.replace("\\text", "")
    s = s.replace("{", "")
    s = s.replace("}", "")
    s = s.replace("$", "")
    s = s.replace("\\", "")
    s = s.replace("~", "")
    s = re.sub(r"\s+", " ", s).strip()
    return s


def _is_missing_value(value: str) -> bool:
    s = _clean_cell(value).strip().lower()
    return s in {"", "-", "n.a.", "n.a", "na", "n/a"}


def _normalize_component(raw: str) -> Optional[str]:
    s = raw.strip()
    if not s:
        return None
    s_clean = _clean_cell(s)
    s_low = s_clean.lower()
    if s_low in {"water", "h2o"}:
        return "water"
    if s_low == "methanol":
        return "methanol"
    if s_low == "ethanol":
        return "ethanol"

    compact = re.sub(r"\s+", "", s_clean)
    if re.search(r"li\^?\+|li\+", compact, re.I):
        return "Li+"
    if re.search(r"na\^?\+|na\+", compact, re.I):
        return "Na+"
    if re.search(r"k\^?\+|k\+", compact, re.I):
        return "K+"
    if re.search(r"cl\^?-|cl-", compact, re.I):
        return "Cl-"
    if re.search(r"br\^?-|br-", compact, re.I):
        return "Br-"
    if re.search(r"i\^?-|i-", compact, re.I):
        return "I-"
    return None


def _normalize_value(raw: str) -> Optional[str]:
    if _is_missing_value(raw):
        return None
    s = _clean_cell(raw)
    s = re.sub(r"\^\([^)]+\)", "", s)
    s = re.sub(r"\^\d+", "", s)
    s = s.strip()
    if not s:
        return None
    # Keep formula text as-is if it has T, letters, or operators.
    if re.search(r"[A-Za-z]", s):
        return s
    try:
        num = float(s)
    except ValueError:
        return s
    if abs(num) < 1e-15:
        return "0"
    return format(num, ".15g")


def _markdown_table_after_anchor(lines: Sequence[str], anchor: str) -> List[List[str]]:
    anchor_idx = None
    anchor_low = anchor.lower()
    for idx, line in enumerate(lines):
        if line.strip().lower().startswith(anchor_low):
            anchor_idx = idx
            break
    if anchor_idx is None:
        for idx, line in enumerate(lines):
            if anchor_low in line.lower():
                anchor_idx = idx
                break
    if anchor_idx is None:
        raise ValueError(f"Could not find anchor: {anchor}")

    start_idx = None
    for idx in range(anchor_idx + 1, len(lines)):
        line = lines[idx].strip()
        if line.startswith("|"):
            start_idx = idx
            break
    if start_idx is None:
        raise ValueError(f"Could not find markdown table after anchor: {anchor}")

    rows: List[List[str]] = []
    for idx in range(start_idx, len(lines)):
        line = lines[idx].strip()
        if not line.startswith("|"):
            break
        cells = [c.strip() for c in line.strip("|").split("|")]
        if all(re.fullmatch(r":?-{2,}:?", c.replace(" ", "")) for c in cells):
            continue
        rows.append(cells)
    return rows


def _empty_matrix() -> Dict[str, Dict[str, str]]:
    matrix: Dict[str, Dict[str, str]] = {}
    for i in COMPONENT_ORDER:
        matrix[i] = {j: "0" for j in COMPONENT_ORDER}
    return matrix


def _is_zero(value: Optional[str]) -> bool:
    if value is None:
        return True
    s = value.strip()
    if s == "":
        return True
    try:
        return abs(float(s)) < 1e-15
    except ValueError:
        return False


def _set_pair(
    matrix: Dict[str, Dict[str, str]],
    comp_i: str,
    comp_j: str,
    value: Optional[str],
    allow_override: bool,
) -> None:
    if comp_i not in TARGET_COMPONENTS or comp_j not in TARGET_COMPONENTS:
        return
    if comp_i == comp_j or value is None:
        return
    current = matrix[comp_i][comp_j]
    if _is_zero(current) or allow_override:
        matrix[comp_i][comp_j] = value
        return
    if current != value and not _is_zero(value):
        raise ValueError(f"Conflicting values for {comp_i}/{comp_j}: {current} vs {value}")


def _enforce_symmetry(matrix: Dict[str, Dict[str, str]]) -> None:
    for i, comp_i in enumerate(COMPONENT_ORDER):
        for comp_j in COMPONENT_ORDER[i + 1 :]:
            a = matrix[comp_i][comp_j]
            b = matrix[comp_j][comp_i]
            if _is_zero(a) and not _is_zero(b):
                matrix[comp_i][comp_j] = b
            elif _is_zero(b) and not _is_zero(a):
                matrix[comp_j][comp_i] = a
            elif not _is_zero(a) and not _is_zero(b) and a != b:
                raise ValueError(f"Asymmetric conflict for {comp_i}/{comp_j}: {a} vs {b}")
    for comp in COMPONENT_ORDER:
        matrix[comp][comp] = "0"


def _display_component(component: str) -> str:
    return DISPLAY_NAME_MAP.get(component, component[:1].upper() + component[1:])


def _dataset_dir(dataset_key: str) -> Path:
    return OUT_BASE_DIR / dataset_key


def _dataset_pure_dir(dataset_key: str) -> Path:
    return _dataset_dir(dataset_key) / "pure"


def _dataset_pure_path(dataset_key: str) -> Optional[Path]:
    filename = PURE_FILENAME_BY_DATASET.get(dataset_key, "any_solvent.csv")
    if filename is None:
        return None
    return _dataset_pure_dir(dataset_key) / filename


def _dataset_binary_dir(dataset_key: str) -> Path:
    return _dataset_dir(dataset_key) / "mixed" / "binary_interaction"


def _matrix_has_parameter_data(matrix: Dict[str, Dict[str, str]]) -> bool:
    for i, comp_i in enumerate(COMPONENT_ORDER):
        for comp_j in COMPONENT_ORDER[i + 1 :]:
            if not _is_zero(matrix[comp_i][comp_j]):
                return True
    return False


def _set_cation_anion_kij_one(matrix: Dict[str, Dict[str, str]]) -> None:
    cations = ["Li+", "Na+", "K+"]
    anions = ["Cl-", "Br-", "I-"]
    for cat in cations:
        for an in anions:
            _set_pair(matrix, cat, an, "1.0", allow_override=True)

def _resolve_runtime_payload(canonical_options: Dict[str, Any]) -> Dict[str, Any]:
    from epcsaft.parameters import _resolve_runtime_options

    resolved = _resolve_runtime_options(canonical_options)
    return {
        "preset_key": resolved["preset_key"],
        "runtime_options": resolved["runtime"],
    }


def _build_user_options_payload(dataset_key: str, canonical_options: Dict[str, Any]) -> Dict[str, Any]:
    from epcsaft.parameters import minimize_user_options

    if dataset_key not in PAPER_KEY_BY_DATASET:
        raise ValueError(f"Missing paper key mapping for dataset {dataset_key}")
    return minimize_user_options(canonical_options)


def _assert_runtime_sentinels(dataset_key: str, runtime_options: Dict[str, Any]) -> None:
    for key, expected_value in RUNTIME_SENTINELS[dataset_key].items():
        actual_value = runtime_options.get(key)
        if actual_value != expected_value:
            raise ValueError(
                f"{dataset_key}: runtime sentinel mismatch for {key}: {actual_value} != {expected_value}"
            )


def _validate_user_options_payload(dataset_key: str, payload: Dict[str, Any]) -> None:
    if not isinstance(payload, dict):
        raise ValueError(f"{dataset_key}: user_options payload must be a dict")
    runtime_options = _resolve_runtime_payload(payload)["runtime_options"]
    missing = RUNTIME_REQUIRED_KEYS - set(runtime_options.keys())
    if missing:
        raise ValueError(
            f"{dataset_key}: runtime_options missing required keys: {sorted(missing)}"
        )
    _assert_runtime_sentinels(dataset_key, runtime_options)


def _write_user_options_json(path: Path, payload: Dict[str, Any]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2)
        handle.write("\n")


def _write_binary_csv(path: Path, matrix: Dict[str, Dict[str, str]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["component"] + DISPLAY_COMPONENT_ORDER)
        for row_comp in COMPONENT_ORDER:
            writer.writerow(
                [_display_component(row_comp)]
                + [matrix[row_comp][col_comp] for col_comp in COMPONENT_ORDER]
            )


def _write_pure_csv(path: Path, rows: Sequence[PureRow]) -> None:
    # Streamlined wide schema regardless of paper source.
    out_rows: Dict[str, Dict[str, str]] = {}
    for comp in COMPONENT_ORDER:
        display = _display_component(comp)
        out_rows[display] = {col: "" for col in PURE_WIDE_COLUMNS}
        out_rows[display]["component"] = display

    for row in rows:
        comp = row.component
        if comp not in TARGET_COMPONENTS:
            continue
        display = _display_component(comp)
        mapped = PURE_PARAM_MAP.get(row.parameter)
        if mapped:
            out_rows[display][mapped] = row.value
        if row.parameter == "association_sites" and row.value.strip() == "2":
            out_rows[display]["assoc_scheme"] = "2B"
        if row.parameter in {"association_energy", "association_volume"} and row.value.strip():
            if not out_rows[display]["assoc_scheme"]:
                out_rows[display]["assoc_scheme"] = "2B"

    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=PURE_WIDE_COLUMNS)
        writer.writeheader()
        for display in DISPLAY_COMPONENT_ORDER:
            writer.writerow(out_rows[display])


def _extract_2005(lines: Sequence[str], paper_key: str) -> Tuple[Dict[str, Dict[str, str]], List[PureRow]]:
    matrix = _empty_matrix()
    pure_rows: List[PureRow] = []

    t1 = _markdown_table_after_anchor(lines, "Table 1. ePC-SAFT Parameters for Water")
    # Vertical table: parameter | value
    param_map = {
        "segment number": ("m_seg", "-"),
        "segment diameter [aa]": ("sigma", "A"),
        "dispersion energy [k]": ("u_over_kb", "K"),
        "association sites": ("association_sites", "-"),
        "association energy [k]": ("association_energy", "K"),
        "association volume": ("association_volume", "-"),
    }
    for row in t1:
        if len(row) < 2:
            continue
        p_raw = _clean_cell(row[0]).lower()
        v_raw = _normalize_value(row[1])
        if v_raw is None:
            continue
        key = p_raw.replace("å", "a").replace("Å", "a")
        if key in param_map:
            p_name, unit = param_map[key]
            pure_rows.append(
                PureRow("water", p_name, v_raw, unit, "Table 1", paper_key, "")
            )

    t2 = _markdown_table_after_anchor(lines, "Table 2. Optimized ePC-SAFT Parameters for Alkali Halide Ions")
    for row in t2[1:]:
        if len(row) < 4:
            continue
        comp = _normalize_component(row[0])
        if comp not in TARGET_COMPONENTS:
            continue
        sigma = _normalize_value(row[2])
        uval = _normalize_value(row[3])
        if sigma is not None:
            pure_rows.append(PureRow(comp, "sigma", sigma, "A", "Table 2", paper_key, ""))
        if uval is not None:
            pure_rows.append(PureRow(comp, "u_over_kb", uval, "K", "Table 2", paper_key, ""))

    _set_cation_anion_kij_one(matrix)
    return matrix, pure_rows


def _extract_2008(lines: Sequence[str], paper_key: str) -> Tuple[Dict[str, Dict[str, str]], List[PureRow]]:
    matrix = _empty_matrix()
    pure_rows: List[PureRow] = []

    t1 = _markdown_table_after_anchor(lines, "Table 1")
    for row in t1[1:]:
        if len(row) < 4:
            continue
        param = _clean_cell(row[0]).lower()
        abbrev = _clean_cell(row[2]).lower()
        value = _normalize_value(row[3])
        if value is None:
            continue
        if "segment number" in param:
            pure_rows.append(PureRow("water", "m_seg", value, "-", "Table 1", paper_key, ""))
        elif "segment diameter" in param and "sigma" in abbrev:
            pure_rows.append(PureRow("water", "sigma", value, "A", "Table 1", paper_key, ""))
        elif abbrev in {"t_dep, 1", "t_dep, 2", "t_dep, 3", "t_dep, 4"}:
            pure_rows.append(PureRow("water", abbrev.replace(", ", "_"), value, row[1].strip(), "Table 1", paper_key, ""))
        elif "dispersion energy" in param:
            pure_rows.append(PureRow("water", "u_over_kb", value, "K", "Table 1", paper_key, ""))
        elif "association sites" in param:
            pure_rows.append(PureRow("water", "association_sites", value, "-", "Table 1", paper_key, ""))
        elif "association energy" in param:
            pure_rows.append(PureRow("water", "association_energy", value, "K", "Table 1", paper_key, ""))
        elif "association volume" in param:
            pure_rows.append(PureRow("water", "association_volume", value, "-", "Table 1", paper_key, ""))

    t2 = _markdown_table_after_anchor(lines, "Table 2")
    # Dual-sided table in one row: Ion sigma u | Ion sigma u
    for row in t2[2:]:
        if len(row) < 6:
            continue
        left_comp = _normalize_component(row[0])
        left_sigma = _normalize_value(row[1])
        left_u = _normalize_value(row[2])
        right_comp = _normalize_component(row[3])
        right_sigma = _normalize_value(row[4])
        right_u = _normalize_value(row[5])

        if left_comp in TARGET_COMPONENTS:
            if left_sigma is not None:
                pure_rows.append(PureRow(left_comp, "sigma", left_sigma, "A", "Table 2", paper_key, ""))
            if left_u is not None:
                pure_rows.append(PureRow(left_comp, "u_over_kb", left_u, "K", "Table 2", paper_key, ""))

        if right_comp in TARGET_COMPONENTS:
            if right_sigma is not None:
                pure_rows.append(PureRow(right_comp, "sigma", right_sigma, "A", "Table 2", paper_key, ""))
            if right_u is not None:
                pure_rows.append(PureRow(right_comp, "u_over_kb", right_u, "K", "Table 2", paper_key, ""))

    _set_cation_anion_kij_one(matrix)
    return matrix, pure_rows


def _extract_2014(text: str, paper_key: str) -> Tuple[Dict[str, Dict[str, str]], List[PureRow]]:
    matrix = _empty_matrix()
    pure_rows: List[PureRow] = []

    # Table 1 water parameters
    m_match = re.search(r"\\hline \$m_\{i\}\^\{\\text \{seg \}\}\$ & ([^&]+) &", text)
    u_match = re.search(r"\\hline \$u_\{i\} .*?& ([^&]+) &", text)
    n_match = re.search(r"\\hline \$N_\{i\}\$ & ([^&]+) &", text)
    e_assoc_match = re.search(r"\\hline \$\\varepsilon\^\{A_\{i\} B_\{i\}\} .*?& ([^&]+) &", text)
    v_assoc_match = re.search(r"\\hline \$k\^\{A_\{i\} B_\{i\}\}\$ & ([^&]+) &", text)
    sigma_expr_match = re.search(r"The expression \$\\sigma=([^$]+)\$ was used", text)

    if m_match:
        pure_rows.append(PureRow("water", "m_seg", _normalize_value(m_match.group(1)) or "", "-", "Table 1", paper_key, ""))
    if u_match:
        pure_rows.append(PureRow("water", "u_over_kb", _normalize_value(u_match.group(1)) or "", "K", "Table 1", paper_key, ""))
    if n_match:
        pure_rows.append(PureRow("water", "association_sites", _normalize_value(n_match.group(1)) or "", "-", "Table 1", paper_key, ""))
    if e_assoc_match:
        pure_rows.append(PureRow("water", "association_energy", _normalize_value(e_assoc_match.group(1)) or "", "K", "Table 1", paper_key, ""))
    if v_assoc_match:
        pure_rows.append(PureRow("water", "association_volume", _normalize_value(v_assoc_match.group(1)) or "", "-", "Table 1", paper_key, ""))
    if sigma_expr_match:
        pure_rows.append(PureRow("water", "sigma_expr", _clean_cell(sigma_expr_match.group(1)), "A", "Table 1", paper_key, "footnote a"))

    lines = text.splitlines()

    # Table 2 rows for strategy 1 + strategy 2
    in_table2 = False
    for line in lines:
        s = line.strip()
        if "Table 2 - ePG-SAFT parameters for ions" in s:
            in_table2 = True
            continue
        if in_table2 and s.startswith("\\caption{Table 3"):
            in_table2 = False
            break
        if not in_table2:
            continue
        if not (s.startswith("\\hline") and "&" in s and "$" in s):
            continue
        if "Univalent" in s or "Bivalent" in s or "Bi-/trivalent" in s or "Ion &" in s:
            continue
        row = s
        if row.startswith("\\hline"):
            row = row[len("\\hline") :].strip()
        if row.endswith("\\\\"):
            row = row[:-2].strip()
        cells = [c.strip() for c in row.split("&")]
        if len(cells) < 6:
            continue
        comp = _normalize_component(cells[0])
        if comp not in TARGET_COMPONENTS:
            continue
        s1_sigma = _normalize_value(cells[1])
        s1_u = _normalize_value(cells[2])
        s2_sigma = _normalize_value(cells[3])
        s2_u = _normalize_value(cells[4])
        k_ion_water = _normalize_value(cells[5])
        if s1_sigma:
            pure_rows.append(PureRow(comp, "sigma", s1_sigma, "A", "Table 2", paper_key, "strategy 1"))
        if s1_u:
            pure_rows.append(PureRow(comp, "u_over_kb", s1_u, "K", "Table 2", paper_key, "strategy 1"))
        if s2_sigma:
            pure_rows.append(PureRow(comp, "sigma", s2_sigma, "A", "Table 2", paper_key, "strategy 2"))
        if s2_u:
            pure_rows.append(PureRow(comp, "u_over_kb", s2_u, "K", "Table 2", paper_key, "strategy 2"))
        _set_pair(matrix, comp, "water", k_ion_water, allow_override=True)

    # Footnotes for temperature-dependent kij
    na_formula = re.search(
        r"For water\s*\$?\s*/\s*\\mathrm\{Na\}\^\{\+\}\$?.*?k.*?=([^$]+)\$",
        text,
        re.I,
    )
    k_formula = re.search(
        r"For water\s*\$?\s*/\s*\\mathrm\{K\}\^\{\+\}\$?.*?k.*?=([^$]+)\$",
        text,
        re.I,
    )
    if na_formula:
        _set_pair(matrix, "Na+", "water", _clean_cell(na_formula.group(1)), allow_override=True)
    if k_formula:
        _set_pair(matrix, "K+", "water", _clean_cell(k_formula.group(1)), allow_override=True)

    # Table 3 cation-anion matrix
    in_table3 = False
    cations: List[str] = []
    for line in lines:
        s = line.strip()
        if "Table 3 - Binary" in s:
            in_table3 = True
            continue
        if in_table3 and s.startswith("\\caption{Table 4"):
            in_table3 = False
            break
        if not in_table3:
            continue
        if not (s.startswith("\\hline") and "&" in s and "$" in s):
            continue
        row = s
        if row.startswith("\\hline"):
            row = row[len("\\hline") :].strip()
        if row.endswith("\\\\"):
            row = row[:-2].strip()
        cells = [c.strip() for c in row.split("&")]
        if not cells:
            continue
        first = cells[0]
        if "k_{\\text {cation-anion }}" in first:
            cations = [(_normalize_component(c) or "") for c in cells[1:]]
            continue
        anion = _normalize_component(first)
        if anion not in TARGET_COMPONENTS:
            continue
        for idx, val in enumerate(cells[1:]):
            if idx >= len(cations):
                break
            cation = cations[idx]
            if cation not in TARGET_COMPONENTS:
                continue
            _set_pair(matrix, cation, anion, _normalize_value(val), allow_override=True)

    _enforce_symmetry(matrix)
    return matrix, pure_rows


def _extract_2020(lines: Sequence[str], paper_key: str) -> Tuple[Dict[str, Dict[str, str]], List[PureRow]]:
    matrix = _empty_matrix()
    pure_rows: List[PureRow] = []

    t1 = _markdown_table_after_anchor(lines, "Table 1")
    for row in t1[1:]:
        if len(row) < 7:
            continue
        comp = _normalize_component(row[0])
        if comp not in TARGET_COMPONENTS:
            continue
        m_seg = _normalize_value(row[1])
        sigma = _normalize_value(row[2])
        u = _normalize_value(row[3])
        e_assoc = _normalize_value(row[4])
        v_assoc = _normalize_value(row[5])
        if m_seg:
            pure_rows.append(PureRow(comp, "m_seg", m_seg, "-", "Table 1", paper_key, ""))
        if sigma and sigma != "*":
            pure_rows.append(PureRow(comp, "sigma", sigma, "A", "Table 1", paper_key, ""))
        if u:
            pure_rows.append(PureRow(comp, "u_over_kb", u, "K", "Table 1", paper_key, ""))
        if e_assoc:
            pure_rows.append(PureRow(comp, "association_energy", e_assoc, "K", "Table 1", paper_key, ""))
        if v_assoc:
            pure_rows.append(PureRow(comp, "association_volume", v_assoc, "-", "Table 1", paper_key, ""))

    for line in lines:
        if line.strip().startswith("* ") and "sigma=" in line:
            pure_rows.append(PureRow("water", "sigma_expr", _clean_cell(line.strip()[2:]), "A", "Table 1", paper_key, "footnote"))
            break

    t2 = _markdown_table_after_anchor(lines, "Table 2")
    for row in t2[1:]:
        if len(row) < 4:
            continue
        comp = _normalize_component(row[0])
        if comp not in TARGET_COMPONENTS:
            continue
        sigma = _normalize_value(row[1])
        u = _normalize_value(row[2])
        k_ion_water = _normalize_value(row[3])
        if sigma:
            pure_rows.append(PureRow(comp, "sigma", sigma, "A", "Table 2", paper_key, ""))
        if u:
            pure_rows.append(PureRow(comp, "u_over_kb", u, "K", "Table 2", paper_key, ""))
        _set_pair(matrix, comp, "water", k_ion_water, allow_override=True)

    t3 = _markdown_table_after_anchor(lines, "Table 3")
    header = t3[0]
    cations = [_normalize_component(c) for c in header[1:]]
    for row in t3[1:]:
        anion = _normalize_component(row[0])
        if anion not in TARGET_COMPONENTS:
            continue
        for idx, cell in enumerate(row[1:]):
            if idx >= len(cations):
                break
            cation = cations[idx]
            if cation not in TARGET_COMPONENTS:
                continue
            _set_pair(matrix, cation, anion, _normalize_value(cell), allow_override=True)

    _enforce_symmetry(matrix)
    return matrix, pure_rows


def _extract_2021(lines: Sequence[str], paper_key: str) -> Tuple[Dict[str, Dict[str, str]], List[PureRow]]:
    matrix = _empty_matrix()
    pure_rows: List[PureRow] = []

    t1 = _markdown_table_after_anchor(lines, "Table 1")
    for row in t1[1:]:
        if len(row) < 7:
            continue
        comp = _normalize_component(row[0])
        if comp not in TARGET_COMPONENTS:
            continue
        m_seg = _normalize_value(row[1])
        sigma = _normalize_value(row[2])
        u = _normalize_value(row[3])
        e_assoc = _normalize_value(row[4])
        v_assoc = _normalize_value(row[5])
        if m_seg:
            pure_rows.append(PureRow(comp, "m_seg", m_seg, "-", "Table 1", paper_key, ""))
        if sigma:
            pure_rows.append(PureRow(comp, "sigma", sigma, "A", "Table 1", paper_key, ""))
        if u:
            pure_rows.append(PureRow(comp, "u_over_kb", u, "K", "Table 1", paper_key, ""))
        if e_assoc:
            pure_rows.append(PureRow(comp, "association_energy", e_assoc, "K", "Table 1", paper_key, ""))
        if v_assoc:
            pure_rows.append(PureRow(comp, "association_volume", v_assoc, "-", "Table 1", paper_key, ""))

    t3 = _markdown_table_after_anchor(lines, "Table 3")
    for row in t3[1:]:
        if len(row) < 4:
            continue
        comp = _normalize_component(row[0])
        if comp not in TARGET_COMPONENTS:
            continue
        sigma = _normalize_value(row[1])
        u = _normalize_value(row[2])
        k_ion_water = _normalize_value(row[3])
        if sigma:
            pure_rows.append(PureRow(comp, "sigma", sigma, "A", "Table 3", paper_key, ""))
        if u:
            pure_rows.append(PureRow(comp, "u_over_kb", u, "K", "Table 3", paper_key, ""))
        _set_pair(matrix, comp, "water", k_ion_water, allow_override=True)

    t4 = _markdown_table_after_anchor(lines, "Table 4")
    cations = [_normalize_component(c) for c in t4[0][1:]]
    for row in t4[1:]:
        anion = _normalize_component(row[0])
        if anion not in TARGET_COMPONENTS:
            continue
        for idx, cell in enumerate(row[1:]):
            if idx >= len(cations):
                break
            cation = cations[idx]
            if cation not in TARGET_COMPONENTS:
                continue
            _set_pair(matrix, cation, anion, _normalize_value(cell), allow_override=True)

    t8 = _markdown_table_after_anchor(lines, "Table 8")
    ions = [_normalize_component(c) for c in t8[0][1:]]
    for row in t8[1:]:
        solvent = _normalize_component(row[0])
        if solvent not in TARGET_COMPONENTS:
            continue
        for idx, cell in enumerate(row[1:]):
            if idx >= len(ions):
                break
            ion = ions[idx]
            if ion not in TARGET_COMPONENTS:
                continue
            _set_pair(matrix, solvent, ion, _normalize_value(cell), allow_override=True)

    _enforce_symmetry(matrix)
    return matrix, pure_rows


def _extract_2025(lines: Sequence[str], paper_key: str) -> Tuple[Dict[str, Dict[str, str]], List[PureRow]]:
    matrix = _empty_matrix()
    pure_rows: List[PureRow] = []

    t2 = _markdown_table_after_anchor(lines, "Table 2.")
    header = [h.strip().lower() for h in t2[0][1:]]
    params = [row[0].strip() for row in t2[1:]]
    values_by_param = {row[0].strip(): row[1:] for row in t2[1:]}
    solvents = [c for c in header if c in {"water", "methanol", "ethanol"}]
    for solvent in solvents:
        idx = header.index(solvent)
        for p in params:
            val = _normalize_value(values_by_param[p][idx])
            if val is None:
                continue
            p_clean = _clean_cell(p).lower()
            if p_clean.startswith("m_"):
                pure_rows.append(PureRow(solvent, "m_seg", val, "-", "Table 2", paper_key, ""))
            elif p_clean.startswith("sigma"):
                if val != "b":
                    pure_rows.append(PureRow(solvent, "sigma", val, "A", "Table 2", paper_key, ""))
            elif p_clean.startswith("u_"):
                pure_rows.append(PureRow(solvent, "u_over_kb", val, "K", "Table 2", paper_key, ""))
            elif "varepsilon" in p_clean:
                pure_rows.append(PureRow(solvent, "association_energy", val, "K", "Table 2", paper_key, ""))
            elif "kappa" in p_clean:
                pure_rows.append(PureRow(solvent, "association_volume", val, "-", "Table 2", paper_key, ""))
            elif p_clean == "f_k":
                pure_rows.append(PureRow(solvent, "f_k", val, "-", "Table 2", paper_key, ""))

    for line in lines:
        if "sigma=2.7927+" in line:
            pure_rows.append(PureRow("water", "sigma_expr", _clean_cell(line), "A", "Table 2", paper_key, "footnote b"))
            break

    t3 = _markdown_table_after_anchor(lines, "Table 3.")
    for row in t3[1:]:
        if len(row) < 4:
            continue
        comp = _normalize_component(row[0])
        if comp not in TARGET_COMPONENTS:
            continue
        sigma = _normalize_value(row[1])
        u = _normalize_value(row[2])
        d_born = _normalize_value(row[3])
        if sigma:
            pure_rows.append(PureRow(comp, "sigma", sigma, "A", "Table 3", paper_key, ""))
        if u:
            pure_rows.append(PureRow(comp, "u_over_kb", u, "K", "Table 3", paper_key, ""))
        if d_born:
            pure_rows.append(PureRow(comp, "d_born", d_born, "A", "Table 3", paper_key, ""))

    t4 = _markdown_table_after_anchor(lines, "Table 4.")
    cations = [_normalize_component(c) for c in t4[0][1:]]
    for row in t4[1:]:
        anion = _normalize_component(row[0])
        if anion not in TARGET_COMPONENTS:
            continue
        for idx, cell in enumerate(row[1:]):
            if idx >= len(cations):
                break
            cation = cations[idx]
            if cation not in TARGET_COMPONENTS:
                continue
            _set_pair(matrix, cation, anion, _normalize_value(cell), allow_override=True)

    t5 = _markdown_table_after_anchor(lines, "Table 5.")
    cols = [_normalize_component(c) for c in t5[0][1:]]
    for row in t5[1:]:
        row_comp = _normalize_component(row[0])
        if row_comp not in TARGET_COMPONENTS:
            continue
        for idx, cell in enumerate(row[1:]):
            if idx >= len(cols):
                break
            col_comp = cols[idx]
            if col_comp not in TARGET_COMPONENTS:
                continue
            _set_pair(matrix, row_comp, col_comp, _normalize_value(cell), allow_override=True)

    _enforce_symmetry(matrix)
    return matrix, pure_rows


def _extract_optional_interactions(
    paper_key: str,
    text: str,
) -> Dict[str, Dict[str, Dict[str, str]]]:
    # Placeholder routing for optional interaction tables.
    # Under the current scoped components, no l_ij or khb_ij matrices are emitted.
    _ = paper_key
    _ = text
    return {}


def _extract_for_paper(
    paper_key: str,
    text: str,
) -> Tuple[Dict[str, Dict[str, str]], List[PureRow], Dict[str, Dict[str, Dict[str, str]]]]:
    lines = text.splitlines()
    if paper_key == "2005_Cameretti":
        matrix, pure_rows = _extract_2005(lines, paper_key)
    elif paper_key == "2008_Held":
        matrix, pure_rows = _extract_2008(lines, paper_key)
    elif paper_key == "2014_Held":
        matrix, pure_rows = _extract_2014(text, paper_key)
    elif paper_key == "2020_Bulow":
        matrix, pure_rows = _extract_2020(lines, paper_key)
    elif paper_key == "2021_Bulow":
        matrix, pure_rows = _extract_2021(lines, paper_key)
    elif paper_key == "2025_Figiel":
        matrix, pure_rows = _extract_2025(lines, paper_key)
    else:
        raise ValueError(f"Unsupported paper key: {paper_key}")

    optional_interactions = _extract_optional_interactions(paper_key, text)
    return matrix, pure_rows, optional_interactions

def _validate_outputs(
    data_by_paper: Dict[str, Tuple[Dict[str, Dict[str, str]], List[PureRow], Dict[str, Dict[str, Dict[str, str]]]]],
    user_options_by_dataset: Dict[str, Dict[str, Any]],
    wrote_files: bool,
) -> None:
    # In-memory validations.
    for paper_key, (matrix, rows, optional_interactions) in data_by_paper.items():
        if len(matrix) != len(COMPONENT_ORDER):
            raise ValueError(f"{paper_key}: binary matrix row count mismatch")
        for row_comp in COMPONENT_ORDER:
            if row_comp not in matrix:
                raise ValueError(f"{paper_key}: missing row component {row_comp}")
            if len(matrix[row_comp]) != len(COMPONENT_ORDER):
                raise ValueError(f"{paper_key}: binary matrix column count mismatch for {row_comp}")
            if matrix[row_comp][row_comp] != "0":
                raise ValueError(f"{paper_key}: diagonal is not zero for {row_comp}")
        for i, comp_i in enumerate(COMPONENT_ORDER):
            for comp_j in COMPONENT_ORDER[i + 1 :]:
                if matrix[comp_i][comp_j] != matrix[comp_j][comp_i]:
                    raise ValueError(f"{paper_key}: asymmetry at {comp_i}/{comp_j}")
        if paper_key in {"2005_Cameretti", "2008_Held"}:
            for cat in ("Li+", "Na+", "K+"):
                for an in ("Cl-", "Br-", "I-"):
                    if matrix[cat][an] != "1.0" or matrix[an][cat] != "1.0":
                        raise ValueError(f"{paper_key}: expected cation-anion k_ij=1.0 for {cat}/{an}")

        expected_cols = {"component", "parameter", "value", "unit", "source_table", "paper", "notes"}
        for row in rows:
            fields = set(PureRow.__dataclass_fields__.keys())  # type: ignore[attr-defined]
            if fields != expected_cols:
                raise ValueError("PureRow schema mismatch")
            if row.component not in TARGET_COMPONENTS:
                raise ValueError(f"{paper_key}: pure row has non-target component {row.component}")

        for interaction_type, opt_matrix in optional_interactions.items():
            if interaction_type not in OPTIONAL_INTERACTION_TYPES:
                raise ValueError(f"{paper_key}: unsupported optional interaction type {interaction_type}")
            if not _matrix_has_parameter_data(opt_matrix):
                raise ValueError(f"{paper_key}: optional {interaction_type} present but has no scoped data")

    for dataset_key in EXPECTED_BASENAMES:
        if dataset_key not in user_options_by_dataset:
            raise ValueError(f"Missing user_options payload for dataset {dataset_key}")
        _validate_user_options_payload(dataset_key, user_options_by_dataset[dataset_key])

    # Sentinel checks.
    m2020 = data_by_paper["2020_Bulow"][0]
    if m2020["Cl-"]["Li+"] != "0.669":
        raise ValueError(f"2020 sentinel mismatch Cl-/Li+: {m2020['Cl-']['Li+']}")
    m2014 = data_by_paper["2014_Held"][0]
    if m2014["Cl-"]["Li+"] != "0.669":
        raise ValueError(f"2014 sentinel mismatch Cl-/Li+: {m2014['Cl-']['Li+']}")
    m2021 = data_by_paper["2021_Bulow"][0]
    if m2021["methanol"]["Na+"] != "-0.31":
        raise ValueError(f"2021 sentinel mismatch methanol/Na+: {m2021['methanol']['Na+']}")
    m2025 = data_by_paper["2025_Figiel"][0]
    if m2025["water"]["methanol"] != "-0.0878":
        raise ValueError(f"2025 sentinel mismatch water/methanol: {m2025['water']['methanol']}")

    if wrote_files:
        for dataset_key in EXPECTED_BASENAMES:
            ds_dir = _dataset_dir(dataset_key)
            if not ds_dir.is_dir():
                raise ValueError(f"Missing dataset directory: {ds_dir}")
            pure_dir = _dataset_pure_dir(dataset_key)
            if not pure_dir.is_dir():
                raise ValueError(f"Missing pure directory for {dataset_key}")
            pure_path = _dataset_pure_path(dataset_key)
            if pure_path is not None and not pure_path.exists():
                raise ValueError(f"Missing pure parameter file for {dataset_key}: {pure_path.name}")
            if pure_path is None and not any(pure_dir.glob("*.csv")):
                raise ValueError(f"Missing explicit pure parameter files for {dataset_key}")
            bin_dir = _dataset_binary_dir(dataset_key)
            if not bin_dir.is_dir():
                raise ValueError(f"Missing mixed/binary_interaction directory for {dataset_key}")
            kij_path = bin_dir / "k_ij.csv"
            if not kij_path.exists():
                raise ValueError(f"Missing k_ij.csv for {dataset_key}")
            user_options_path = ds_dir / "user_options.json"
            if not user_options_path.exists():
                raise ValueError(f"Missing user_options.json for {dataset_key}")

            optional_present = data_by_paper[dataset_key][2]
            for interaction_type in OPTIONAL_INTERACTION_TYPES:
                path = bin_dir / f"{interaction_type}.csv"
                should_exist = interaction_type in optional_present and _matrix_has_parameter_data(optional_present[interaction_type])
                if should_exist and not path.exists():
                    raise ValueError(f"Missing optional {interaction_type}.csv for {dataset_key}")
                if (not should_exist) and path.exists():
                    raise ValueError(f"Unexpected optional {interaction_type}.csv for {dataset_key}")

            file_payload = json.loads(user_options_path.read_text(encoding="utf-8"))
            _validate_user_options_payload(dataset_key, file_payload)
            if file_payload != user_options_by_dataset[dataset_key]:
                raise ValueError(f"{dataset_key}: on-disk user_options.json differs from generated payload")

        for legacy_dir in LEGACY_FLAT_DIRS:
            if legacy_dir.exists():
                raise ValueError(f"Legacy flat directory should be removed: {legacy_dir}")


def build(check_only: bool) -> None:
    OUT_BASE_DIR.mkdir(parents=True, exist_ok=True)
    if set(EXPECTED_BASENAMES) != set(DATASET_USER_OPTIONS.keys()):
        raise ValueError("DATASET_USER_OPTIONS keys must match EXPECTED_BASENAMES exactly")

    user_options_by_dataset: Dict[str, Dict[str, Any]] = {}
    for dataset_key in EXPECTED_BASENAMES:
        user_options_by_dataset[dataset_key] = _build_user_options_payload(
            dataset_key,
            DATASET_USER_OPTIONS[dataset_key],
        )

    extracted: Dict[str, Tuple[Dict[str, Dict[str, str]], List[PureRow], Dict[str, Dict[str, Dict[str, str]]]]] = {}
    for base_name, md_file in PAPER_FILES.items():
        src = PAPER_DIR / md_file
        if not src.exists():
            raise FileNotFoundError(f"Missing paper markdown file: {src}")
        text = src.read_text(encoding="utf-8")
        matrix, pure_rows, optional_interactions = _extract_for_paper(base_name, text)
        _enforce_symmetry(matrix)
        for opt_matrix in optional_interactions.values():
            _enforce_symmetry(opt_matrix)
        extracted[base_name] = (matrix, pure_rows, optional_interactions)

    if not check_only:
        for base_name, (matrix, pure_rows, optional_interactions) in extracted.items():
            ds_dir = _dataset_dir(base_name)
            pure_dir = _dataset_pure_dir(base_name)
            bin_dir = _dataset_binary_dir(base_name)
            ds_dir.mkdir(parents=True, exist_ok=True)
            pure_dir.mkdir(parents=True, exist_ok=True)
            bin_dir.mkdir(parents=True, exist_ok=True)

            pure_path = _dataset_pure_path(base_name)
            if pure_path is not None:
                _write_pure_csv(pure_path, pure_rows)
            _write_binary_csv(bin_dir / "k_ij.csv", matrix)
            _write_user_options_json(ds_dir / "user_options.json", user_options_by_dataset[base_name])

            for interaction_type in OPTIONAL_INTERACTION_TYPES:
                maybe_existing = bin_dir / f"{interaction_type}.csv"
                if maybe_existing.exists():
                    maybe_existing.unlink()

            for interaction_type, opt_matrix in optional_interactions.items():
                if interaction_type not in OPTIONAL_INTERACTION_TYPES:
                    continue
                if _matrix_has_parameter_data(opt_matrix):
                    _write_binary_csv(bin_dir / f"{interaction_type}.csv", opt_matrix)

        for legacy_dir in LEGACY_FLAT_DIRS:
            if legacy_dir.exists():
                shutil.rmtree(legacy_dir)

    _validate_outputs(extracted, user_options_by_dataset, wrote_files=not check_only)

    mode = "check-only" if check_only else "write"
    print(f"Extraction complete ({mode}). Processed {len(extracted)} papers.")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Extract per-paper binary interaction matrices and pure-component parameter CSVs."
    )
    parser.add_argument(
        "--check",
        action="store_true",
        help="Run extraction and validations without writing CSV files.",
    )
    args = parser.parse_args()
    build(check_only=args.check)


if __name__ == "__main__":
    main()


