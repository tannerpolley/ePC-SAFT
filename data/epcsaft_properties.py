"""Dataset-driven parameter loader for PC-SAFT/ePC-SAFT.

This v2 module reads parameter sets from
``data/pcsaft_parameters/<dataset>/`` and returns dictionaries in the same
shape expected by ``pcsaft`` runtime calls.

Public API:
    - get_prop_dict(dataset_name, species, x, T, user_options=None)
    - available_datasets()
    - _resolve_runtime_options(user_options)
    - molality_to_molefraction(...)
    - molefraction_to_molality(...)
"""

from __future__ import annotations

import copy
import csv
import json
import math
import re
from pathlib import Path
from typing import Dict, Iterable

import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[1]
DATASET_ROOT = REPO_ROOT / "data" / "pcsaft_parameters"

BASE_KEYS = ["MW", "m", "s", "e", "e_assoc", "vol_a", "assoc_scheme", "dipm", "dip_num", "z", "dielc"]
OPTIONAL_KEYS = ["d_born", "f_solv"]

COMPONENT_ALIASES = {
    "H2O-2B-Li": "H2O",
    "H2O-2B-NaCl": "H2O",
    "H2O-Salt-2001": "H2O",
    "Water": "H2O",
    "water": "H2O",
    "methanol": "Methanol",
    "ethanol": "Ethanol",
    "butanol": "Butanol",
    "propanol": "Propanol",
}

_COMPONENT_DEFAULTS = {
    "H+": {"MW": 1.008e-3, "z": 1.0, "d_born": 1.0},
    "Li+": {"MW": 6.94e-3, "z": 1.0, "d_born": 2.784},
    "Na+": {"MW": 22.98e-3, "z": 1.0, "d_born": 3.445},
    "K+": {"MW": 39.10e-3, "z": 1.0, "d_born": 4.150},
    "NH4+": {"MW": 18.038e-3, "z": 1.0, "d_born": 3.0},
    "Cl-": {"MW": 35.45e-3, "z": -1.0, "d_born": 4.100},
    "Br-": {"MW": 79.90e-3, "z": -1.0, "d_born": 4.480},
    "I-": {"MW": 126.90e-3, "z": -1.0, "d_born": 4.985},
    "H2O": {
        "MW": 18.01528e-3,
        "m": 1.2047,
        "s": lambda T: 2.7927 + 10.11 * np.exp(-0.01775 * T) - 1.417 * np.exp(-0.01146 * T),
        "e": 353.95,
        "e_assoc": 2425.7,
        "vol_a": 0.04509,
        "assoc_scheme": "2B",
        "dipm": 0.0,
        "dip_num": 1.0,
        "z": 0.0,
        "dielc": 78.09,
        "d_born": 0.0,
        "f_solv": 1.5,
    },
    "Methanol": {
        "MW": 32.04e-3,
        "m": 1.5255,
        "s": 3.2300,
        "e": 188.90,
        "e_assoc": 2899.5,
        "vol_a": 0.03518,
        "assoc_scheme": "2B",
        "dipm": 0.0,
        "dip_num": 1.0,
        "z": 0.0,
        "dielc": 33.05,
        "d_born": 0.0,
        "f_solv": 1.4,
    },
    "Ethanol": {
        "MW": 46.068e-3,
        "m": 2.3827,
        "s": 3.1771,
        "e": 198.24,
        "e_assoc": 2653.4,
        "vol_a": 0.03238,
        "assoc_scheme": "2B",
        "dipm": 0.0,
        "dip_num": 1.0,
        "z": 0.0,
        "dielc": 24.88,
        "d_born": 0.0,
        "f_solv": 1.6,
    },
    "Propanol": {
        "MW": 60.095e-3,
        "m": 3.0,
        "s": 3.2522,
        "e": 233.40,
        "e_assoc": 2276.78,
        "vol_a": 0.01527,
        "assoc_scheme": "2B",
        "dipm": 0.0,
        "dip_num": 1.0,
        "z": 0.0,
        "dielc": 20.47,
        "d_born": 0.0,
        "f_solv": 1.0,
    },
    "Butanol": {
        "MW": 74.1216e-3,
        "m": 2.7510,
        "s": 3.6139,
        "e": 259.59,
        "e_assoc": 2544.56,
        "vol_a": 0.00669,
        "assoc_scheme": "2B",
        "dipm": 0.0,
        "dip_num": 1.0,
        "z": 0.0,
        "dielc": 17.51,
        "d_born": 0.0,
        "f_solv": 1.0,
    },
    "Hexane": {
        "MW": 86.178e-3,
        "m": 3.0578,
        "s": 3.7983,
        "e": 236.77,
        "e_assoc": 0.0,
        "vol_a": 0.0,
        "assoc_scheme": "",
        "dipm": 0.0,
        "dip_num": 1.0,
        "z": 0.0,
        "dielc": 2.0,
        "d_born": 0.0,
        "f_solv": 1.0,
    },
    "Dodecane": {
        "MW": 170.334e-3,
        "m": 5.3060,
        "s": 3.8959,
        "e": 249.21,
        "e_assoc": 0.0,
        "vol_a": 0.0,
        "assoc_scheme": "",
        "dipm": 0.0,
        "dip_num": 1.0,
        "z": 0.0,
        "dielc": 2.0,
        "d_born": 0.0,
        "f_solv": 1.0,
    },
}

_LINEAR_T_RE = re.compile(
    r"^\s*([+\-]?\d*\.?\d+(?:[eE][+\-]?\d+)?(?:\*10\^[+\-]?\d+)?)\s*\*?\s*T(?:\s*/\s*K)?\s*([+\-]\s*\d*\.?\d+(?:[eE][+\-]?\d+)?)\s*$",
    flags=re.IGNORECASE,
)
_FLOAT_RE = re.compile(r"[+\-]?\d*\.?\d+(?:[eE][+\-]?\d+)?")

_RULE_ALIASES = {
    "constant": 0,
    "rule0": 0,
    "linear": 1,
    "linear-mixing-mole": 1,
    "rule1": 1,
    "linear-mixing-weight": 2,
    "rule2": 2,
    "combined": 3,
    "rule3": 3,
    "empirical": 4,
    "rule4": 4,
    "rule5": 5,
    "rule6": 6,
}
_DIFF_MODE_ALIASES = {"analytic": 0, "numeric": 1}
_BORN_RADIUS_ALIASES = {
    "sigma": 1,
    "sigma_const": 2,
    "sigma_temp": 3,
    "dborn": 4,
    "dborn_delta": 5,
}
_BORN_TERM_DEFAULT = {
    "numerical": False,
    "sum_term": True,
    "deps_dx_term": True,
    "d_born_mode": 1,
}
_CANONICAL_ELEC_MODEL = {
    "born_model": 1,
    "dielc_rule": "linear-mixing-mole",
    "dielc_diff_rule": "same",
    "dielc_diff_mode": "analytic",
    "born_term_options": copy.deepcopy(_BORN_TERM_DEFAULT),
    "eps_r_bulk": "mix",
    "bjeruum_treatment": False,
}

_DATASET_CACHE: Dict[str, dict] = {}
_MISSING = object()


def _coerce_bool(value) -> bool:
    if isinstance(value, bool):
        return value
    if isinstance(value, (int, np.integer)):
        return bool(value)
    if isinstance(value, str):
        v = value.strip().lower()
        if v in {"1", "true", "yes", "y", "on"}:
            return True
        if v in {"0", "false", "no", "n", "off"}:
            return False
    raise ValueError(f"Cannot coerce value '{value}' to bool.")


def _deep_update(base: dict, updates: dict) -> dict:
    merged = copy.deepcopy(base)
    for key, value in updates.items():
        if isinstance(value, dict) and isinstance(merged.get(key), dict):
            merged[key] = _deep_update(merged[key], value)
        else:
            merged[key] = value
    return merged


def _normalize_component(name: str) -> str:
    return COMPONENT_ALIASES.get(name, name)


def available_datasets() -> list[str]:
    if not DATASET_ROOT.exists():
        return []
    return sorted(p.name for p in DATASET_ROOT.iterdir() if p.is_dir())


def _dataset_dir(dataset_name: str) -> Path:
    path = DATASET_ROOT / dataset_name
    if not path.exists():
        raise FileNotFoundError(f"Unknown dataset '{dataset_name}'. Available datasets: {available_datasets()}")
    return path


def _read_csv(path: Path) -> list[dict]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def _load_matrix(path: Path) -> dict[tuple[str, str], str]:
    if not path.exists():
        return {}
    rows = _read_csv(path)
    if not rows:
        return {}
    columns = [c for c in rows[0].keys() if c and c != "component"]
    matrix = {}
    for row in rows:
        row_comp = _normalize_component(str(row.get("component", "")).strip())
        if not row_comp:
            continue
        for col in columns:
            col_comp = _normalize_component(col.strip())
            matrix[(row_comp, col_comp)] = str(row.get(col, "") or "").strip()
    return matrix


def _strip_preset_keys(canonical: dict) -> dict:
    cleaned = copy.deepcopy(canonical)
    elec = cleaned.get("elec_model")
    if isinstance(elec, dict):
        elec.pop("preset", None)
        elec.pop("base", None)
    return cleaned


def _load_canonical_user_options(dataset_dir: Path) -> dict:
    path = dataset_dir / "user_options.json"
    if not path.exists():
        return {}
    payload = json.loads(path.read_text(encoding="utf-8"))
    canonical = payload.get("canonical_user_options", {})
    if not isinstance(canonical, dict):
        return {}
    return _strip_preset_keys(canonical)


def _load_dataset(dataset_name: str) -> dict:
    if dataset_name in _DATASET_CACHE:
        return _DATASET_CACHE[dataset_name]

    dataset_dir = _dataset_dir(dataset_name)
    pure_path = dataset_dir / "pure.csv"
    if not pure_path.exists():
        raise FileNotFoundError(f"Missing pure.csv for dataset '{dataset_name}' at {pure_path}")

    pure_rows = _read_csv(pure_path)
    pure_map = {}
    for row in pure_rows:
        comp = _normalize_component(str(row.get("component", "")).strip())
        if not comp:
            continue
        pure_map[comp] = {k: str(v or "").strip() for k, v in row.items()}

    bi_dir = dataset_dir / "binary_interaction"
    data = {
        "dataset_name": dataset_name,
        "dataset_dir": dataset_dir,
        "pure": pure_map,
        "k_ij": _load_matrix(bi_dir / "k_ij.csv"),
        "l_ij": _load_matrix(bi_dir / "l_ij.csv"),
        "k_hb": _load_matrix(bi_dir / "k_hb_ij.csv") or _load_matrix(bi_dir / "khb_ij.csv"),
        "canonical_user_options": _load_canonical_user_options(dataset_dir),
    }
    _DATASET_CACHE[dataset_name] = data
    return data


def _maybe_float(raw) -> float | None:
    if raw is None:
        return None
    if isinstance(raw, (float, int, np.floating, np.integer)):
        value = float(raw)
        return value if math.isfinite(value) else None
    text = str(raw).strip()
    if not text or text.lower() in {"nan", "none", "null"}:
        return None
    try:
        value = float(text)
        return value if math.isfinite(value) else None
    except ValueError:
        return None


def _parse_t_coefficient(token: str) -> float:
    text = token.replace(" ", "")
    if "*10^" in text:
        mantissa, exponent = text.split("*10^", 1)
        return float(mantissa) * (10 ** int(exponent))
    return float(text)


def _parse_linear_t_expression(raw: str, T: float) -> float | None:
    text = raw.replace(" ", "")
    text = text.replace("/K", "").replace("/k", "")
    match = _LINEAR_T_RE.match(text)
    if not match:
        return None
    slope = _parse_t_coefficient(match.group(1))
    intercept = float(match.group(2).replace(" ", ""))
    return slope * T + intercept


def _parse_water_sigma_expression(raw: str, T: float) -> float | None:
    text = raw.replace(" ", "").lower()
    if "2.7927" in text and "10.11" in text and "1.417" in text:
        return 2.7927 + 10.11 * np.exp(-0.01775 * T) - 1.417 * np.exp(-0.01146 * T)
    return None


def _parse_association_scheme(raw: str) -> str | None:
    text = str(raw or "").strip()
    if not text or text.lower() in {"nan", "none", "null"}:
        return None
    return text


def _deterministic_default(component: str, prop: str, T: float):
    entry = _COMPONENT_DEFAULTS.get(component)
    if entry is not None and prop in entry:
        value = entry[prop]
        return value(T) if callable(value) else value

    is_ion = component.endswith("+") or component.endswith("-")
    if is_ion:
        if prop == "z":
            return 1.0 if component.endswith("+") else -1.0
        if prop == "m":
            return 1.0
        if prop in {"e_assoc", "vol_a", "dipm", "d_born"}:
            return 0.0
        if prop == "dip_num":
            return 1.0
        if prop == "assoc_scheme":
            return None
        if prop == "dielc":
            return 8.0
        if prop == "f_solv":
            return 1.0

    if prop == "assoc_scheme":
        return None

    # Generic charge inference for unknown ions.
    if prop == "z":
        if component.endswith("+"):
            return 1.0
        if component.endswith("-"):
            return -1.0
    return _MISSING


def _parse_cell_value(raw, *, dataset: str, component: str, field: str, T: float):
    if field == "assoc_scheme":
        return _parse_association_scheme(raw)

    numeric = _maybe_float(raw)
    if numeric is not None:
        return numeric

    text = str(raw or "").strip()
    if not text or text.lower() in {"nan", "none", "null"}:
        return None

    if field == "s" and component == "H2O":
        parsed = _parse_water_sigma_expression(text, T)
        if parsed is not None:
            return parsed

    linear = _parse_linear_t_expression(text, T)
    if linear is not None:
        return linear

    if "=" in text:
        rhs = text.split("=")[-1].strip()
        rhs_numeric = _maybe_float(rhs)
        if rhs_numeric is not None:
            return rhs_numeric

    numbers = _FLOAT_RE.findall(text)
    if len(numbers) == 1:
        return float(numbers[0])

    raise ValueError(
        f"Unsupported value in dataset '{dataset}', component '{component}', field '{field}': '{text}'."
    )


def _resolve_component_field(dataset: dict, component: str, field: str, T: float):
    row = dataset["pure"].get(component)
    if row is None:
        raise KeyError(f"Component '{component}' is missing in dataset '{dataset['dataset_name']}' pure.csv.")

    parsed = _parse_cell_value(
        row.get(field, ""),
        dataset=dataset["dataset_name"],
        component=component,
        field=field,
        T=T,
    )
    if parsed is not None:
        return parsed

    fallback = _deterministic_default(component, field, T)
    if fallback is not _MISSING:
        return fallback

    raise KeyError(
        f"Missing required value in dataset '{dataset['dataset_name']}', component '{component}', field '{field}'."
    )


def _matrix_value(dataset: dict, matrix_name: str, c1: str, c2: str, T: float) -> float:
    matrix = dataset[matrix_name]
    raw = matrix.get((c1, c2))
    if raw is None:
        raw = matrix.get((c2, c1))
    if raw is None or not str(raw).strip():
        return 0.0
    value = _parse_cell_value(
        raw,
        dataset=dataset["dataset_name"],
        component=f"{c1}|{c2}",
        field=matrix_name,
        T=T,
    )
    if value is None:
        return 0.0
    if isinstance(value, str):
        raise ValueError(
            f"Non-numeric matrix value in dataset '{dataset['dataset_name']}' for pair '{c1}|{c2}' in {matrix_name}."
        )
    return float(value)


def _as_rule_number(value, aliases: dict[str, int]) -> int:
    if isinstance(value, (int, np.integer)):
        return int(value)
    if isinstance(value, str):
        key = value.strip().lower()
        if key.isdigit() or (key.startswith("-") and key[1:].isdigit()):
            return int(key)
        if key in aliases:
            return int(aliases[key])
    raise ValueError(f"Unknown rule option '{value}'. Supported aliases: {sorted(aliases.keys())}.")


def _resolve_born_eps_mode(value) -> int:
    if isinstance(value, (int, np.integer)):
        mode = int(value)
        if mode in (0, 1):
            return mode
    key = str(value).strip().lower()
    if key in {"mix", "bulk", "eps_r_mix"}:
        return 0
    if key in {"solvent", "eps_r_solvent"}:
        return 1
    raise ValueError("born_rel_perm/eps_r_bulk must be 'mix' or 'solvent'.")


def _resolve_born_radius_mode(value) -> int:
    if isinstance(value, (int, np.integer)):
        mode = int(value)
    else:
        key = str(value).strip().lower()
        if key.isdigit() or (key.startswith("-") and key[1:].isdigit()):
            mode = int(key)
        elif key in _BORN_RADIUS_ALIASES:
            mode = int(_BORN_RADIUS_ALIASES[key])
        else:
            raise ValueError(
                f"Unknown born d_born mode '{value}'. Supported values are 1,2,3,4,5 and aliases {sorted(_BORN_RADIUS_ALIASES.keys())}."
            )
    if mode < 1 or mode > 5:
        raise ValueError("born d_born mode must be in 1..5.")
    return mode


def _normalize_elec_model(model) -> dict:
    out = copy.deepcopy(_CANONICAL_ELEC_MODEL)
    if model is None:
        return out
    if not isinstance(model, dict):
        raise TypeError(f"elec_model must be a dict, got {type(model).__name__}.")

    normalized = {}
    for key, value in model.items():
        k = key
        if key == "born_rel_perm":
            k = "eps_r_bulk"
        elif key == "emperical":
            k = "empirical"
        normalized[k] = value

    if "born_term_options" in normalized and not isinstance(normalized["born_term_options"], dict):
        raise TypeError("elec_model['born_term_options'] must be a dict.")

    out = _deep_update(out, normalized)
    out["born_model"] = int(out["born_model"])
    if out["born_model"] not in (0, 1, 2):
        raise ValueError("born_model must be one of 0, 1, 2.")

    term = out["born_term_options"]
    term["numerical"] = _coerce_bool(term["numerical"])
    term["sum_term"] = _coerce_bool(term["sum_term"])
    term["deps_dx_term"] = _coerce_bool(term["deps_dx_term"])
    term["d_born_mode"] = _resolve_born_radius_mode(term["d_born_mode"])

    if out["born_model"] == 1 and term["d_born_mode"] not in (1, 2, 3, 4):
        raise ValueError("born_model=1 requires born_term_options['d_born_mode'] in {1,2,3,4}.")
    return out


def _flatten_model_to_runtime(model: dict) -> dict:
    dielc_rule = _as_rule_number(model["dielc_rule"], _RULE_ALIASES)
    dielc_diff_mode = _as_rule_number(model["dielc_diff_mode"], _DIFF_MODE_ALIASES)

    born_model = int(model["born_model"])
    term = model.get("born_term_options", _BORN_TERM_DEFAULT)
    numerical = _coerce_bool(term.get("numerical", False))
    sum_term = _coerce_bool(term.get("sum_term", True))
    deps_dx_term = _coerce_bool(term.get("deps_dx_term", True))
    d_born_mode = _resolve_born_radius_mode(term.get("d_born_mode", 1))

    if born_model == 0:
        born_diff_mode = 0
        born_radius_model = 1
    elif born_model == 1:
        born_radius_model = d_born_mode
        if numerical:
            born_diff_mode = 1
        elif not deps_dx_term:
            born_diff_mode = 3
        elif not sum_term:
            born_diff_mode = 2
        else:
            born_diff_mode = 0
    elif born_model == 2:
        born_radius_model = 5
        born_diff_mode = 0
    else:
        raise ValueError("born_model must be one of 0, 1, 2.")

    if born_model == 2 and born_radius_model != 5:
        raise ValueError("born_model=2 requires born_radius_model=5.")

    dh_model = 2 if _coerce_bool(model.get("bjeruum_treatment", False)) else 1
    born_eps_mode = _resolve_born_eps_mode(model.get("eps_r_bulk", "mix"))

    return {
        "born_model": int(born_model),
        "born_radius_model": int(born_radius_model),
        "born_diff_mode": int(born_diff_mode),
        "born_eps_mode": int(born_eps_mode),
        "dielc_rule": int(dielc_rule),
        "dielc_diff_mode": int(dielc_diff_mode),
        "DH_model": int(dh_model),
        "bjeruum_treatment": bool(dh_model == 2),
    }


def _resolve_runtime_options(user_options=None) -> dict:
    if user_options is None:
        user_options = {}
    if not isinstance(user_options, dict):
        raise TypeError("user_options must be a dict.")

    rejected_legacy = {"born_diff_model", "born_diff_options", "born_diff_mode"} & set(user_options)
    if rejected_legacy:
        raise KeyError(
            "Legacy Born option key(s) {} are no longer accepted. Use user_options['elec_model']['born_term_options'].".format(
                sorted(rejected_legacy)
            )
        )

    allowed = {
        "elec_model",
        "bjeruum_treatment",
        "dielc_rule",
        "dielc_diff_mode",
        "born_rel_perm",
        "eps_r_bulk",
        "DH_model",
        "debug",
    }
    unknown = set(user_options) - allowed
    if unknown:
        raise KeyError(f"Unknown user_options key(s): {sorted(unknown)}")

    model = _normalize_elec_model(user_options.get("elec_model"))
    if "bjeruum_treatment" in user_options:
        model["bjeruum_treatment"] = _coerce_bool(user_options["bjeruum_treatment"])

    runtime = _flatten_model_to_runtime(model)

    for key in ("dielc_rule", "dielc_diff_mode", "DH_model"):
        if key in user_options:
            runtime[key] = int(user_options[key])

    if "born_rel_perm" in user_options and "eps_r_bulk" in user_options:
        raise ValueError("Use only one Born permittivity selector: born_rel_perm or eps_r_bulk.")
    if "born_rel_perm" in user_options:
        runtime["born_eps_mode"] = _resolve_born_eps_mode(user_options["born_rel_perm"])
    elif "eps_r_bulk" in user_options:
        runtime["born_eps_mode"] = _resolve_born_eps_mode(user_options["eps_r_bulk"])

    if "DH_model" in user_options:
        runtime["bjeruum_treatment"] = int(user_options["DH_model"]) == 2

    runtime["debug"] = bool(user_options.get("debug", False))

    return {
        "runtime": runtime,
        "model": model,
        "preset_key": "dataset_canonical",
        "preset": {},
        "catalog": None,
    }


def _default_species_entry(species_name: str) -> dict:
    comp = _normalize_component(species_name)
    entry = _COMPONENT_DEFAULTS.get(comp)
    if entry is None:
        raise KeyError(f"No default data for species '{species_name}'.")
    resolved = {}
    for key, value in entry.items():
        resolved[key] = value(298.15) if callable(value) else value
    return resolved


def molality_to_molefraction(molality, species=None, solvent=None, basis_mass_kg=1.0):
    """Convert salt molality (mol/kg solvent) to species mole-fraction vector."""
    if species is None:
        raise ValueError("species must be provided.")

    species = list(species)
    molality = float(molality)
    basis_mass_kg = float(basis_mass_kg)

    cations = [sp for sp in species if sp.endswith("+")]
    anions = [sp for sp in species if sp.endswith("-")]
    if len(cations) != 1 or len(anions) != 1:
        cations = [sp for sp in species if _default_species_entry(sp).get("z", 0.0) > 0.0]
        anions = [sp for sp in species if _default_species_entry(sp).get("z", 0.0) < 0.0]
    if len(cations) != 1 or len(anions) != 1:
        raise ValueError("Expected exactly one cation and one anion in species list.")

    cation = cations[0]
    anion = anions[0]

    if solvent is None:
        neutrals = [sp for sp in species if _default_species_entry(sp).get("z", 0.0) == 0.0]
        if len(neutrals) != 1:
            raise ValueError("Expected exactly one neutral solvent species when solvent is not specified.")
        solvent = neutrals[0]
    elif solvent not in species:
        raise ValueError(f"Solvent '{solvent}' is not present in species list.")

    z_cat = float(_default_species_entry(cation).get("z", 0.0))
    z_an = float(_default_species_entry(anion).get("z", 0.0))
    if z_cat <= 0.0 or z_an >= 0.0:
        raise ValueError("Invalid cation/anion charges in species list.")

    z_cat_abs = int(round(abs(z_cat)))
    z_an_abs = int(round(abs(z_an)))
    gcd_z = math.gcd(z_cat_abs, z_an_abs)
    v_cat = z_an_abs // gcd_z
    v_an = z_cat_abs // gcd_z

    mw_solvent = float(_default_species_entry(solvent).get("MW", np.nan))
    if not np.isfinite(mw_solvent) or mw_solvent <= 0.0:
        raise ValueError(f"Invalid MW for solvent '{solvent}'.")

    n_solvent = basis_mass_kg / mw_solvent
    n_cation = molality * basis_mass_kg * v_cat
    n_anion = molality * basis_mass_kg * v_an

    n_totals = {sp: 0.0 for sp in species}
    n_totals[solvent] += n_solvent
    n_totals[cation] += n_cation
    n_totals[anion] += n_anion

    total = sum(n_totals.values())
    if total <= 0.0:
        raise ValueError("Computed total moles is non-positive.")

    return np.array([n_totals[sp] / total for sp in species], dtype=float)


def molefraction_to_molality(x, species):
    """Convert mole fractions to molality for 1:1 salt systems."""
    x_arr = np.asarray(x, dtype=float)
    if x_arr.ndim != 1:
        raise ValueError("x must be a 1D array-like vector.")
    if len(species) != x_arr.size:
        raise ValueError("species and x length mismatch.")

    charges = np.asarray([float(_default_species_entry(sp).get("z", 0.0)) for sp in species], dtype=float)
    cation_idx = [i for i, z in enumerate(charges) if z > 0.0]
    solvent_idx = [i for i, z in enumerate(charges) if abs(z) < 1e-12]
    if not cation_idx or not solvent_idx:
        raise ValueError("Need at least one cation and one solvent component.")

    solvent_i = solvent_idx[-1]
    mw_solvent = float(_default_species_entry(species[solvent_i]).get("MW", np.nan))
    if not np.isfinite(mw_solvent) or mw_solvent <= 0.0:
        raise ValueError("Could not resolve solvent MW.")
    if x_arr[solvent_i] <= 0.0:
        raise ValueError("Solvent mole fraction must be > 0 to compute molality.")

    return float(x_arr[cation_idx[0]] / (x_arr[solvent_i] * mw_solvent))


def get_prop_dict(dataset_name: str, species: Iterable[str], x, T: float, user_options: dict | None = None) -> dict:
    """Build a runtime parameter dictionary from a named dataset."""
    del x  # API compatibility placeholder

    dataset = _load_dataset(dataset_name)
    species = list(species)
    components = [_normalize_component(s) for s in species]

    merged_options = _deep_update(dataset["canonical_user_options"], user_options or {})
    resolved = _resolve_runtime_options(merged_options)
    runtime = resolved["runtime"]

    prop_dic: dict = {}
    for field in BASE_KEYS:
        values = [_resolve_component_field(dataset, comp, field, T) for comp in components]
        if field == "assoc_scheme":
            prop_dic[field] = list(values)
        else:
            prop_dic[field] = np.asarray(values, dtype=float)

    for field in OPTIONAL_KEYS:
        values = []
        for comp in components:
            try:
                values.append(_resolve_component_field(dataset, comp, field, T))
            except KeyError:
                values.append(0.0)
        prop_dic[field] = np.asarray(values, dtype=float)

    n = len(species)
    k_ij = np.zeros((n, n), dtype=float)
    l_ij = np.zeros((n, n), dtype=float)
    k_hb = np.zeros((n, n), dtype=float)
    for i, c1 in enumerate(components):
        for j, c2 in enumerate(components):
            k_ij[i, j] = _matrix_value(dataset, "k_ij", c1, c2, T)
            l_ij[i, j] = _matrix_value(dataset, "l_ij", c1, c2, T)
            k_hb[i, j] = _matrix_value(dataset, "k_hb", c1, c2, T)

    prop_dic["k_ij"] = k_ij
    prop_dic["l_ij"] = l_ij
    prop_dic["k_hb"] = k_hb

    if np.all(np.abs(prop_dic["z"]) < 1e-12):
        prop_dic["z"] = np.array([])

    prop_dic["elec_model"] = copy.deepcopy(resolved["model"])
    prop_dic["elec_model_dataset"] = dataset_name
    prop_dic["bjeruum_treatment"] = bool(runtime["bjeruum_treatment"])
    prop_dic["born_model"] = int(runtime["born_model"])
    prop_dic["born_radius_model"] = int(runtime["born_radius_model"])
    prop_dic["born_diff_mode"] = int(runtime["born_diff_mode"])
    prop_dic["born_eps_mode"] = int(runtime["born_eps_mode"])
    prop_dic["DH_model"] = int(runtime["DH_model"])
    prop_dic["dielc_rule"] = int(runtime["dielc_rule"])
    prop_dic["dielc_diff_mode"] = int(runtime["dielc_diff_mode"])
    prop_dic["debug"] = bool(runtime["debug"])
    return prop_dic
