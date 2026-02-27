"""Build ePC-SAFT electrolyte parameter catalog JSON + CSV reference tables.

This script curates paper tables into a set-aware catalog for:
- component parameters (component -> parameter -> set -> value)
- kij parameter sets
- elec_model presets/aliases
"""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any, Dict, List, Tuple


REPO_ROOT = Path(__file__).resolve().parents[2]
PAPER_DIR = REPO_ROOT / "docs" / "papers" / "md"
OUT_DIR = REPO_ROOT / "data" / "pcsaft_parameters" / "catalog"
OUT_JSON = OUT_DIR / "pcsaft-parameter-catalog.json"


SOURCE_NOTES = {
    "2001": {
        "file": "Gross, Sadowski - 2001 - PC-SAFT An equation of state based on a perturbation theory for chain molec.md",
        "tables": "non-electrolyte baseline only",
        "notes": "No ion parameter set; scanned for provenance completeness.",
    },
    "2002": {
        "file": "Gross, Sadowski - 2002 - Application of the PC-SAFT equation of state to associating systems.md",
        "tables": "Table 1",
        "notes": "Used for associating-solvent parameters (water/methanol/ethanol baseline).",
    },
    "2005": {
        "file": "Cameretti, Sadowski, Mollerup - 2005 - Modeling of Aqueous Electrolyte Solutions with Perturbed-Chai.md",
        "tables": "Table 1, Table 2",
        "notes": "Original ePC-SAFT aqueous electrolyte set; kij=0 convention.",
    },
    "2008": {
        "file": "Held, Cameretti, Sadowski - 2008 - Modeling aqueous electrolyte solutions. Part 1. Fully dissociated.md",
        "tables": "Table 1, Table 2",
        "notes": "Strategy-1 ion parameters and temperature-dependent water sigma.",
    },
    "2014_s1": {
        "file": "Held et al. - 2014 - ePC-SAFT Revised.md",
        "tables": "Table 2 (strategy 1), Table 1 (water)",
        "notes": "Strategy-1 values reused from 2008 family.",
    },
    "2014_s2": {
        "file": "Held et al. - 2014 - ePC-SAFT Revised.md",
        "tables": "Table 2 (strategy 2), Table 3",
        "notes": "Revised ion set + ion-ion and water-ion kij.",
    },
    "2019": {
        "file": "Bülow, Ji, Held - 2019 - Incorporating a concentration-dependent dielectric constant into ePC-SAFT (1).md",
        "tables": "modeling concept",
        "notes": "No new pure-ion set; dielectric concentration dependency focus.",
    },
    "2020": {
        "file": "Bülow, Ascani, Held - 2020 - ePC-SAFT advanced - Part I Physical meaning of including a concentratio.md",
        "tables": "Table 1, Table 2, Table 3, Table 4",
        "notes": "Advanced electrolyte baseline used for transfer to organics.",
    },
    "2021": {
        "file": "Bülow, Ascani, Held - 2021 - ePC-SAFT advanced – Part II Application to Salt Solubility in Ionic and.md",
        "tables": "Table 8",
        "notes": "Bjerrum-related ion-organic kij correlation set at 298.15 K.",
    },
    "2025": {
        "file": "Figiel, Yu, Held - 2025 - Predicting Thermodynamic Properties of Ions in Single Solvents and in Mixe.md",
        "tables": "Table 2, Table 3, Table 4, Table 5",
        "notes": "Modified Born treatment with fk and d_born references.",
    },
}


def _water_sigma_expr() -> Dict[str, Any]:
    return {"type": "water_sigma_2008"}


def _linear_t(a: float, b: float) -> Dict[str, Any]:
    return {"type": "linear_t", "a": a, "b": b}


def _build_catalog() -> Dict[str, Any]:
    component_parameters: Dict[str, Dict[str, Dict[str, Any]]] = {
        "H2O": {
            "MW": {"default": 18.01528e-3},
            "m": {
                "salt_2005": 1.09528,
                "salt_2008": 1.2047,
                "salt_2014_s1": 1.2047,
                "salt_2014_s2": 1.2047,
                "salt_2020": 1.2047,
                "salt_2025": 1.2047,
                "default": 1.2047,
            },
            "s": {
                "salt_2005": 2.88980,
                "salt_2008": _water_sigma_expr(),
                "salt_2014_s1": _water_sigma_expr(),
                "salt_2014_s2": _water_sigma_expr(),
                "salt_2020": _water_sigma_expr(),
                "salt_2025": _water_sigma_expr(),
                "default": _water_sigma_expr(),
            },
            "e": {
                "salt_2005": 365.956,
                "salt_2008": 353.9449,
                "salt_2014_s1": 353.9449,
                "salt_2014_s2": 353.95,
                "salt_2020": 353.95,
                "salt_2025": 353.95,
                "default": 353.9449,
            },
            "e_assoc": {
                "salt_2005": 2515.6706,
                "salt_2008": 2425.6714,
                "salt_2014_s1": 2425.6714,
                "salt_2014_s2": 2425.7,
                "salt_2020": 2425.7,
                "salt_2025": 2425.7,
                "default": 2425.7,
            },
            "vol_a": {
                "salt_2005": 0.0348679836,
                "salt_2008": 0.04509,
                "salt_2014_s1": 0.04509,
                "salt_2014_s2": 0.04509,
                "salt_2020": 0.04509,
                "salt_2025": 0.04509,
                "default": 0.04509,
            },
            "assoc_scheme": {"default": "2B"},
            "dipm": {"default": 0.0},
            "dip_num": {"default": 1},
            "z": {"default": 0.0},
            "dielc": {"default": 78.09},
            "d_born": {"default": 0.0},
            "f_solv": {"default": 1.5, "2025": 1.5},
        },
        "Methanol": {
            "MW": {"default": 32.04e-3},
            "m": {"default": 1.5255},
            "s": {"default": 3.2300},
            "e": {"default": 188.90},
            "e_assoc": {"default": 2899.5},
            "vol_a": {"default": 0.03518},
            "assoc_scheme": {"default": "2B"},
            "dipm": {"default": 0.0},
            "dip_num": {"default": 1},
            "z": {"default": 0.0},
            "dielc": {"default": 33.05},
            "d_born": {"default": 0.0},
            "f_solv": {"default": 1.0, "2025": 1.4},
        },
        "Ethanol": {
            "MW": {"default": 46.068e-3},
            "m": {"default": 2.3827},
            "s": {"default": 3.1771},
            "e": {"default": 198.24},
            "e_assoc": {"default": 2653.4},
            "vol_a": {"default": 0.03238},
            "assoc_scheme": {"default": "2B"},
            "dipm": {"default": 0.0},
            "dip_num": {"default": 1},
            "z": {"default": 0.0},
            "dielc": {"default": 24.88},
            "d_born": {"default": 0.0},
            "f_solv": {"default": 1.0, "2025": 1.6},
        },
    }

    ion_rows = {
        "Li+": {"2005": (1.8059, 1110.9261), "2008": (1.8177, 2697.2795), "2014_s1": (1.8177, 2697.2795), "2014_s2": (2.8449, 360.0), "2020": (2.8449, 360.0), "2025": (2.8449, 360.0), "z": 1.0, "MW": 6.941e-3, "d_born": 2.784},
        "Na+": {"2005": (1.6262, 119.8060), "2008": (2.4122, 646.0504), "2014_s1": (2.4122, 646.0504), "2014_s2": (2.8232, 230.0), "2020": (2.8232, 230.0), "2025": (2.8232, 230.0), "z": 1.0, "MW": 22.98e-3, "d_born": 3.445},
        "K+": {"2005": (2.7602, 8.8773), "2008": (2.9698, 271.0518), "2014_s1": (2.9698, 271.0518), "2014_s2": (3.3417, 200.0), "2020": (3.3417, 200.0), "2025": (3.3417, 200.0), "z": 1.0, "MW": 39.0983e-3, "d_born": 4.150},
        "Cl-": {"2005": (3.5991, 359.6604), "2008": (3.0575, 47.2878), "2014_s1": (3.0575, 47.2878), "2014_s2": (2.7560, 170.0), "2020": (2.7560, 170.0), "2025": (2.7560, 170.0), "z": -1.0, "MW": 35.453e-3, "d_born": 4.100},
        "Br-": {"2005": (3.8225, 524.0636), "2008": (3.4573, 60.2216), "2014_s1": (3.4573, 60.2216), "2014_s2": (3.0707, 190.0), "2020": (3.0707, 190.0), "2025": (3.0707, 190.0), "z": -1.0, "MW": 79.904e-3, "d_born": 4.480},
        "I-": {"2005": (4.1766, 413.0494), "2008": (3.9319, 80.4347), "2014_s1": (3.9319, 80.4347), "2014_s2": (3.6672, 200.0), "2020": (3.6672, 200.0), "2025": (3.6672, 200.0), "z": -1.0, "MW": 126.90447e-3, "d_born": 4.985},
    }

    for ion, row in ion_rows.items():
        component_parameters[ion] = {
            "MW": {"default": row["MW"]},
            "m": {"default": 1.0},
            "s": {k: row[k][0] for k in ("2005", "2008", "2014_s1", "2014_s2", "2020", "2025")} | {"default": row["2020"][0]},
            "e": {k: row[k][1] for k in ("2005", "2008", "2014_s1", "2014_s2", "2020", "2025")} | {"default": row["2020"][1]},
            "e_assoc": {"default": 0.0},
            "vol_a": {"default": 0.0},
            "assoc_scheme": {"default": None},
            "dipm": {"default": 0.0},
            "dip_num": {"default": 1},
            "z": {"default": row["z"]},
            "dielc": {"default": 8.0},
            "d_born": {"default": row["d_born"], "2025": row["d_born"]},
            "f_solv": {"default": 1.0},
        }

    kij_ion_pairs_2014 = {
        "Li+|Cl-": 0.669, "Li+|Br-": 0.591, "Li+|I-": 0.002,
        "Na+|Cl-": 0.317, "Na+|Br-": 0.290, "Na+|I-": 0.018,
        "K+|Cl-": 0.064, "K+|Br-": -0.102, "K+|I-": -0.312,
    }

    kij_parameters = {
        "2005": {"default_zero": True, "ion_pair_default": 0.0, "pairs": {}},
        "2008": {"default_zero": True, "ion_pair_default": 1.0, "pairs": {}},
        "2014_s1": {"default_zero": True, "ion_pair_default": 1.0, "pairs": {}},
        "2014_s2": {
            "default_zero": True,
            "pairs": {
                **kij_ion_pairs_2014,
                "Na+|H2O": _linear_t(-0.007981, 2.37999),
                "K+|H2O": _linear_t(-0.004012, 1.3959),
                "Li+|H2O": -0.25,
                "Cl-|H2O": -0.25,
                "Br-|H2O": -0.25,
                "I-|H2O": -0.25,
            },
        },
        "2020": {
            "default_zero": True,
            "pairs": {
                **kij_ion_pairs_2014,
                "Na+|H2O": 0.0045,
                "K+|H2O": 0.1997,
                "Li+|H2O": -0.25,
                "Cl-|H2O": -0.25,
                "Br-|H2O": -0.25,
                "I-|H2O": -0.25,
            },
        },
        "2021": {
            "default_zero": True,
            "pairs": {
                "Na+|Methanol": -0.31,
                "K+|Methanol": 0.47,
                "Cl-|Methanol": -0.21,
                "Br-|Methanol": -0.42,
                "Na+|Ethanol": 0.42,
                "K+|Ethanol": -0.2,
                "Cl-|Ethanol": -0.15,
                "Br-|Ethanol": -0.35,
                "I-|Ethanol": -0.38,
            },
        },
        "2025": {
            "default_zero": True,
            "pairs": {
                "H+|Cl-": -0.9, "H+|Br-": -0.7, "Li+|Cl-": 0.8, "Li+|Br-": 0.5,
                "Na+|Cl-": 0.8, "Na+|Br-": 0.65, "Na+|I-": 0.45,
                "K+|Cl-": 0.0, "K+|Br-": -0.35, "K+|I-": 0.0,
                "Li+|H2O": -0.4, "Na+|H2O": -0.3, "K+|H2O": -0.1,
                "Cl-|H2O": -0.3, "Br-|H2O": -0.3, "I-|H2O": -0.05,
                "H+|Methanol": -0.3, "Li+|Methanol": -0.9, "Na+|Methanol": -0.25, "K+|Methanol": 0.32,
                "Cl-|Methanol": 0.5, "Br-|Methanol": 0.15, "I-|Methanol": 0.37,
                "H+|Ethanol": -0.6, "Li+|Ethanol": -0.8, "Na+|Ethanol": 0.05, "K+|Ethanol": 0.53,
                "Cl-|Ethanol": 0.8, "Br-|Ethanol": 0.0, "I-|Ethanol": 0.18,
                "H2O|Methanol": -0.0878, "H2O|Ethanol": -0.0617,
            },
        },
    }

    presets = {
        "legacy_default": {
            "parameter_set_key": "default",
            "component_set_key": {"default": "default", "H2O": "salt_2020"},
            "kij_set": "2025",
            "model": {
                "born_model": 1,
                "dielc_rule": "linear-mixing-mole", "dielc_diff_rule": "same",
                "dielc_diff_mode": "analytic",
                "born_term_options": {"numerical": False, "sum_term": True, "deps_dx_term": True, "d_born_mode": 1},
                "eps_r_bulk": "mix", "bjeruum_treatment": False,
            },
        },
        "2005": {
            "parameter_set_key": "2005",
            "component_set_key": {"default": "2005", "H2O": "salt_2005", "Methanol": "default", "Ethanol": "default"},
            "kij_set": "2005",
            "model": {
                "born_model": 0,
                "dielc_rule": "constant", "dielc_diff_rule": "same",
                "dielc_diff_mode": "analytic",
                "born_term_options": {"numerical": False, "sum_term": True, "deps_dx_term": True, "d_born_mode": 1},
                "eps_r_bulk": "mix", "bjeruum_treatment": False,
            },
        },
        "2008": {
            "parameter_set_key": "2008",
            "component_set_key": {"default": "2008", "H2O": "salt_2008"},
            "kij_set": "2008",
            "model": {
                "born_model": 0,
                "dielc_rule": "constant", "dielc_diff_rule": "same",
                "dielc_diff_mode": "analytic",
                "born_term_options": {"numerical": False, "sum_term": True, "deps_dx_term": True, "d_born_mode": 1},
                "eps_r_bulk": "mix", "bjeruum_treatment": False,
            },
        },
        "2014_s1": {
            "parameter_set_key": "2014_s1",
            "component_set_key": {"default": "2014_s1", "H2O": "salt_2014_s1"},
            "kij_set": "2014_s1",
            "model": {
                "born_model": 0,
                "dielc_rule": "constant", "dielc_diff_rule": "same",
                "dielc_diff_mode": "analytic",
                "born_term_options": {"numerical": False, "sum_term": True, "deps_dx_term": True, "d_born_mode": 1},
                "eps_r_bulk": "mix", "bjeruum_treatment": False,
            },
        },
        "2014_s2": {
            "parameter_set_key": "2014_s2",
            "component_set_key": {"default": "2014_s2", "H2O": "salt_2014_s2"},
            "kij_set": "2014_s2",
            "model": {
                "born_model": 0,
                "dielc_rule": "constant", "dielc_diff_rule": "same",
                "dielc_diff_mode": "analytic",
                "born_term_options": {"numerical": False, "sum_term": True, "deps_dx_term": True, "d_born_mode": 1},
                "eps_r_bulk": "mix", "bjeruum_treatment": False,
            },
        },
        "2020": {
            "parameter_set_key": "2020",
            "component_set_key": {"default": "2020", "H2O": "salt_2020"},
            "kij_set": "2020",
            "model": {
                "born_model": 1,
                "dielc_rule": "linear-mixing-mole", "dielc_diff_rule": "same",
                "dielc_diff_mode": "analytic",
                "born_term_options": {"numerical": False, "sum_term": True, "deps_dx_term": True, "d_born_mode": 1},
                "eps_r_bulk": "mix", "bjeruum_treatment": False,
            },
        },
        "2025": {
            "parameter_set_key": "2025",
            "component_set_key": {"default": "2025", "H2O": "salt_2025"},
            "kij_set": "2025",
            "model": {
                "born_model": 2,
                "dielc_rule": "empirical", "dielc_diff_rule": "same",
                "dielc_diff_mode": "analytic",
                "born_term_options": {"numerical": False, "sum_term": True, "deps_dx_term": True, "d_born_mode": 1},
                "eps_r_bulk": "mix", "bjeruum_treatment": False,
            },
        },
    }

    aliases = {
        "preset_int": {"1": "2005", "2": "2008", "3": "2014_s1", "4": "2014_s2", "5": "2020", "6": "2025"},
        "preset_str": {"2005": "2005", "2008": "2008", "2014_s1": "2014_s1", "2014_s2": "2014_s2", "2020": "2020", "2025": "2025"},
        "component_aliases": {
            "H2O-2B-Li": "H2O",
            "H2O-2B-NaCl": "H2O",
            "H2O-Salt-2001": "H2O",
            "H2O-Salt-": "H2O",
            "H2O": "H2O",
        },
    }

    return {
        "metadata": {"schema_version": 1, "description": "ePC-SAFT electrolyte parameter set catalog"},
        "component_parameters": component_parameters,
        "kij_parameters": kij_parameters,
        "elec_model_presets": presets,
        "aliases": aliases,
        "provenance": SOURCE_NOTES,
    }


def _eval_csv_rows(catalog: Dict[str, Any]) -> Dict[str, List[Dict[str, Any]]]:
    rows_by_set: Dict[str, List[Dict[str, Any]]] = {k: [] for k in SOURCE_NOTES}
    for comp, params in catalog["component_parameters"].items():
        for pname, set_values in params.items():
            for set_key, value in set_values.items():
                if set_key == "default":
                    continue
                if set_key.startswith("salt_"):
                    ref_set = set_key.replace("salt_", "")
                else:
                    ref_set = set_key
                if ref_set not in rows_by_set:
                    continue
                rows_by_set[ref_set].append(
                    {
                        "entry_type": "component",
                        "name": comp,
                        "parameter": pname,
                        "set_key": set_key,
                        "value": json.dumps(value, ensure_ascii=True) if isinstance(value, dict) else value,
                    }
                )
    for set_key, cfg in catalog["kij_parameters"].items():
        if set_key not in rows_by_set:
            continue
        for pair, value in cfg.get("pairs", {}).items():
            rows_by_set[set_key].append(
                {
                    "entry_type": "kij_pair",
                    "name": pair,
                    "parameter": "k_ij",
                    "set_key": set_key,
                    "value": json.dumps(value, ensure_ascii=True) if isinstance(value, dict) else value,
                }
            )
    return rows_by_set


def build() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    catalog = _build_catalog()

    OUT_JSON.write_text(json.dumps(catalog, indent=2, sort_keys=True), encoding="utf-8")

    rows_by_set = _eval_csv_rows(catalog)
    for set_key, rows in rows_by_set.items():
        out_csv = OUT_DIR / f"{set_key}.csv"
        with out_csv.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=["entry_type", "name", "parameter", "set_key", "value"])
            writer.writeheader()
            for row in sorted(rows, key=lambda r: (r["entry_type"], r["name"], r["parameter"], r["set_key"])):
                writer.writerow(row)

    index_csv = OUT_DIR / "index.csv"
    with index_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["set_key", "source_file", "source_exists", "tables", "notes"])
        writer.writeheader()
        for set_key, meta in SOURCE_NOTES.items():
            file_path = PAPER_DIR / meta["file"]
            writer.writerow(
                {
                    "set_key": set_key,
                    "source_file": str(file_path),
                    "source_exists": int(file_path.exists()),
                    "tables": meta["tables"],
                    "notes": meta["notes"],
                }
            )

    print(f"Wrote catalog: {OUT_JSON}")
    print(f"Wrote reference CSVs to: {OUT_DIR}")


if __name__ == "__main__":
    build()
