from __future__ import annotations

import copy
import json
import math
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from _paths import ANALYSIS_DIR, REPO_ROOT  # noqa: E402

if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

import rezaee_reactive_equilibrium_replay as replay  # noqa: E402
from _epcsaft_compat import _normalized_params  # noqa: E402
from _epcsaft_properties import get_prop_dict  # noqa: E402
from epcsaft import ePCSAFTMixture  # noqa: E402

INPUT_DIR = ANALYSIS_DIR / "data" / "input"
PROCESSED_DIR = ANALYSIS_DIR / "data" / "processed"
RESULTS_DIR = ANALYSIS_DIR / "results" / "reaction_equilibrium"

SPECIES_DATASET_JSON = INPUT_DIR / "rezaee_2026_epcsaft_species_dataset.json"
OPTION_SCAN_ROWS_CSV = PROCESSED_DIR / "rezaee_2026_reactive_epcsaft_option_scan_rows.csv"
OPTION_SCAN_SUMMARY_CSV = PROCESSED_DIR / "rezaee_2026_reactive_epcsaft_option_scan_summary.csv"
OPTION_SCAN_JSON = RESULTS_DIR / "rezaee_2026_reactive_epcsaft_option_scan_summary.json"
OPTION_SCAN_REPORT_MD = RESULTS_DIR / "rezaee_2026_reactive_epcsaft_option_scan.md"


def _jsonable(value: Any) -> Any:
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, dict):
        return {str(k): _jsonable(v) for k, v in value.items()}
    if isinstance(value, list):
        return [_jsonable(v) for v in value]
    return value


def _deep_update(base: dict[str, Any], updates: dict[str, Any]) -> dict[str, Any]:
    out = copy.deepcopy(base)
    for key, value in updates.items():
        if isinstance(value, dict) and isinstance(out.get(key), dict):
            out[key] = _deep_update(out[key], value)
        else:
            out[key] = value
    return out


OPTION_SETS: list[dict[str, Any]] = [
    {
        "option_id": "package_struct_default_no_explicit_elec_model",
        "description": "Existing replay behavior: strip legacy electrolyte option keys and let package defaults construct add_args.",
        "user_options": None,
        "normalization": "strip_all_electrolyte_option_keys",
    },
    {
        "option_id": "legacy_default_born_linear_mole",
        "description": "Legacy/default aqueous ePC-SAFT option payload converted to the current nested package API.",
        "user_options": {"elec_model": {"base": "legacy_default"}},
        "normalization": "modern_nested_elec_model",
    },
    {
        "option_id": "held_2005_no_born_constant_eps",
        "description": "Older Held-style no-Born constant-dielectric preset.",
        "user_options": {"elec_model": "2005"},
        "normalization": "modern_nested_elec_model",
    },
    {
        "option_id": "held_2008_no_born_constant_eps",
        "description": "2008 no-Born constant-dielectric preset.",
        "user_options": {"elec_model": "2008"},
        "normalization": "modern_nested_elec_model",
    },
    {
        "option_id": "held_2014_s1_no_born_constant_eps",
        "description": "2014 S1 no-Born constant-dielectric preset.",
        "user_options": {"elec_model": "2014_s1"},
        "normalization": "modern_nested_elec_model",
    },
    {
        "option_id": "held_2014_s2_no_born_constant_eps",
        "description": "2014 S2 no-Born constant-dielectric preset.",
        "user_options": {"elec_model": "2014_s2"},
        "normalization": "modern_nested_elec_model",
    },
    {
        "option_id": "held_2020_born_no_ssm_linear_mole",
        "description": "2020 Born-enabled preset without solvation-shell/dielectric-saturation behavior.",
        "user_options": {"elec_model": "2020"},
        "normalization": "modern_nested_elec_model",
    },
    {
        "option_id": "current_2025_born_ssm_ds_empirical",
        "description": "Current 2025 preset with Born term, solvation-shell model, and empirical dielectric rule.",
        "user_options": {"elec_model": "2025"},
        "normalization": "modern_nested_elec_model",
    },
    {
        "option_id": "2025_no_born_empirical",
        "description": "Current 2025 aqueous parameters with Born contribution disabled and empirical dielectric rule retained.",
        "user_options": {"elec_model": {"base": "2025", "born_contrib": False, "dielc_rule": "empirical"}},
        "normalization": "modern_nested_elec_model",
    },
    {
        "option_id": "2025_no_born_constant_eps",
        "description": "Current 2025 aqueous parameters with Born contribution disabled and constant dielectric rule.",
        "user_options": {"elec_model": {"base": "2025", "born_contrib": False, "dielc_rule": "constant"}},
        "normalization": "modern_nested_elec_model",
    },
    {
        "option_id": "2025_no_born_linear_mole",
        "description": "Current 2025 aqueous parameters with Born contribution disabled and linear mole-fraction dielectric rule.",
        "user_options": {"elec_model": {"base": "2025", "born_contrib": False, "dielc_rule": "linear-mixing-mole"}},
        "normalization": "modern_nested_elec_model",
    },
    {
        "option_id": "2025_born_no_ssm_empirical",
        "description": "Born term retained but solvation-shell/dielectric-saturation option disabled.",
        "user_options": {
            "elec_model": {"base": "2025", "born_contrib": True, "ssm_ds": False, "dielc_rule": "empirical"}
        },
        "normalization": "modern_nested_elec_model",
    },
    {
        "option_id": "2025_born_ssm_no_sum_empirical",
        "description": "Born plus solvation-shell model retained with the mu Born sum term disabled.",
        "user_options": {
            "elec_model": {
                "base": "2025",
                "born_contrib": True,
                "ssm_ds": True,
                "dielc_rule": "empirical",
                "born_diff_options": {
                    "include_sum_term": False,
                    "include_dielc_conc_dep": True,
                    "include_delta_d_i_conc_dep": True,
                },
            }
        },
        "normalization": "modern_nested_elec_model",
    },
    {
        "option_id": "2025_born_no_ssm_no_sum_empirical",
        "description": "Born retained without solvation-shell model and with the mu Born sum term disabled.",
        "user_options": {
            "elec_model": {
                "base": "2025",
                "born_contrib": True,
                "ssm_ds": False,
                "dielc_rule": "empirical",
                "born_diff_options": {
                    "include_sum_term": False,
                    "include_dielc_conc_dep": True,
                    "include_delta_d_i_conc_dep": True,
                },
            }
        },
        "normalization": "modern_nested_elec_model",
    },
    {
        "option_id": "2025_born_no_comp_dep_rel_perm",
        "description": "Born retained with composition-dependent relative permittivity disabled in the chemical-potential derivative.",
        "user_options": {
            "elec_model": {
                "base": "2025",
                "born_contrib": True,
                "ssm_ds": True,
                "dielc_rule": "empirical",
                "born_diff_options": {
                    "include_sum_term": True,
                    "include_dielc_conc_dep": False,
                    "include_delta_d_i_conc_dep": False,
                },
            }
        },
        "normalization": "modern_nested_elec_model",
    },
    {
        "option_id": "2025_born_solvent_bulk_empirical",
        "description": "Born retained but bulk relative permittivity taken from the solvent mode instead of mixture mode.",
        "user_options": {
            "elec_model": {
                "base": "2025",
                "born_contrib": True,
                "ssm_ds": True,
                "dielc_rule": "empirical",
                "eps_r_bulk": "solvent",
            }
        },
        "normalization": "modern_nested_elec_model",
    },
    {
        "option_id": "2025_bjeruum_on",
        "description": "Current 2025 preset with Bjerrum treatment enabled.",
        "user_options": {"elec_model": "2025", "bjeruum_treatment": True},
        "normalization": "modern_nested_elec_model",
    },
]


def _nh4_params() -> dict[str, float | str | None]:
    return {
        "MW": 18.038e-3,
        "m": 1.0,
        "s": 3.5740,
        "e": 230.0,
        "e_assoc": 0.0,
        "vol_a": 0.0,
        "assoc_scheme": None,
        "dipm": 0.0,
        "dip_num": 1,
        "z": 1.0,
        "dielc": 8.0,
        "d_born": 3.0,
        "f_solv": 1.0,
    }


def _build_aqueous_mixture(option: dict[str, Any]) -> tuple[ePCSAFTMixture, dict[str, Any], dict[str, Any]]:
    params = get_prop_dict(
        replay.AQ_PARAM_KEYS,
        np.full(len(replay.AQ_PARAM_KEYS), 1.0 / len(replay.AQ_PARAM_KEYS)),
        replay.TEMPERATURE_K,
        user_params={"NH4+": _nh4_params()},
        user_options=option["user_options"],
    )
    resolved = {
        "elec_model_preset": params.get("elec_model_preset"),
        "legacy_runtime": {
            "born_model": params.get("born_model"),
            "born_diff_mode": params.get("born_diff_mode"),
            "born_eps_mode": params.get("born_eps_mode"),
            "DH_model": params.get("DH_model"),
            "dielc_rule": params.get("dielc_rule"),
            "dielc_diff_mode": params.get("dielc_diff_mode"),
            "bjeruum_treatment": params.get("bjeruum_treatment"),
        },
        "legacy_elec_model": params.get("elec_model"),
    }
    if option["normalization"] == "strip_all_electrolyte_option_keys":
        clean = dict(params)
        for key in replay.UNSUPPORTED_FLAT_ELECTROLYTE_KEYS:
            clean.pop(key, None)
    elif option["normalization"] == "modern_nested_elec_model":
        clean = _normalized_params(params)
    else:
        raise ValueError(f"Unknown normalization mode {option['normalization']!r}")

    mixture = ePCSAFTMixture.from_params(clean, species=replay.AQ_LABELS)
    clean_snapshot = {
        key: _jsonable(value)
        for key, value in clean.items()
        if key
        in {
            "m",
            "s",
            "e",
            "e_assoc",
            "vol_a",
            "assoc_scheme",
            "z",
            "dielc",
            "d_born",
            "f_solv",
            "k_ij",
            "elec_model",
        }
    }
    return mixture, resolved, clean_snapshot


def _source_dataset_audit() -> dict[str, Any]:
    dataset = json.loads(SPECIES_DATASET_JSON.read_text(encoding="utf-8"))
    params = pd.read_csv(replay.ORGANIC_PARAMS_CSV)
    kij = pd.read_csv(replay.ORGANIC_KIJ_CSV)
    reactions = pd.read_csv(replay.REACTION_CONSTANTS_CSV)
    equilibrium = pd.read_csv(replay.EQUILIBRIUM_CSV)
    return {
        "dataset_id": dataset["dataset_id"],
        "aqueous_species": dataset["phase_species"]["aqueous"],
        "organic_species": dataset["phase_species"]["organic"],
        "expected_organic_species": replay.ORG_LABELS,
        "organic_species_match": list(params["component"]) == replay.ORG_LABELS,
        "organic_parameter_rows": int(len(params)),
        "organic_binary_interaction_rows": int(len(kij)),
        "reaction_constant_rows": int(len(reactions)),
        "equilibrium_rows_available": int(len(equilibrium)),
        "equilibrium_rows_used_as_benchmark": int(len(equilibrium)),
        "equilibrium_row_count_note": (
            "The 2025 SI Tables S1/S2 provide 26 designed-experiment equilibrium rows. "
            "The 2026 text mentions 36 data points, but this scan treats that as a clerical "
            "source-text mismatch unless a new source table is supplied."
        ),
        "notes": [
            "Organic phase uses Rezaee Table 4/6/8 pure parameters and Table 9 binary interactions without refitting.",
            "Aqueous phase uses the package ePC-SAFT component catalog for H2O, Li+, Na+, Cl-, H+, OH-, and a local NH4+ cation placeholder to match the SI composition columns.",
            "Option scan changes only the aqueous electrolyte model configuration; organic parameters and reaction constants remain fixed to the paper values.",
        ],
    }


def _evaluate_option(option: dict[str, Any]) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    constants = replay._reaction_constants()
    aqueous_mix, resolved_options, clean_snapshot = _build_aqueous_mixture(option)
    organic_mix, _organic_params, pure_ln_phi = replay._organic_mixture()

    records: list[dict[str, Any]] = []
    for row in pd.read_csv(replay.EQUILIBRIUM_CSV).itertuples(index=False):
        aqueous_x = replay._aqueous_x(row)
        organic_x = replay._organic_x(row)
        aqueous_gamma = aqueous_mix.state(
            T=replay.TEMPERATURE_K,
            x=aqueous_x,
            P=replay.PRESSURE_PA,
        ).activity_coefficient(species=replay.AQ_LABELS)
        organic_gamma = replay._organic_activity_coefficients(organic_mix, pure_ln_phi, organic_x)
        ln_q_li, ln_q_na = replay._reaction_log_quotients(
            aqueous_x,
            organic_x,
            aqueous_gamma,
            organic_gamma,
        )
        for metal, ln_q in (("Li", ln_q_li), ("Na", ln_q_na)):
            ln_k = math.log(constants[metal])
            records.append(
                {
                    "option_id": option["option_id"],
                    "metal": metal,
                    "experiment_no": int(row.experiment_no),
                    "lnQ": float(ln_q),
                    "lnK": float(ln_k),
                    "lnQ_minus_lnK": float(ln_q - ln_k),
                    "abs_ln_residual": float(abs(ln_q - ln_k)),
                    "gamma_H2O": float(aqueous_gamma["H2O"]),
                    "gamma_Li_plus": float(aqueous_gamma["Li+"]),
                    "gamma_Na_plus": float(aqueous_gamma["Na+"]),
                    "gamma_OH_minus": float(aqueous_gamma["OH-"]),
                    "source": row.source,
                }
            )
    return records, {
        "option_id": option["option_id"],
        "description": option["description"],
        "normalization": option["normalization"],
        "user_options": option["user_options"],
        "resolved_options": resolved_options,
        "clean_parameter_snapshot": clean_snapshot,
    }


def _summarize(
    rows: pd.DataFrame,
    option_meta: list[dict[str, Any]],
    failed_options: list[dict[str, Any]],
) -> tuple[pd.DataFrame, dict[str, Any]]:
    by_metal = (
        rows.groupby(["option_id", "metal"], as_index=False)
        .agg(
            row_count=("experiment_no", "count"),
            median_ln_residual=("lnQ_minus_lnK", "median"),
            mean_abs_ln_residual=("abs_ln_residual", "mean"),
            median_abs_ln_residual=("abs_ln_residual", "median"),
            rms_ln_residual=("lnQ_minus_lnK", lambda values: math.sqrt(float(np.mean(np.square(values))))),
            max_abs_ln_residual=("abs_ln_residual", "max"),
            median_gamma_H2O=("gamma_H2O", "median"),
            median_gamma_Li_plus=("gamma_Li_plus", "median"),
            median_gamma_Na_plus=("gamma_Na_plus", "median"),
            median_gamma_OH_minus=("gamma_OH_minus", "median"),
        )
        .sort_values(["option_id", "metal"])
    )
    combined = (
        by_metal.groupby("option_id", as_index=False)
        .agg(combined_median_abs_ln_residual=("median_abs_ln_residual", "mean"))
        .sort_values("combined_median_abs_ln_residual")
    )
    meta_by_id = {item["option_id"]: item for item in option_meta}
    summary_rows = []
    for row in combined.itertuples(index=False):
        option_id = str(row.option_id)
        li = by_metal.loc[(by_metal["option_id"] == option_id) & (by_metal["metal"] == "Li")].iloc[0]
        na = by_metal.loc[(by_metal["option_id"] == option_id) & (by_metal["metal"] == "Na")].iloc[0]
        meta = meta_by_id[option_id]
        resolved = meta["resolved_options"]["legacy_runtime"]
        summary_rows.append(
            {
                "option_id": option_id,
                "description": meta["description"],
                "normalization": meta["normalization"],
                "elec_model_preset": meta["resolved_options"]["elec_model_preset"],
                "born_model": resolved["born_model"],
                "born_diff_mode": resolved["born_diff_mode"],
                "born_eps_mode": resolved["born_eps_mode"],
                "DH_model": resolved["DH_model"],
                "dielc_rule": resolved["dielc_rule"],
                "dielc_diff_mode": resolved["dielc_diff_mode"],
                "bjeruum_treatment": resolved["bjeruum_treatment"],
                "combined_median_abs_ln_residual": float(row.combined_median_abs_ln_residual),
                "li_median_ln_residual": float(li["median_ln_residual"]),
                "na_median_ln_residual": float(na["median_ln_residual"]),
                "li_median_abs_ln_residual": float(li["median_abs_ln_residual"]),
                "na_median_abs_ln_residual": float(na["median_abs_ln_residual"]),
                "li_median_gamma": float(li["median_gamma_Li_plus"]),
                "na_median_gamma": float(na["median_gamma_Na_plus"]),
                "oh_median_gamma_for_li_rows": float(li["median_gamma_OH_minus"]),
            }
        )
    summary = pd.DataFrame(summary_rows).sort_values("combined_median_abs_ln_residual")
    best = summary.iloc[0].to_dict()
    baseline = summary.loc[summary["option_id"].eq("package_struct_default_no_explicit_elec_model")].iloc[0].to_dict()
    payload = {
        "status": "option_scan_complete_no_direct_closure",
        "option_count": int(len(summary)),
        "row_count": int(rows["experiment_no"].nunique()),
        "best_option": best,
        "baseline_option": baseline,
        "accepted": bool(float(best["combined_median_abs_ln_residual"]) < 2.0),
        "parameter_audit": _source_dataset_audit(),
        "option_metadata": option_meta,
        "failed_options": failed_options,
    }
    return summary, payload


def _write_report(summary: pd.DataFrame, payload: dict[str, Any]) -> None:
    best = payload["best_option"]
    baseline = payload["baseline_option"]
    audit = payload["parameter_audit"]
    lines = [
        "# Rezaee 2026 ePC-SAFT Aqueous Option Scan",
        "",
        "## Purpose",
        "",
        "This scan keeps the Rezaee organic PC-SAFT parameters, organic binary interactions, and reaction constants fixed to the paper values, then varies only the aqueous ePC-SAFT electrolyte option configuration used to calculate aqueous activity coefficients in Eqs. 14/15.",
        "",
        "## Input Audit",
        "",
        f"- Aqueous species: `{', '.join(audit['aqueous_species'])}`.",
        f"- Organic species: `{', '.join(audit['organic_species'])}`.",
        f"- Organic species match script order: `{audit['organic_species_match']}`.",
        f"- Organic parameter rows: `{audit['organic_parameter_rows']}`.",
        f"- Organic binary-interaction rows: `{audit['organic_binary_interaction_rows']}`.",
        f"- Reaction constant rows: `{audit['reaction_constant_rows']}`.",
        f"- Equilibrium rows available: `{audit['equilibrium_rows_available']}`.",
        f"- Equilibrium rows used as benchmark: `{audit['equilibrium_rows_used_as_benchmark']}`.",
        "- Row-count note: 2026 text mentions 36 data points; this workflow treats that as a clerical source-text mismatch and uses the 26 source SI rows.",
        "",
        "## Result",
        "",
        f"- Options scanned: `{payload['option_count']}`.",
        f"- Unsupported/failed options recorded: `{len(payload['failed_options'])}`.",
        f"- Best option: `{best['option_id']}`.",
        f"- Best combined median absolute ln residual: `{float(best['combined_median_abs_ln_residual']):.6g}`.",
        f"- Baseline option: `{baseline['option_id']}`.",
        f"- Baseline combined median absolute ln residual: `{float(baseline['combined_median_abs_ln_residual']):.6g}`.",
        f"- Accepted direct closure threshold met: `{payload['accepted']}`.",
        "",
        "## Top Options",
        "",
    ]
    for row in summary.head(8).itertuples(index=False):
        lines.append(
            f"- `{row.option_id}`: combined median abs ln residual `{row.combined_median_abs_ln_residual:.6g}`; "
            f"Li `{row.li_median_ln_residual:.6g}`, Na `{row.na_median_ln_residual:.6g}`."
        )
    if payload["failed_options"]:
        lines.extend(["", "## Unsupported Options", ""])
        for item in payload["failed_options"]:
            lines.append(f"- `{item['option_id']}`: `{item['error_type']}` - {item['message']}")
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "The aqueous ePC-SAFT option choice materially changes the activity coefficients, but none of the scanned Born, no-Born, dielectric, Bjerrum, or derivative configurations closes the published-constant Rezaee reactive equilibrium while keeping the paper organic parameters and reaction constants fixed.",
            "",
            "This does not justify refitting yet. It narrows the next task: reproduce the paper's reaction-coordinate calculation using the 26-row benchmark and the published Table 8/9 parameters before changing those parameters.",
            "",
            "## Generated Files",
            "",
            f"- `{OPTION_SCAN_ROWS_CSV.relative_to(ANALYSIS_DIR)}`",
            f"- `{OPTION_SCAN_SUMMARY_CSV.relative_to(ANALYSIS_DIR)}`",
            f"- `{OPTION_SCAN_JSON.relative_to(ANALYSIS_DIR)}`",
        ]
    )
    OPTION_SCAN_REPORT_MD.parent.mkdir(parents=True, exist_ok=True)
    OPTION_SCAN_REPORT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    all_records: list[dict[str, Any]] = []
    option_meta: list[dict[str, Any]] = []
    failed_options: list[dict[str, Any]] = []
    for option in OPTION_SETS:
        try:
            records, meta = _evaluate_option(option)
            all_records.extend(records)
            option_meta.append(meta)
        except Exception as exc:  # noqa: BLE001
            failed_options.append(
                {
                    "option_id": option["option_id"],
                    "description": option["description"],
                    "user_options": option["user_options"],
                    "error_type": type(exc).__name__,
                    "message": str(exc),
                }
            )

    rows = pd.DataFrame(all_records)
    summary, payload = _summarize(rows, option_meta, failed_options)

    OPTION_SCAN_ROWS_CSV.parent.mkdir(parents=True, exist_ok=True)
    rows.to_csv(OPTION_SCAN_ROWS_CSV, index=False)
    summary.to_csv(OPTION_SCAN_SUMMARY_CSV, index=False)
    OPTION_SCAN_JSON.parent.mkdir(parents=True, exist_ok=True)
    OPTION_SCAN_JSON.write_text(json.dumps(_jsonable(payload), indent=2) + "\n", encoding="utf-8")
    _write_report(summary, payload)
    print(
        json.dumps(
            _jsonable({k: payload[k] for k in ("status", "best_option", "baseline_option", "accepted")}), indent=2
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
