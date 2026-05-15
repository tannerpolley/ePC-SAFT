from __future__ import annotations

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

import rezaee_reactive_epcsaft_option_scan as option_scan  # noqa: E402
import rezaee_reactive_equilibrium_replay as replay  # noqa: E402

ANALYSIS_DIR = Path(__file__).resolve().parents[1]
INPUT_DIR = ANALYSIS_DIR / "data" / "input"
PROCESSED_DIR = ANALYSIS_DIR / "data" / "processed"
RESULTS_DIR = ANALYSIS_DIR / "results" / "reaction_equilibrium"

DOE_CSV = INPUT_DIR / "rezaee_2025_doe_extraction_responses.csv"
BENCHMARK_ROWS_CSV = PROCESSED_DIR / "rezaee_2026_paper_basis_reaction_coordinate_rows.csv"
BENCHMARK_SUMMARY_JSON = RESULTS_DIR / "rezaee_2026_paper_basis_reaction_coordinate_summary.json"
BENCHMARK_REPORT_MD = RESULTS_DIR / "rezaee_2026_paper_basis_reaction_coordinate.md"

LI_MW_KG_MOL = 6.94e-3
NA_MW_KG_MOL = 22.989769e-3
CL_MW_KG_MOL = 35.45e-3
H2O_MW_KG_MOL = 18.01528e-3
DES_MW_KG_MOL = 0.20748
TOPO_MW_KG_MOL = 0.38664
DOE_LI_CONCENTRATION_KG_PER_KG_AQUEOUS = 1000.0e-6
PHASE_MASS_KG = 1.0
PRESSURE_PA = 101325.0


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


def _initial_moles(row: Any) -> dict[str, float]:
    li_mass = DOE_LI_CONCENTRATION_KG_PER_KG_AQUEOUS * PHASE_MASS_KG
    na_mass = float(row.na_li_mass_ratio) * li_mass
    n_li = li_mass / LI_MW_KG_MOL
    n_na = na_mass / NA_MW_KG_MOL
    n_cl = n_li + n_na
    cl_mass = n_cl * CL_MW_KG_MOL
    water_mass = max(PHASE_MASS_KG - li_mass - na_mass - cl_mass, 0.5)
    n_h2o = water_mass / H2O_MW_KG_MOL
    n_h = 10.0 ** (-float(row.pH))
    n_oh = 10.0 ** (float(row.pH) - 14.0)
    n_nh4 = max(n_cl + n_oh - n_li - n_na - n_h, 0.0)

    topo_mass = PHASE_MASS_KG * float(row.topo_wt_pct) / 100.0
    des_mass = PHASE_MASS_KG - topo_mass
    return {
        "H2O": n_h2o,
        "Li+": n_li,
        "Na+": n_na,
        "Cl-": n_cl,
        "H+": n_h,
        "OH-": n_oh,
        "NH4+": n_nh4,
        "DES": des_mass / DES_MW_KG_MOL,
        "TOPO": topo_mass / TOPO_MW_KG_MOL,
        "RLi": 0.0,
        "RNa": 0.0,
    }


def _phase_x(moles: dict[str, float], labels: list[str]) -> np.ndarray:
    values = np.asarray([max(float(moles[label]), 1.0e-300) for label in labels], dtype=float)
    return values / float(np.sum(values))


def _extent_to_moles(initial: dict[str, float], xi_li: float, xi_na: float) -> dict[str, float]:
    out = dict(initial)
    out["Li+"] -= xi_li
    out["Na+"] -= xi_na
    out["OH-"] -= xi_li + xi_na
    out["DES"] -= xi_li + xi_na
    out["H2O"] += xi_li + xi_na
    out["RLi"] += xi_li
    out["RNa"] += xi_na
    return out


def _activity_coefficients(
    mode: str,
    aqueous_mix: Any,
    organic_mix: Any,
    pure_ln_phi: np.ndarray,
    aqueous_x: np.ndarray,
    organic_x: np.ndarray,
    temperature_k: float,
) -> tuple[dict[str, float], dict[str, float]]:
    if mode == "ideal":
        return ({label: 1.0 for label in replay.AQ_LABELS}, {label: 1.0 for label in replay.ORG_LABELS})
    aqueous_gamma = aqueous_mix.state(
        T=temperature_k,
        x=aqueous_x,
        P=PRESSURE_PA,
    ).activity_coefficient(species=replay.AQ_LABELS)
    organic_gamma = replay._organic_activity_coefficients(organic_mix, pure_ln_phi, organic_x)
    return aqueous_gamma, organic_gamma


def _solve_row(row: Any, mode: str, aqueous_mix: Any, organic_mix: Any, pure_ln_phi: np.ndarray) -> dict[str, Any]:
    constants = replay._reaction_constants()
    temperature_k = float(row.T_C) + 273.15
    initial = _initial_moles(row)
    organic_total = initial["DES"] + initial["TOPO"]
    xi_li = 1.0e-12
    xi_na = 1.0e-12
    max_extent = max(min(initial["OH-"], initial["DES"], initial["Li+"] + initial["Na+"]), 1.0e-300)
    converged = False
    delta = math.inf

    for iteration in range(1, 300):
        moles = _extent_to_moles(initial, xi_li, xi_na)
        aqueous_x = _phase_x(moles, replay.AQ_LABELS)
        organic_x = _phase_x(moles, replay.ORG_LABELS)
        aqueous_gamma, organic_gamma = _activity_coefficients(
            mode,
            aqueous_mix,
            organic_mix,
            pure_ln_phi,
            aqueous_x,
            organic_x,
            temperature_k,
        )
        x_rli = (
            constants["Li"]
            * (
                (aqueous_x[1] * aqueous_gamma["Li+"])
                * (aqueous_x[5] * aqueous_gamma["OH-"])
                * (organic_x[0] * organic_gamma["DES"])
            )
            / ((organic_gamma["RLi"]) * (aqueous_x[0] * aqueous_gamma["H2O"]))
        )
        x_rna = (
            constants["Na"]
            * (
                (aqueous_x[2] * aqueous_gamma["Na+"])
                * (aqueous_x[5] * aqueous_gamma["OH-"])
                * (organic_x[0] * organic_gamma["DES"])
            )
            / ((organic_gamma["RNa"]) * (aqueous_x[0] * aqueous_gamma["H2O"]))
        )
        next_li = min(max(float(x_rli) * organic_total, 0.0), initial["Li+"], initial["OH-"], initial["DES"])
        next_na = min(
            max(float(x_rna) * organic_total, 0.0),
            initial["Na+"],
            max(initial["OH-"] - next_li, 0.0),
            max(initial["DES"] - next_li, 0.0),
        )
        delta = max(abs(next_li - xi_li), abs(next_na - xi_na))
        xi_li = 0.5 * xi_li + 0.5 * next_li
        xi_na = 0.5 * xi_na + 0.5 * next_na
        if delta < 1.0e-12:
            converged = True
            break

    li_extraction = 100.0 * xi_li / initial["Li+"]
    na_extraction = 100.0 * xi_na / initial["Na+"]
    return {
        "experiment_no": int(row.experiment_no),
        "mode": mode,
        "T_C": float(row.T_C),
        "pH": float(row.pH),
        "topo_wt_pct": float(row.topo_wt_pct),
        "na_li_mass_ratio": float(row.na_li_mass_ratio),
        "initial_Li_mol": float(initial["Li+"]),
        "initial_Na_mol": float(initial["Na+"]),
        "initial_OH_mol_from_pH": float(initial["OH-"]),
        "initial_DES_mol": float(initial["DES"]),
        "oh_limited_max_li_plus_na_extraction_mol": float(max_extent),
        "li_extraction_pct_calc": float(li_extraction),
        "na_extraction_pct_calc": float(na_extraction),
        "li_extraction_pct_exp": float(row.li_extraction_pct_exp),
        "na_extraction_pct_exp": float(row.na_extraction_pct_exp),
        "li_abs_pct_error": float(abs(li_extraction - float(row.li_extraction_pct_exp))),
        "na_abs_pct_error": float(abs(na_extraction - float(row.na_extraction_pct_exp))),
        "selectivity_calc": float(li_extraction / max(na_extraction, 1.0e-300)),
        "selectivity_exp": float(row.li_na_selectivity_exp),
        "converged": bool(converged),
        "iterations": int(iteration),
        "final_extent_delta": float(delta),
        "source": row.source,
    }


def _run_mode(mode: str) -> list[dict[str, Any]]:
    organic_mix, _organic_params, pure_ln_phi = replay._organic_mixture()
    aqueous_mix = None
    if mode != "ideal":
        option = next(item for item in option_scan.OPTION_SETS if item["option_id"] == mode)
        aqueous_mix, _resolved, _snapshot = option_scan._build_aqueous_mixture(option)
    return [
        _solve_row(row, mode, aqueous_mix, organic_mix, pure_ln_phi)
        for row in pd.read_csv(DOE_CSV).itertuples(index=False)
    ]


def _summarize(rows: pd.DataFrame) -> dict[str, Any]:
    grouped = rows.groupby("mode", as_index=False).agg(
        row_count=("experiment_no", "count"),
        mean_li_abs_pct_error=("li_abs_pct_error", "mean"),
        mean_na_abs_pct_error=("na_abs_pct_error", "mean"),
        median_li_extraction_pct_calc=("li_extraction_pct_calc", "median"),
        median_na_extraction_pct_calc=("na_extraction_pct_calc", "median"),
        median_initial_oh_mol=("initial_OH_mol_from_pH", "median"),
        median_initial_li_mol=("initial_Li_mol", "median"),
        median_oh_to_li_mol_ratio=(
            "initial_OH_mol_from_pH",
            lambda values: float(np.median(np.asarray(values) / rows.loc[values.index, "initial_Li_mol"].to_numpy())),
        ),
    )
    best = grouped.sort_values("mean_li_abs_pct_error").iloc[0].to_dict()
    return {
        "status": "paper_basis_reaction_coordinate_benchmark_complete",
        "mode_count": int(grouped.shape[0]),
        "row_count": int(rows["experiment_no"].nunique()),
        "best_mode_by_li_abs_error": best,
        "summary_rows": grouped.to_dict(orient="records"),
        "interpretation": (
            "Literal pH-derived OH- stoichiometry makes Reaction A/B OH-limited by orders of magnitude. "
            "This benchmark therefore does not reproduce the reported extraction percentages with the published "
            "K values before any parameter refit."
        ),
    }


def _write_report(summary: dict[str, Any]) -> None:
    best = summary["best_mode_by_li_abs_error"]
    lines = [
        "# Rezaee 2026 Paper-Basis Reaction-Coordinate Benchmark",
        "",
        "## Purpose",
        "",
        "This benchmark follows the Section 3.2 reaction-coordinate idea directly from the 2025 DOE basis: equal aqueous/organic masses, 1000 ppm Li in the aqueous feed, Na/Li mass ratio from Table 5, pH-derived H+/OH-, TOPO wt% in the organic phase, and published Rezaee reaction constants and organic PC-SAFT parameters. It is deliberately separate from the package phase-equilibrium smoke.",
        "",
        "## Result",
        "",
        f"- Modes tested: `{summary['mode_count']}`.",
        f"- Rows tested: `{summary['row_count']}`.",
        f"- Best mode by Li extraction absolute error: `{best['mode']}`.",
        f"- Mean Li extraction absolute error: `{float(best['mean_li_abs_pct_error']):.6g}` percentage points.",
        f"- Mean Na extraction absolute error: `{float(best['mean_na_abs_pct_error']):.6g}` percentage points.",
        f"- Median calculated Li extraction: `{float(best['median_li_extraction_pct_calc']):.6g}%`.",
        f"- Median calculated Na extraction: `{float(best['median_na_extraction_pct_calc']):.6g}%`.",
        f"- Median initial OH/Li mole ratio from pH basis: `{float(best['median_oh_to_li_mol_ratio']):.6g}`.",
        "",
        "## Interpretation",
        "",
        summary["interpretation"],
        "",
        "The pH-derived hydroxide amount is far too small to support 30-50% lithium extraction if Eq. 14/15 is interpreted as a literal one-OH-minus stoichiometric reaction extent. This points to an apparent-activity/reference-state interpretation in the paper workflow, not a simple ePC-SAFT Born/dielectric option issue.",
        "",
        "## Generated Files",
        "",
        f"- `{BENCHMARK_ROWS_CSV.relative_to(ANALYSIS_DIR)}`",
        f"- `{BENCHMARK_SUMMARY_JSON.relative_to(ANALYSIS_DIR)}`",
    ]
    BENCHMARK_REPORT_MD.parent.mkdir(parents=True, exist_ok=True)
    BENCHMARK_REPORT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    modes = ["ideal", "package_struct_default_no_explicit_elec_model", "2025_born_no_ssm_empirical"]
    records = []
    for mode in modes:
        records.extend(_run_mode(mode))
    rows = pd.DataFrame(records)
    summary = _summarize(rows)
    BENCHMARK_ROWS_CSV.parent.mkdir(parents=True, exist_ok=True)
    rows.to_csv(BENCHMARK_ROWS_CSV, index=False)
    BENCHMARK_SUMMARY_JSON.parent.mkdir(parents=True, exist_ok=True)
    BENCHMARK_SUMMARY_JSON.write_text(json.dumps(_jsonable(summary), indent=2) + "\n", encoding="utf-8")
    _write_report(summary)
    print(json.dumps(_jsonable(summary), indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
