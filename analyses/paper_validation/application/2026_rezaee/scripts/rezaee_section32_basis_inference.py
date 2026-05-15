from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

ANALYSIS_DIR = Path(__file__).resolve().parents[1]
INPUT_DIR = ANALYSIS_DIR / "data" / "input"
PROCESSED_DIR = ANALYSIS_DIR / "data" / "processed"
RESULTS_DIR = ANALYSIS_DIR / "results" / "reaction_equilibrium"

DOE_CSV = INPUT_DIR / "rezaee_2025_doe_extraction_responses.csv"
EQUILIBRIUM_CSV = INPUT_DIR / "rezaee_2025_extraction_equilibrium_mole_fractions.csv"
ROWS_CSV = PROCESSED_DIR / "rezaee_2026_section32_basis_inference_rows.csv"
SUMMARY_JSON = RESULTS_DIR / "rezaee_2026_section32_basis_inference_summary.json"
REPORT_MD = RESULTS_DIR / "rezaee_2026_section32_basis_inference.md"

LI_MW_KG_MOL = 6.94e-3
NA_MW_KG_MOL = 22.989769e-3
DOE_LI_CONCENTRATION_KG_PER_KG_AQUEOUS = 1000.0e-6
WATER_MOLARITY_APPROX = 55.508
BENCHMARK_ROW_COUNT = 26


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


def _safe_div(numerator: float, denominator: float) -> float:
    return float(numerator) / max(float(denominator), 1.0e-300)


def _rows() -> pd.DataFrame:
    doe = pd.read_csv(DOE_CSV)
    equilibrium = pd.read_csv(EQUILIBRIUM_CSV)
    if doe["experiment_no"].nunique() != BENCHMARK_ROW_COUNT:
        raise ValueError(f"Expected {BENCHMARK_ROW_COUNT} DOE rows in {DOE_CSV}")
    if equilibrium["experiment_no"].nunique() != BENCHMARK_ROW_COUNT:
        raise ValueError(f"Expected {BENCHMARK_ROW_COUNT} equilibrium rows in {EQUILIBRIUM_CSV}")
    merged = doe.merge(equilibrium, on="experiment_no", suffixes=("_doe", "_eq"))

    records: list[dict[str, Any]] = []
    for row in merged.itertuples(index=False):
        initial_li_mol = DOE_LI_CONCENTRATION_KG_PER_KG_AQUEOUS / LI_MW_KG_MOL
        initial_na_mol = float(row.na_li_mass_ratio) * DOE_LI_CONCENTRATION_KG_PER_KG_AQUEOUS / NA_MW_KG_MOL
        li_aq_mol = initial_li_mol * (1.0 - float(row.li_extraction_pct_exp) / 100.0)
        na_aq_mol = initial_na_mol * (1.0 - float(row.na_extraction_pct_exp) / 100.0)
        li_org_mol = initial_li_mol - li_aq_mol
        na_org_mol = initial_na_mol - na_aq_mol

        n_aq_from_li = _safe_div(li_aq_mol, row.aqueous_x_Li_plus)
        n_aq_from_na = _safe_div(na_aq_mol, row.aqueous_x_Na_plus)
        n_aq_mean = 0.5 * (n_aq_from_li + n_aq_from_na)
        n_org_from_rli = _safe_div(li_org_mol, row.organic_x_RLi)
        n_org_from_rna = _safe_div(na_org_mol, row.organic_x_RNa)

        aqueous_charge_residual = (
            row.aqueous_x_Li_plus
            + row.aqueous_x_Na_plus
            - row.aqueous_x_Cl_minus
            + row.aqueous_x_H_plus
            - row.aqueous_x_OH_minus
            + row.aqueous_x_NH4_plus
        )
        x_oh_from_reported_ph = 10.0 ** (float(row.pH) - 14.0) / WATER_MOLARITY_APPROX
        h_mol_from_si = float(row.aqueous_x_H_plus) * n_aq_mean
        oh_mol_from_si = float(row.aqueous_x_OH_minus) * n_aq_mean
        nh4_mol_from_si = float(row.aqueous_x_NH4_plus) * n_aq_mean

        records.append(
            {
                "experiment_no": int(row.experiment_no),
                "T_C": float(row.T_C),
                "pH": float(row.pH),
                "topo_wt_pct": float(row.topo_wt_pct),
                "na_li_mass_ratio": float(row.na_li_mass_ratio),
                "Li_extraction_pct_exp": float(row.li_extraction_pct_exp),
                "Na_extraction_pct_exp": float(row.na_extraction_pct_exp),
                "selectivity_exp": float(row.li_na_selectivity_exp),
                "initial_Li_mol_from_1000ppm_per_kg_aqueous": float(initial_li_mol),
                "initial_Na_mol_from_mass_ratio": float(initial_na_mol),
                "Li_aqueous_mol_from_extraction": float(li_aq_mol),
                "Na_aqueous_mol_from_extraction": float(na_aq_mol),
                "Li_organic_mol_from_extraction": float(li_org_mol),
                "Na_organic_mol_from_extraction": float(na_org_mol),
                "N_aq_total_from_Li": float(n_aq_from_li),
                "N_aq_total_from_Na": float(n_aq_from_na),
                "N_aq_total_Li_over_Na": _safe_div(n_aq_from_li, n_aq_from_na),
                "N_org_total_from_RLi": float(n_org_from_rli),
                "N_org_total_from_RNa": float(n_org_from_rna),
                "N_org_total_RLi_over_RNa": _safe_div(n_org_from_rli, n_org_from_rna),
                "aqueous_charge_residual": float(aqueous_charge_residual),
                "xOH_SI": float(row.aqueous_x_OH_minus),
                "xOH_from_reported_pH_water_molarity": float(x_oh_from_reported_ph),
                "xOH_SI_over_pH_estimate": _safe_div(row.aqueous_x_OH_minus, x_oh_from_reported_ph),
                "H_mol_from_SI_mean_aq_total": float(h_mol_from_si),
                "OH_mol_from_SI_mean_aq_total": float(oh_mol_from_si),
                "NH4_mol_from_SI_mean_aq_total": float(nh4_mol_from_si),
                "H_mol_over_extracted_Li_plus_Na": _safe_div(h_mol_from_si, li_org_mol + na_org_mol),
                "OH_mol_over_extracted_Li_plus_Na": _safe_div(oh_mol_from_si, li_org_mol + na_org_mol),
                "source": row.source_eq,
            }
        )
    return pd.DataFrame(records)


def _summarize(rows: pd.DataFrame) -> dict[str, Any]:
    return {
        "status": "section32_basis_inference_complete",
        "row_count": int(rows["experiment_no"].nunique()),
        "conservation_diagnostics": {
            "median_N_aq_total_Li_over_Na": float(rows["N_aq_total_Li_over_Na"].median()),
            "min_N_aq_total_Li_over_Na": float(rows["N_aq_total_Li_over_Na"].min()),
            "max_N_aq_total_Li_over_Na": float(rows["N_aq_total_Li_over_Na"].max()),
            "median_N_org_total_RLi_over_RNa": float(rows["N_org_total_RLi_over_RNa"].median()),
            "min_N_org_total_RLi_over_RNa": float(rows["N_org_total_RLi_over_RNa"].min()),
            "max_N_org_total_RLi_over_RNa": float(rows["N_org_total_RLi_over_RNa"].max()),
            "max_abs_aqueous_charge_residual": float(np.max(np.abs(rows["aqueous_charge_residual"]))),
        },
        "basis_diagnostics": {
            "median_xOH_SI_over_pH_estimate": float(rows["xOH_SI_over_pH_estimate"].median()),
            "median_H_mol_over_extracted_Li_plus_Na": float(rows["H_mol_over_extracted_Li_plus_Na"].median()),
            "median_OH_mol_over_extracted_Li_plus_Na": float(rows["OH_mol_over_extracted_Li_plus_Na"].median()),
            "median_N_aq_total_from_Li": float(rows["N_aq_total_from_Li"].median()),
            "median_N_aq_total_from_Na": float(rows["N_aq_total_from_Na"].median()),
            "median_N_org_total_from_RLi": float(rows["N_org_total_from_RLi"].median()),
            "median_N_org_total_from_RNa": float(rows["N_org_total_from_RNa"].median()),
        },
        "interpretation": (
            "The SI aqueous rows are nearly charge balanced and the OH- mole fractions track the reported pH. "
            "However, combining Table 5 extraction percentages with SI organic RLi/RNa mole fractions does not "
            "produce one consistent organic phase total for one-metal RLi/RNa complexes. This is a source-basis "
            "blocker for exact Section 3.2 reproduction until the original initial-mole/reaction-coordinate "
            "worksheet or a documented phase-amount convention is available."
        ),
    }


def _write_report(summary: dict[str, Any]) -> None:
    conservation = summary["conservation_diagnostics"]
    basis = summary["basis_diagnostics"]
    lines = [
        "# Rezaee Section 3.2 Initial-Basis Inference",
        "",
        "## Purpose",
        "",
        "This diagnostic tests whether the 2025 Table 5 extraction percentages and 2025 SI aqueous/organic equilibrium mole fractions define a self-consistent mole basis for the Rezaee 2026 Eq. 17 reaction-coordinate calculation.",
        "",
        "The check does not fit parameters and does not use package equilibrium solvers. It only applies metal conservation to infer phase totals from Li, Na, RLi, and RNa rows.",
        "",
        "## Result",
        "",
        f"- Rows evaluated: `{summary['row_count']}`.",
        f"- Median aqueous total consistency, Li-derived / Na-derived: `{conservation['median_N_aq_total_Li_over_Na']}`.",
        f"- Median organic total consistency, RLi-derived / RNa-derived: `{conservation['median_N_org_total_RLi_over_RNa']}`.",
        f"- Max absolute aqueous charge residual: `{conservation['max_abs_aqueous_charge_residual']}`.",
        f"- Median SI xOH / pH-derived xOH estimate: `{basis['median_xOH_SI_over_pH_estimate']}`.",
        f"- Median SI H moles / extracted Li+Na moles: `{basis['median_H_mol_over_extracted_Li_plus_Na']}`.",
        f"- Median SI OH moles / extracted Li+Na moles: `{basis['median_OH_mol_over_extracted_Li_plus_Na']}`.",
        "",
        "## Interpretation",
        "",
        summary["interpretation"],
        "",
        "The direct replication should therefore keep using the 26 rows as the benchmark set, but it should not claim the paper's post-Table-9 AARD values until the hidden initial-mole/phase-amount convention is supplied or reconstructed from a source-backed worksheet.",
        "",
        "## Generated Files",
        "",
        f"- `{ROWS_CSV.relative_to(ANALYSIS_DIR)}`",
        f"- `{SUMMARY_JSON.relative_to(ANALYSIS_DIR)}`",
    ]
    REPORT_MD.parent.mkdir(parents=True, exist_ok=True)
    REPORT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    rows = _rows()
    summary = _summarize(rows)
    ROWS_CSV.parent.mkdir(parents=True, exist_ok=True)
    rows.to_csv(ROWS_CSV, index=False)
    SUMMARY_JSON.parent.mkdir(parents=True, exist_ok=True)
    SUMMARY_JSON.write_text(json.dumps(_jsonable(summary), indent=2) + "\n", encoding="utf-8")
    _write_report(summary)
    print(json.dumps(_jsonable(summary), indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
