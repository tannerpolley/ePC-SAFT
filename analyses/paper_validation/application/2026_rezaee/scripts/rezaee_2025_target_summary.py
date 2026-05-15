from __future__ import annotations

import csv
from pathlib import Path
from typing import Any

ANALYSIS_DIR = Path(__file__).resolve().parents[1]
INPUT_DIR = ANALYSIS_DIR / "data" / "input"
PROCESSED_DIR = ANALYSIS_DIR / "data" / "processed"
RESULTS_DIR = ANALYSIS_DIR / "results" / "targets"

DOE_CSV = INPUT_DIR / "rezaee_2025_doe_extraction_responses.csv"
OPTIMUM_CSV = INPUT_DIR / "rezaee_2025_optimum_neighborhood.csv"
REAL_BRINE_CSV = INPUT_DIR / "rezaee_2025_real_brine_extraction.csv"
EQUILIBRIUM_CSV = INPUT_DIR / "rezaee_2025_extraction_equilibrium_mole_fractions.csv"

SUMMARY_CSV = PROCESSED_DIR / "rezaee_2025_extraction_target_summary.csv"
EQUILIBRIUM_SUMMARY_CSV = (
    PROCESSED_DIR / "rezaee_2025_extraction_equilibrium_summary.csv"
)
REPORT_MD = RESULTS_DIR / "rezaee_2025_extraction_target_summary.md"


def _rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _float(row: dict[str, str], key: str) -> float | None:
    value = row.get(key, "").strip()
    if not value:
        return None
    return float(value.replace(",", ""))


def _write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        raise ValueError(f"No rows to write for {path}")
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def _summarize_targets() -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []

    for row in _rows(DOE_CSV):
        rows.append(
            {
                "target_group": "doe_table5",
                "target_id": f"table5_run_{int(float(row['experiment_no'])):02d}",
                "T_C": _float(row, "T_C"),
                "pH": _float(row, "pH"),
                "topo_wt_pct": _float(row, "topo_wt_pct"),
                "na_li_mass_ratio": _float(row, "na_li_mass_ratio"),
                "li_extraction_pct_exp": _float(row, "li_extraction_pct_exp"),
                "na_extraction_pct_exp": _float(row, "na_extraction_pct_exp"),
                "li_na_selectivity_exp": _float(row, "li_na_selectivity_exp"),
                "source": row["source"],
                "notes": "Direct Li/Na extraction response target from Table 5.",
            }
        )

    for row in _rows(OPTIMUM_CSV):
        li = _float(row, "li_extraction_pct_exp")
        selectivity = _float(row, "li_na_selectivity_exp")
        inferred_na = li / selectivity if li is not None and selectivity else None
        rows.append(
            {
                "target_group": "optimum_table7",
                "target_id": f"table7_run_{int(float(row['experiment_no'])):02d}",
                "T_C": _float(row, "T_C"),
                "pH": _float(row, "pH"),
                "topo_wt_pct": _float(row, "topo_wt_pct"),
                "na_li_mass_ratio": _float(row, "na_li_mass_ratio"),
                "li_extraction_pct_exp": li,
                "na_extraction_pct_exp": inferred_na,
                "li_na_selectivity_exp": selectivity,
                "source": row["source"],
                "notes": "Table 7 reports Li extraction and selectivity; Na extraction is inferred as Li extraction divided by selectivity.",
            }
        )

    real_rows = {row["ion"]: row for row in _rows(REAL_BRINE_CSV)}
    li_row = real_rows["Li"]
    na_row = real_rows["Na"]
    li_eff = _float(li_row, "extraction_efficiency_pct")
    na_eff = _float(na_row, "extraction_efficiency_pct")
    rows.append(
        {
            "target_group": "real_brine_table8",
            "target_id": "table8_li_na",
            "T_C": "",
            "pH": "",
            "topo_wt_pct": 10.0,
            "na_li_mass_ratio": "",
            "li_extraction_pct_exp": li_eff,
            "na_extraction_pct_exp": na_eff,
            "li_na_selectivity_exp": (
                li_eff / na_eff if li_eff is not None and na_eff else ""
            ),
            "source": "Rezaee2025_Table8",
            "notes": "Li/Na extraction target from synthetic brine modeled on an Iran source brine under the selected operating conditions.",
        }
    )

    return rows


def _summarize_equilibrium_rows() -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for row in _rows(EQUILIBRIUM_CSV):
        aq_li = _float(row, "aqueous_x_Li_plus")
        aq_na = _float(row, "aqueous_x_Na_plus")
        org_li = _float(row, "organic_x_RLi")
        org_na = _float(row, "organic_x_RNa")
        aq_sum = sum(
            _float(row, key) or 0.0 for key in row if key.startswith("aqueous_x_")
        )
        org_sum = sum(
            _float(row, key) or 0.0 for key in row if key.startswith("organic_x_")
        )
        li_distribution = org_li / aq_li if org_li is not None and aq_li else None
        na_distribution = org_na / aq_na if org_na is not None and aq_na else None
        rows.append(
            {
                "experiment_no": int(float(row["experiment_no"])),
                "aqueous_mole_fraction_sum": aq_sum,
                "organic_mole_fraction_sum": org_sum,
                "li_distribution_molfrac_ratio": li_distribution,
                "na_distribution_molfrac_ratio": na_distribution,
                "li_na_selectivity_molfrac_proxy": (
                    li_distribution / na_distribution
                    if li_distribution is not None and na_distribution
                    else None
                ),
                "source": row["source"],
                "notes": (
                    "Diagnostic mole-fraction ratio from SI Tables S1/S2. Use the raw phase mole fractions for fitting; "
                    "this proxy is not the paper extraction percentage because phase amounts are not included here."
                ),
            }
        )
    return rows


def _write_report(
    summary_rows: list[dict[str, Any]], equilibrium_rows: list[dict[str, Any]]
) -> None:
    doe_rows = [row for row in summary_rows if row["target_group"] == "doe_table5"]
    optimum_rows = [
        row for row in summary_rows if row["target_group"] == "optimum_table7"
    ]
    real_row = next(
        row for row in summary_rows if row["target_group"] == "real_brine_table8"
    )
    selected = next(row for row in optimum_rows if row["target_id"] == "table7_run_03")

    lines = [
        "# Rezaee 2025 Extraction Target Summary",
        "",
        "## Scope",
        "",
        "This report summarizes source-backed extraction-percentage targets that can be used to compare a Rezaee-style reactive extraction calculation against the paper.",
        "",
        "The SI equilibrium-composition data are now tracked separately in `data/input/rezaee_2025_extraction_equilibrium_mole_fractions.csv` and summarized in `data/processed/rezaee_2025_extraction_equilibrium_summary.csv`.",
        "",
        "## Target Counts",
        "",
        f"- Table 5 DOE Li/Na response rows: `{len(doe_rows)}`.",
        f"- Table 7 optimum-neighborhood rows: `{len(optimum_rows)}`.",
        "- Table 8 real-brine Li/Na comparison rows: `1` combined Li/Na target.",
        f"- SI Tables S1/S2 equilibrium mole-fraction rows: `{len(equilibrium_rows)}`.",
        "",
        "## Selected Operating Point",
        "",
        f"- `T_C`: `{selected['T_C']}`.",
        f"- `pH`: `{selected['pH']}`.",
        f"- `topo_wt_pct`: `{selected['topo_wt_pct']}`.",
        f"- `na_li_mass_ratio`: `{selected['na_li_mass_ratio']}`.",
        f"- `li_extraction_pct_exp`: `{selected['li_extraction_pct_exp']}`.",
        f"- `na_extraction_pct_exp_inferred`: `{selected['na_extraction_pct_exp']}`.",
        f"- `li_na_selectivity_exp`: `{selected['li_na_selectivity_exp']}`.",
        "",
        "## Iran-Source Model Brine Check",
        "",
        f"- `li_extraction_pct_exp`: `{real_row['li_extraction_pct_exp']}`.",
        f"- `na_extraction_pct_exp`: `{real_row['na_extraction_pct_exp']}`.",
        f"- `li_na_selectivity_exp`: `{real_row['li_na_selectivity_exp']}`.",
        "",
        "## Remaining Model Gate",
        "",
        "The response and equilibrium-composition rows are enough to start a source-backed Li/Na fitting target. The remaining gap is model closure: connect the aqueous ions, organic DES/TOPO pseudo-components, and RLi/RNa complex species through a calibrated reactive/electrolyte LLE objective before claiming extraction-percent reproduction.",
    ]
    REPORT_MD.parent.mkdir(parents=True, exist_ok=True)
    REPORT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    summary_rows = _summarize_targets()
    equilibrium_rows = _summarize_equilibrium_rows()
    _write_csv(SUMMARY_CSV, summary_rows)
    _write_csv(EQUILIBRIUM_SUMMARY_CSV, equilibrium_rows)
    _write_report(summary_rows, equilibrium_rows)
    print(f"Wrote {SUMMARY_CSV}")
    print(f"Wrote {EQUILIBRIUM_SUMMARY_CSV}")
    print(f"Wrote {REPORT_MD}")
    print({"target_rows": len(summary_rows), "equilibrium_rows": len(equilibrium_rows)})
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
