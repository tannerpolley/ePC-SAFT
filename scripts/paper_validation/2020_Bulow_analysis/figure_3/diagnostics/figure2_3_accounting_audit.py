from __future__ import annotations

import csv
import math
import sys
from pathlib import Path

import numpy as np


SCRIPT_DIR = Path(__file__).resolve().parent
FIGURE_DIR = SCRIPT_DIR.parent
ANALYSIS_ROOT = FIGURE_DIR.parent
REPO_ROOT = ANALYSIS_ROOT.parents[2]
OUTPUT_DIR = SCRIPT_DIR

if str(ANALYSIS_ROOT) not in sys.path:
    sys.path.insert(0, str(ANALYSIS_ROOT))
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

import _model_overlay as overlay
import _plot_common as common


R_GAS = 8.31446261815324
T_REF = 298.15
FIG2_PATH = ANALYSIS_ROOT / "figure_2" / "data" / "water_comparisons.csv"
FIG3_PATH = ANALYSIS_ROOT / "figure_3" / "data" / "water_contributions.csv"
TERMS = ("hc", "disp", "assoc", "born")
TERM_KEY_MAP = {
    "hc": "hc",
    "disp": "disp",
    "assoc": "assoc",
    "born": "born",
}
PAPER_ROW_KEYS = {
    "hc": ("hc avg", "hc"),
    "disp": ("disp avg", "disp"),
    "assoc": ("assoc avg", "assoc"),
    "born": ("born avg", "born"),
}


def _safe_corr(a: np.ndarray, b: np.ndarray) -> float:
    if a.size < 2 or b.size < 2:
        return math.nan
    if np.allclose(a, a[0]) or np.allclose(b, b[0]):
        return math.nan
    return float(np.corrcoef(a, b)[0, 1])


def _fit_line(x: np.ndarray, y: np.ndarray) -> tuple[float, float, float]:
    if (not np.all(np.isfinite(x))) or (not np.all(np.isfinite(y))):
        return math.nan, math.nan, math.nan
    if x.size < 2 or np.allclose(x, x[0]):
        return math.nan, math.nan, math.nan
    coeffs = np.polyfit(x, y, 1)
    slope = float(coeffs[0])
    intercept = float(coeffs[1])
    fit = slope * x + intercept
    rmse = float(np.sqrt(np.mean((fit - y) ** 2)))
    return slope, intercept, rmse


def _stable_logz_over_zminus1(z_total: float) -> float:
    if abs(z_total - 1.0) < 1.0e-10:
        return 1.0
    return math.log(z_total) / (z_total - 1.0)


def _paper_values(frame: common.IndexedCsv, term_key: str, ions: list[str]) -> np.ndarray:
    for row_key in PAPER_ROW_KEYS[term_key]:
        if row_key in frame.index:
            return frame.values(row_key, ions)
    raise KeyError(f"Missing paper row for {term_key}")


def _model_state(ion: str) -> tuple[dict[str, object], int]:
    return overlay._infinite_dilution_terms("advanced", ion, "water")


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    fig2 = common.load_indexed_csv(FIG2_PATH)
    fig3 = common.load_indexed_csv(FIG3_PATH)
    ions = list(fig3.columns)

    detail_rows: list[dict[str, object]] = []
    term_summary_rows: list[dict[str, object]] = []
    ion_summary_rows: list[dict[str, object]] = []

    paper_total = fig2.values("advanced", ions)
    paper_born = _paper_values(fig3, "born", ions)
    paper_short = (
        _paper_values(fig3, "hc", ions)
        + _paper_values(fig3, "disp", ions)
        + _paper_values(fig3, "assoc", ions)
    )
    paper_sum4 = paper_short + paper_born
    paper_missing_to_total = paper_total - paper_sum4
    paper_total_minus_born = paper_total - paper_born

    model_mu_by_term: dict[str, np.ndarray] = {term: np.empty(len(ions), dtype=float) for term in TERMS}
    model_lnf_by_term: dict[str, np.ndarray] = {term: np.empty(len(ions), dtype=float) for term in TERMS}
    model_corr_by_term: dict[str, np.ndarray] = {term: np.empty(len(ions), dtype=float) for term in TERMS}
    model_total = np.empty(len(ions), dtype=float)
    model_zcorr_total = np.empty(len(ions), dtype=float)

    for idx, ion in enumerate(ions):
        terms, ion_idx = _model_state(ion)
        z_total = float(terms["z_total"])
        z_factor = _stable_logz_over_zminus1(z_total)
        for term in TERMS:
            key = TERM_KEY_MAP[term]
            mu = float(R_GAS * T_REF * np.asarray(terms[f"mu_{key}"], dtype=float)[ion_idx] / 1000.0)
            lnf = float(R_GAS * T_REF * np.asarray(terms[f"lnfugcoef_{key}"], dtype=float)[ion_idx] / 1000.0)
            z_alpha = float(terms[f"z_{key}"])
            zcorr = float(R_GAS * T_REF * z_alpha * z_factor / 1000.0)
            model_mu_by_term[term][idx] = mu
            model_lnf_by_term[term][idx] = lnf
            model_corr_by_term[term][idx] = zcorr
            paper_val = float(_paper_values(fig3, term, [ion])[0])
            detail_rows.append(
                {
                    "ion": ion,
                    "term": term,
                    "paper": paper_val,
                    "model_mu": mu,
                    "model_lnfug": lnf,
                    "model_zcorr_share": zcorr,
                    "delta_paper_minus_mu": float(paper_val - mu),
                    "delta_paper_minus_lnfug": float(paper_val - lnf),
                }
            )

        model_total[idx] = float(overlay.gsolv_ion("advanced", ion, "water"))
        model_zcorr_total[idx] = float(np.sum([model_corr_by_term[term][idx] for term in TERMS]))

    model_short_mu = model_mu_by_term["hc"] + model_mu_by_term["disp"] + model_mu_by_term["assoc"]
    model_short_lnf = model_lnf_by_term["hc"] + model_lnf_by_term["disp"] + model_lnf_by_term["assoc"]
    model_short_corr = model_corr_by_term["hc"] + model_corr_by_term["disp"] + model_corr_by_term["assoc"]
    model_sum4_mu = model_short_mu + model_mu_by_term["born"]
    model_sum4_lnf = model_short_lnf + model_lnf_by_term["born"]

    for term in TERMS:
        paper = np.asarray([row["paper"] for row in detail_rows if row["term"] == term], dtype=float)
        candidate_map = {
            "mu": model_mu_by_term[term],
            "lnfug": model_lnf_by_term[term],
            "zcorr_share": model_corr_by_term[term],
            "neg_mu": -model_mu_by_term[term],
            "neg_lnfug": -model_lnf_by_term[term],
        }
        for candidate_name, candidate in candidate_map.items():
            slope, intercept, fit_rmse = _fit_line(candidate, paper)
            term_summary_rows.append(
                {
                    "term": term,
                    "candidate": candidate_name,
                    "rmse": float(np.sqrt(np.mean((candidate - paper) ** 2))),
                    "mae": float(np.mean(np.abs(candidate - paper))),
                    "corr": _safe_corr(candidate, paper),
                    "fit_slope": slope,
                    "fit_intercept": intercept,
                    "fit_rmse": fit_rmse,
                }
            )

    for idx, ion in enumerate(ions):
        ion_summary_rows.append(
            {
                "ion": ion,
                "paper_total": float(paper_total[idx]),
                "paper_sum_hc_disp_assoc_born": float(paper_sum4[idx]),
                "paper_total_minus_born": float(paper_total_minus_born[idx]),
                "paper_short_sum": float(paper_short[idx]),
                "paper_missing_to_total": float(paper_missing_to_total[idx]),
                "model_total_lnfug": float(model_total[idx]),
                "model_sum_hc_disp_assoc_born_lnfug": float(model_sum4_lnf[idx]),
                "model_sum_hc_disp_assoc_born_mu": float(model_sum4_mu[idx]),
                "model_short_sum_lnfug": float(model_short_lnf[idx]),
                "model_short_sum_mu": float(model_short_mu[idx]),
                "model_total_zcorr_share": float(model_zcorr_total[idx]),
            }
        )

    detail_path = OUTPUT_DIR / "figure2_3_accounting_detail.csv"
    with detail_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "ion",
                "term",
                "paper",
                "model_mu",
                "model_lnfug",
                "model_zcorr_share",
                "delta_paper_minus_mu",
                "delta_paper_minus_lnfug",
            ],
        )
        writer.writeheader()
        writer.writerows(detail_rows)

    term_summary_path = OUTPUT_DIR / "figure2_3_accounting_term_summary.csv"
    with term_summary_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["term", "candidate", "rmse", "mae", "corr", "fit_slope", "fit_intercept", "fit_rmse"],
        )
        writer.writeheader()
        writer.writerows(sorted(term_summary_rows, key=lambda row: (row["term"], row["rmse"])))

    ion_summary_path = OUTPUT_DIR / "figure2_3_accounting_ion_summary.csv"
    with ion_summary_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "ion",
                "paper_total",
                "paper_sum_hc_disp_assoc_born",
                "paper_total_minus_born",
                "paper_short_sum",
                "paper_missing_to_total",
                "model_total_lnfug",
                "model_sum_hc_disp_assoc_born_lnfug",
                "model_sum_hc_disp_assoc_born_mu",
                "model_short_sum_lnfug",
                "model_short_sum_mu",
                "model_total_zcorr_share",
            ],
        )
        writer.writeheader()
        writer.writerows(ion_summary_rows)

    summary_path = OUTPUT_DIR / "figure2_3_accounting_audit.md"
    best_by_term: dict[str, dict[str, object]] = {}
    for row in sorted(term_summary_rows, key=lambda entry: (entry["term"], entry["rmse"])):
        best_by_term.setdefault(str(row["term"]), row)

    with summary_path.open("w", encoding="utf-8") as handle:
        handle.write("# Figure 2-3 Accounting Audit\n\n")
        handle.write("This audit compares the paper Figure 3 hydration contributions against several model-side intermediate quantities at the same infinite-dilution state used for Figure 2.\n\n")

        handle.write("## Main Findings\n\n")
        handle.write(
            f"- The paper Figure 3 contributions do **not** sum to the paper Figure 2 advanced totals. The missing amount `paper total - (hc + disp + assoc + born)` ranges from `{float(np.min(paper_missing_to_total)):.3f}` to `{float(np.max(paper_missing_to_total)):.3f}` kJ/mol across the ions.\n"
        )
        handle.write(
            f"- The paper quantity `Figure 2 total - Figure 3 Born` matches the model short-range sum `hc + disp + assoc` on the correct $\\ln\\varphi$ basis very closely: RMSE = `{float(np.sqrt(np.mean((paper_total_minus_born - model_short_lnf) ** 2))):.3f}` kJ/mol.\n"
        )
        handle.write(
            f"- Individually, the paper `hc`, `disp`, and `assoc` bars match the model **raw** $\\tilde{{\\mu}}^\\alpha$ terms far better than the correct $\\ln\\varphi^\\alpha$ contributions. This is strongest evidence that the paper Figure 3 short-range bars were not plotted on the same bookkeeping basis as the Figure 2 totals.\n"
        )
        handle.write(
            f"- The Born bars are insensitive to that distinction here: Born RMSE is `{best_by_term['born']['rmse']:.3f}` kJ/mol for both `mu` and `lnfug`, which explains why Born appears consistent while the other contributions do not.\n\n"
        )

        handle.write("## Candidate Ranking By Term\n\n")
        handle.write("| term | best candidate | RMSE | MAE | corr |\n")
        handle.write("| --- | --- | --- | --- | --- |\n")
        for term in TERMS:
            row = best_by_term[term]
            handle.write(
                f"| {term} | {row['candidate']} | {row['rmse']:.3f} | {row['mae']:.3f} | {row['corr']:.3f} |\n"
            )
        handle.write("\n")

        handle.write("## Mu vs Lnfug Comparison\n\n")
        handle.write("| term | RMSE(mu) | RMSE(lnfug) | fit slope mu->paper | fit slope lnfug->paper |\n")
        handle.write("| --- | --- | --- | --- | --- |\n")
        for term in TERMS:
            mu_row = next(row for row in term_summary_rows if row["term"] == term and row["candidate"] == "mu")
            lnf_row = next(row for row in term_summary_rows if row["term"] == term and row["candidate"] == "lnfug")
            handle.write(
                f"| {term} | {mu_row['rmse']:.3f} | {lnf_row['rmse']:.3f} | {mu_row['fit_slope']:.3f} | {lnf_row['fit_slope']:.3f} |\n"
            )
        handle.write("\n")

        handle.write("## Ion-Level Total Check\n\n")
        handle.write("| ion | paper Fig2 total | paper Fig3 sum | missing to total | paper total - born | model short sum lnfug |\n")
        handle.write("| --- | --- | --- | --- | --- | --- |\n")
        for row in ion_summary_rows:
            handle.write(
                f"| {row['ion']} | {row['paper_total']:.3f} | {row['paper_sum_hc_disp_assoc_born']:.3f} | {row['paper_missing_to_total']:.3f} | {row['paper_total_minus_born']:.3f} | {row['model_short_sum_lnfug']:.3f} |\n"
            )
        handle.write("\n")

        handle.write("## Interpretation\n\n")
        handle.write(
            "- The most plausible explanation is that Figure 2 uses the correct Gibbs-energy quantity $RT\\ln\\varphi_i^\\infty$, while the non-Born Figure 3 bars were generated from a pre-compressibility intermediate closer to the raw $RT\\tilde{\\mu}_i^\\alpha$ terms.\n"
        )
        handle.write(
            "- Converting from $\\tilde{\\mu}_i^\\alpha$ to $\\ln\\varphi_i^\\alpha$ adds the explicit $Z^\\alpha$ share correction from Eq. `lnphi_alpha`. For hydration in water, those short-range $Z$-share corrections are very large and compensating, so the individual `hc`, `disp`, and `assoc` bars swing dramatically while the sum still closes to the correct total.\n"
        )
        handle.write(
            "- The paper Figure 3 short-range bars do not close even on the raw-`mu` basis, so there is likely a second issue as well: either digitization error in the small bars, rounding from the publication, or inconsistent internal plotting/export bookkeeping in the paper.\n"
        )

    print(f"Wrote {detail_path}")
    print(f"Wrote {term_summary_path}")
    print(f"Wrote {ion_summary_path}")
    print(f"Wrote {summary_path}")


if __name__ == "__main__":
    main()
