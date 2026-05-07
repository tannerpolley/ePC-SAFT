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
OUTPUT_DIR = common.analysis_runs_path(__file__, "_placeholder", category=("figure_5", "diagnostics")).parent

if str(ANALYSIS_ROOT) not in sys.path:
    sys.path.insert(0, str(ANALYSIS_ROOT))
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

import _model_overlay as overlay
import _plot_common as common

R_GAS = 8.31446261815324
T_REF = 298.15
FIG5_PATHS = {
    "methanol": common.analysis_data_path(
        ANALYSIS_ROOT / "figure_5", "water-methanol-contributions.csv", kind="processed", category="figure_5"
    ),
    "ethanol": common.analysis_data_path(
        ANALYSIS_ROOT / "figure_5", "water-ethanol-contributions.csv", kind="processed", category="figure_5"
    ),
}
IONS = ("Na+", "Cl-", "I-")
TERM_TO_KEY = {
    "hc": "hc",
    "disp": "disp",
    "assoc": "assoc",
    "dh": "ion",
    "born": "born",
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


def _term_payload(ion: str, solvent: str, term: str) -> dict[str, float]:
    key = TERM_TO_KEY[term]
    terms, idx = overlay._infinite_dilution_terms("advanced", ion, solvent)
    rt = R_GAS * T_REF / 1000.0
    return {
        "a": float(rt * float(terms[f"a_{key}"])),
        "z_raw": float(rt * float(terms[f"z_raw_{key}"])),
        "z_norm": float(float(terms[f"z_{key}"])),
        "dadx": float(rt * np.asarray(terms[f"dadx_{key}"], dtype=float)[idx]),
        "sum_x_dadx": float(rt * float(terms[f"sum_x_dadx_{key}"])),
        "mu": float(rt * np.asarray(terms[f"mu_{key}"], dtype=float)[idx]),
        "lnfug": float(rt * np.asarray(terms[f"lnfugcoef_{key}"], dtype=float)[idx]),
        "z_share": float(
            rt
            * (
                np.asarray(terms[f"mu_{key}"], dtype=float)[idx]
                - np.asarray(terms[f"lnfugcoef_{key}"], dtype=float)[idx]
            )
        ),
        "z_total": float(terms["z_total"]),
    }


def _paper_value(fig5_data, term: str, ion: str) -> float:
    if term == "born":
        return float(fig5_data.scalar("Born", ion))
    if term == "dh":
        return 0.0
    return float(fig5_data.scalar(term, ion))


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    detail_rows: list[dict[str, object]] = []
    summary_rows: list[dict[str, object]] = []

    for solvent in ("ethanol", "methanol"):
        fig5 = common.load_indexed_csv(FIG5_PATHS[solvent])

        for term in ("hc", "disp", "assoc", "dh", "born"):
            paper = np.empty(len(IONS), dtype=float)
            candidates: dict[str, np.ndarray] = {
                "a": np.empty(len(IONS), dtype=float),
                "z_raw": np.empty(len(IONS), dtype=float),
                "dadx": np.empty(len(IONS), dtype=float),
                "neg_sum_x_dadx": np.empty(len(IONS), dtype=float),
                "a_plus_z_raw": np.empty(len(IONS), dtype=float),
                "dadx_minus_sum": np.empty(len(IONS), dtype=float),
                "a_plus_dadx_minus_sum": np.empty(len(IONS), dtype=float),
                "reconstructed_mu": np.empty(len(IONS), dtype=float),
                "mu": np.empty(len(IONS), dtype=float),
                "lnfug": np.empty(len(IONS), dtype=float),
                "z_share": np.empty(len(IONS), dtype=float),
            }

            for idx, ion in enumerate(IONS):
                paper[idx] = _paper_value(fig5, term, ion)
                organic = _term_payload(ion, solvent, term)
                water = _term_payload(ion, "water", term)

                delta = {key: organic[key] - water[key] for key in organic}
                delta["neg_sum_x_dadx"] = -delta["sum_x_dadx"]
                delta["a_plus_z_raw"] = delta["a"] + delta["z_raw"]
                delta["dadx_minus_sum"] = delta["dadx"] - delta["sum_x_dadx"]
                delta["a_plus_dadx_minus_sum"] = delta["a"] + delta["dadx"] - delta["sum_x_dadx"]
                delta["reconstructed_mu"] = delta["a"] + delta["z_raw"] + delta["dadx"] - delta["sum_x_dadx"]

                for name in candidates:
                    candidates[name][idx] = float(delta[name])

                detail_rows.append(
                    {
                        "solvent": solvent,
                        "ion": ion,
                        "term": term,
                        "paper": paper[idx],
                        "water_a": water["a"],
                        "organic_a": organic["a"],
                        "delta_a": delta["a"],
                        "water_z_raw": water["z_raw"],
                        "organic_z_raw": organic["z_raw"],
                        "delta_z_raw": delta["z_raw"],
                        "water_dadx": water["dadx"],
                        "organic_dadx": organic["dadx"],
                        "delta_dadx": delta["dadx"],
                        "water_sum_x_dadx": water["sum_x_dadx"],
                        "organic_sum_x_dadx": organic["sum_x_dadx"],
                        "delta_sum_x_dadx": delta["sum_x_dadx"],
                        "delta_neg_sum_x_dadx": delta["neg_sum_x_dadx"],
                        "delta_a_plus_z_raw": delta["a_plus_z_raw"],
                        "delta_dadx_minus_sum": delta["dadx_minus_sum"],
                        "delta_a_plus_dadx_minus_sum": delta["a_plus_dadx_minus_sum"],
                        "water_mu": water["mu"],
                        "organic_mu": organic["mu"],
                        "delta_mu": delta["mu"],
                        "water_lnfug": water["lnfug"],
                        "organic_lnfug": organic["lnfug"],
                        "delta_lnfug": delta["lnfug"],
                        "delta_reconstructed_mu": delta["reconstructed_mu"],
                        "water_z_share": water["z_share"],
                        "organic_z_share": organic["z_share"],
                        "delta_z_share": delta["z_share"],
                        "water_z_norm": water["z_norm"],
                        "organic_z_norm": organic["z_norm"],
                        "water_z_total": water["z_total"],
                        "organic_z_total": organic["z_total"],
                    }
                )

            for candidate_name, candidate in candidates.items():
                slope, intercept, fit_rmse = _fit_line(candidate, paper)
                summary_rows.append(
                    {
                        "solvent": solvent,
                        "term": term,
                        "candidate": candidate_name,
                        "rmse": float(np.sqrt(np.mean((candidate - paper) ** 2))),
                        "mae": float(np.mean(np.abs(candidate - paper))),
                        "corr": _safe_corr(candidate, paper),
                        "fit_slope": slope,
                        "fit_intercept": intercept,
                        "fit_rmse": fit_rmse,
                        "sign_matches": int(np.sum(np.sign(candidate) == np.sign(paper))),
                    }
                )

    detail_path = OUTPUT_DIR / "figure5_internal_term_detail.csv"
    with detail_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "solvent",
                "ion",
                "term",
                "paper",
                "water_a",
                "organic_a",
                "delta_a",
                "water_z_raw",
                "organic_z_raw",
                "delta_z_raw",
                "water_dadx",
                "organic_dadx",
                "delta_dadx",
                "water_sum_x_dadx",
                "organic_sum_x_dadx",
                "delta_sum_x_dadx",
                "delta_neg_sum_x_dadx",
                "delta_a_plus_z_raw",
                "delta_dadx_minus_sum",
                "delta_a_plus_dadx_minus_sum",
                "water_mu",
                "organic_mu",
                "delta_mu",
                "water_lnfug",
                "organic_lnfug",
                "delta_lnfug",
                "delta_reconstructed_mu",
                "water_z_share",
                "organic_z_share",
                "delta_z_share",
                "water_z_norm",
                "organic_z_norm",
                "water_z_total",
                "organic_z_total",
            ],
        )
        writer.writeheader()
        writer.writerows(detail_rows)

    summary_path = OUTPUT_DIR / "figure5_internal_term_summary.csv"
    with summary_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "solvent",
                "term",
                "candidate",
                "rmse",
                "mae",
                "corr",
                "fit_slope",
                "fit_intercept",
                "fit_rmse",
                "sign_matches",
            ],
        )
        writer.writeheader()
        writer.writerows(sorted(summary_rows, key=lambda row: (row["solvent"], row["term"], row["rmse"])))

    best_rows: dict[tuple[str, str], dict[str, object]] = {}
    for row in sorted(summary_rows, key=lambda entry: (entry["solvent"], entry["term"], entry["rmse"])):
        best_rows.setdefault((str(row["solvent"]), str(row["term"])), row)

    report_path = OUTPUT_DIR / "figure5_internal_term_audit.md"
    with report_path.open("w", encoding="utf-8") as handle:
        handle.write("# Figure 5 Internal Term Audit\n\n")
        handle.write(
            "This audit compares the paper Figure 5 transfer bars against the model-side transfer deltas of the internal contribution pieces exposed from `epcsaft_fugacity_coefficient_terms(...)`.\n\n"
        )
        handle.write(
            "For each contribution $\\alpha$, the report tracks the state-level pieces used in the code path:\n\n"
        )
        handle.write("- $a^\\alpha$\n")
        handle.write("- $Z^\\alpha_{\\mathrm{raw}}$\n")
        handle.write("- $(\\partial a^\\alpha/\\partial x_i)$ at the ion index\n")
        handle.write("- $\\sum_j x_j (\\partial a^\\alpha/\\partial x_j)$\n")
        handle.write("- the final code-level $\\mu_i^\\alpha$\n")
        handle.write("- the final $\\ln\\varphi_i^\\alpha$\n\n")
        handle.write(
            "Transfer deltas are formed as organic minus water at the same infinite-dilution reference states used for Figure 4 and Figure 5.\n\n"
        )

        handle.write("## Best-Matching Candidate By Solvent and Contribution\n\n")
        for solvent in ("ethanol", "methanol"):
            handle.write(f"### {solvent.title()}\n\n")
            for term in ("hc", "disp", "assoc", "dh", "born"):
                row = best_rows[(solvent, term)]
                handle.write(
                    f"- `{term}`: best candidate is `{row['candidate']}` with RMSE `{row['rmse']:.3f} kJ/mol`, "
                    f"MAE `{row['mae']:.3f} kJ/mol`, sign matches `{int(row['sign_matches'])}/3`.\n"
                )
            handle.write("\n")

        handle.write("## Notes\n\n")
        handle.write(
            "- `reconstructed_mu` is $\\Delta a^\\alpha + \\Delta Z^\\alpha_{raw} + \\Delta(\\partial a^\\alpha/\\partial x_i) - \\Delta\\sum_j x_j (\\partial a^\\alpha/\\partial x_j)$.\n"
        )
        handle.write(
            "- For `hc`, `disp`, `polar`, `dh`, and `born`, `reconstructed_mu` should match the exposed `mu` branch to numerical precision.\n"
        )
        handle.write(
            "- The association branch is dumped exactly as the code currently assembles it. In the current implementation, the exposed `dadx_assoc` column is the same direct code-level quantity as `mu_assoc`, not a separate post-processed correction path.\n"
        )
        handle.write(
            "- The detailed per-ion dump is in `figure5_internal_term_detail.csv`, and the candidate ranking table is in `figure5_internal_term_summary.csv`.\n"
        )


if __name__ == "__main__":
    main()
