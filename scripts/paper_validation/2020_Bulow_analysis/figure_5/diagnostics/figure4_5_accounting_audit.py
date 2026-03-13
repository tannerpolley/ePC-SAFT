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
FIG4_PATHS = {
    "methanol": ANALYSIS_ROOT / "figure_4" / "data" / "water-methanol-comparison.csv",
    "ethanol": ANALYSIS_ROOT / "figure_4" / "data" / "water-ethanol-comparison.csv",
}
FIG5_PATHS = {
    "methanol": ANALYSIS_ROOT / "figure_5" / "data" / "water-methanol-contributions.csv",
    "ethanol": ANALYSIS_ROOT / "figure_5" / "data" / "water-ethanol-contributions.csv",
}
IONS = ("Na+", "Cl-", "I-")
TERMS = ("hc", "disp", "assoc", "born")


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
    key = term
    terms, idx = overlay._infinite_dilution_terms("advanced", ion, solvent)
    mu = float(R_GAS * T_REF * np.asarray(terms[f"mu_{key}"], dtype=float)[idx] / 1000.0)
    lnfug = float(R_GAS * T_REF * np.asarray(terms[f"lnfugcoef_{key}"], dtype=float)[idx] / 1000.0)
    z_share = float(mu - lnfug)
    return {
        "mu": mu,
        "lnfug": lnfug,
        "zcorr_share": z_share,
        "z_term": float(terms[f"z_{key}"]),
        "z_total": float(terms["z_total"]),
    }


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    detail_rows: list[dict[str, object]] = []
    summary_rows: list[dict[str, object]] = []
    total_rows: list[dict[str, object]] = []

    for solvent in ("methanol", "ethanol"):
        fig4 = common.load_indexed_csv(FIG4_PATHS[solvent])
        fig5 = common.load_indexed_csv(FIG5_PATHS[solvent])

        for term in TERMS:
            paper = np.asarray(
                [float(fig5.scalar("Born" if term == "born" else term, ion)) for ion in IONS],
                dtype=float,
            )
            delta_mu = np.empty(len(IONS), dtype=float)
            delta_lnfug = np.empty(len(IONS), dtype=float)
            delta_zcorr = np.empty(len(IONS), dtype=float)

            for idx, ion in enumerate(IONS):
                organic = _term_payload(ion, solvent, term)
                water = _term_payload(ion, "water", term)

                delta_mu[idx] = organic["mu"] - water["mu"]
                delta_lnfug[idx] = organic["lnfug"] - water["lnfug"]
                delta_zcorr[idx] = organic["zcorr_share"] - water["zcorr_share"]

                detail_rows.append(
                    {
                        "solvent": solvent,
                        "ion": ion,
                        "term": term,
                        "paper": float(paper[idx]),
                        "transfer_mu": float(delta_mu[idx]),
                        "transfer_lnfug": float(delta_lnfug[idx]),
                        "transfer_zcorr_share": float(delta_zcorr[idx]),
                        "water_mu": water["mu"],
                        "organic_mu": organic["mu"],
                        "water_lnfug": water["lnfug"],
                        "organic_lnfug": organic["lnfug"],
                        "water_zcorr_share": water["zcorr_share"],
                        "organic_zcorr_share": organic["zcorr_share"],
                        "water_z_term": water["z_term"],
                        "organic_z_term": organic["z_term"],
                        "water_z_total": water["z_total"],
                        "organic_z_total": organic["z_total"],
                        "paper_sign": int(np.sign(paper[idx])),
                        "mu_sign": int(np.sign(delta_mu[idx])),
                        "lnfug_sign": int(np.sign(delta_lnfug[idx])),
                    }
                )

            for candidate_name, candidate in {
                "mu": delta_mu,
                "lnfug": delta_lnfug,
                "zcorr_share": delta_zcorr,
                "neg_mu": -delta_mu,
                "neg_lnfug": -delta_lnfug,
            }.items():
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

        for ion in IONS:
            paper_total = float(fig4.scalar("advanced", ion))
            paper_short = float(sum(float(fig5.scalar(key, ion)) for key in ("hc", "disp", "assoc")))
            paper_born = float(fig5.scalar("Born", ion))
            model_transfer = overlay.transfer_breakdown("advanced", ion, solvent)
            model_short_mu = sum(
                _term_payload(ion, solvent, key)["mu"] - _term_payload(ion, "water", key)["mu"]
                for key in ("hc", "disp", "assoc")
            )
            total_rows.append(
                {
                    "solvent": solvent,
                    "ion": ion,
                    "paper_total_figure4": paper_total,
                    "paper_short_sum_figure5": paper_short,
                    "paper_born_figure5": paper_born,
                    "paper_total_minus_born": paper_total - paper_born,
                    "model_total_lnfug": float(overlay.transfer_total("advanced", ion, solvent)),
                    "model_short_sum_lnfug": float(sum(model_transfer[key] for key in ("hc", "disp", "assoc"))),
                    "model_short_sum_mu": float(model_short_mu),
                    "model_born_lnfug": float(model_transfer["born"]),
                }
            )

    detail_path = OUTPUT_DIR / "figure4_5_accounting_detail.csv"
    with detail_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "solvent",
                "ion",
                "term",
                "paper",
                "transfer_mu",
                "transfer_lnfug",
                "transfer_zcorr_share",
                "water_mu",
                "organic_mu",
                "water_lnfug",
                "organic_lnfug",
                "water_zcorr_share",
                "organic_zcorr_share",
                "water_z_term",
                "organic_z_term",
                "water_z_total",
                "organic_z_total",
                "paper_sign",
                "mu_sign",
                "lnfug_sign",
            ],
        )
        writer.writeheader()
        writer.writerows(detail_rows)

    summary_path = OUTPUT_DIR / "figure4_5_accounting_summary.csv"
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

    total_path = OUTPUT_DIR / "figure4_5_total_crosscheck.csv"
    with total_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "solvent",
                "ion",
                "paper_total_figure4",
                "paper_short_sum_figure5",
                "paper_born_figure5",
                "paper_total_minus_born",
                "model_total_lnfug",
                "model_short_sum_lnfug",
                "model_short_sum_mu",
                "model_born_lnfug",
            ],
        )
        writer.writeheader()
        writer.writerows(total_rows)

    best_rows: dict[tuple[str, str], dict[str, object]] = {}
    for row in sorted(summary_rows, key=lambda entry: (entry["solvent"], entry["term"], entry["rmse"])):
        best_rows.setdefault((str(row["solvent"]), str(row["term"])), row)

    report_path = OUTPUT_DIR / "figure4_5_accounting_audit.md"
    with report_path.open("w", encoding="utf-8") as handle:
        handle.write("# Figure 4-5 Accounting Audit\n\n")
        handle.write("This audit compares the paper Figure 5 transfer contributions against model-side transfer values formed from exposed `mu`, `lnfug`, and `Z`-share terms at the same infinite-dilution states used for Figure 4.\n\n")

        handle.write("## C++ Contribution Construction\n\n")
        handle.write("- In [pcsaft_electrolyte.cpp](C:/Users/Tanner/Documents/git/PC-SAFT/pcsaft_electrolyte.cpp), the hard-chain and dispersion branches are built as\n")
        handle.write("  `mu_hc[i] = ares_hc + Zhc + dahc_dx[i] - sum_j x[j]*dahc_dx[j]`\n")
        handle.write("  `mu_disp[i] = ares_disp + Zdisp + dadisp_dx[i] - sum_j x[j]*dadisp_dx[j]`\n")
        handle.write("- The Debye-Huckel and Born branches use the same pattern with `a_DH` / `a_born`, `Z_DH` / `Zborn`, and the corresponding `dadx` vectors.\n")
        handle.write("- The current public API already exposes the downstream `mu_*`, `lnfugcoef_*`, and `z_*` terms, but it does **not** yet expose the internal `a^alpha`, `dadx^alpha`, or `sum_x dadx^alpha` pieces separately. This audit therefore uses the fully exposed downstream pieces first.\n\n")

        handle.write("## Main Findings\n\n")
        for solvent in ("methanol", "ethanol"):
            handle.write(f"### {solvent.capitalize()}\n\n")
            for term in TERMS:
                best = best_rows[(solvent, term)]
                mu_row = next(row for row in summary_rows if row["solvent"] == solvent and row["term"] == term and row["candidate"] == "mu")
                lnf_row = next(row for row in summary_rows if row["solvent"] == solvent and row["term"] == term and row["candidate"] == "lnfug")
                handle.write(
                    f"- `{term}`: best exposed match is `{best['candidate']}` with RMSE `{best['rmse']:.3f}` kJ/mol; `mu` RMSE is `{mu_row['rmse']:.3f}` and `lnfug` RMSE is `{lnf_row['rmse']:.3f}`.\n"
                )
            handle.write("\n")

        handle.write("## Candidate Ranking By Solvent and Term\n\n")
        handle.write("| solvent | term | best candidate | RMSE | MAE | corr | sign matches |\n")
        handle.write("| --- | --- | --- | --- | --- | --- | --- |\n")
        for solvent in ("methanol", "ethanol"):
            for term in TERMS:
                row = best_rows[(solvent, term)]
                handle.write(
                    f"| {solvent} | {term} | {row['candidate']} | {row['rmse']:.3f} | {row['mae']:.3f} | {row['corr']:.3f} | {row['sign_matches']} |\n"
                )
        handle.write("\n")

        handle.write("## Figure 4 vs Figure 5 Cross-Check\n\n")
        handle.write("| solvent | ion | paper total (Fig 4) | paper total - born | model short sum lnfug | model short sum mu | model born lnfug |\n")
        handle.write("| --- | --- | --- | --- | --- | --- | --- |\n")
        for row in total_rows:
            handle.write(
                f"| {row['solvent']} | {row['ion']} | {row['paper_total_figure4']:.3f} | {row['paper_total_minus_born']:.3f} | {row['model_short_sum_lnfug']:.3f} | {row['model_short_sum_mu']:.3f} | {row['model_born_lnfug']:.3f} |\n"
            )
        handle.write("\n")

        handle.write("## Interpretation\n\n")
        handle.write("- For both alcohols, the paper `hc`, `disp`, and `assoc` bars are much closer to the transfer in raw $\\Delta\\tilde{\\mu}^\\alpha$ than to the transfer in $\\Delta\\ln\\varphi^\\alpha$.\n")
        handle.write("- The Born transfer bar is different: it does **not** match the model well even on the raw-`mu` route, which means Figure 5 is not explained by a single bookkeeping switch the way Figure 3 largely is.\n")
        handle.write("- Ethanol is the clearest case: the paper short-range bars are close to the raw `mu` transfer values, while the model totals only close when the correct $\\ln\\varphi$ transfer contributions are used. That is why Figure 4 totals can agree while Figure 5 component bars disagree strongly, including sign flips.\n")
        handle.write("- The detailed CSV contains the water-side and organic-side state values for each branch, so you can inspect exactly where each sign change enters.\n")

    print(f"Wrote {detail_path}")
    print(f"Wrote {summary_path}")
    print(f"Wrote {total_path}")
    print(f"Wrote {report_path}")


if __name__ == "__main__":
    main()
