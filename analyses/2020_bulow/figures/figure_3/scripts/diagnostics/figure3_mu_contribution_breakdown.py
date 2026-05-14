from __future__ import annotations

import csv
import sys
from pathlib import Path

import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
FIGURE_DIR = SCRIPT_DIR.parent
ANALYSIS_ROOT = FIGURE_DIR.parents[2]
REPO_ROOT = ANALYSIS_ROOT.parents[1]

if str(ANALYSIS_ROOT) not in sys.path:
    sys.path.insert(0, str(ANALYSIS_ROOT))
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

import _plot_common as common
import _model_overlay as overlay

DATA_PATH = common.analysis_data_path(FIGURE_DIR, "water_contributions.csv", kind="processed", category="figure_3")
OUTPUT_CSV = common.analysis_runs_path(__file__, "figure3_mu_contribution_breakdown.csv", category=("figure_3", "diagnostics"))
R_GAS = overlay.R_GAS
T_REF = overlay.T_REF
RT_KJMOL = R_GAS * T_REF / 1000.0

CONTRIBUTION_MAP = {
    "hc": {"paper_rows": ("hc avg", "hc"), "suffix": "hc"},
    "disp": {"paper_rows": ("disp avg", "disp"), "suffix": "disp"},
    "assoc": {"paper_rows": ("assoc avg", "assoc"), "suffix": "assoc"},
    "dh": {"paper_rows": (), "suffix": "ion"},
    "born": {"paper_rows": ("born avg", "born"), "suffix": "born"},
}


def _paper_value(frame: common.Table, contribution: str, ion: str) -> float:
    if contribution == "dh":
        return 0.0
    for row_key in CONTRIBUTION_MAP[contribution]["paper_rows"]:
        if row_key in frame.index:
            return frame.scalar(row_key, ion)
    raise KeyError(f"Missing paper row for contribution '{contribution}'.")


def _kjmol(value: float) -> float:
    return float(RT_KJMOL * float(value))


def _build_rows() -> list[dict[str, object]]:
    frame = common.load_indexed_csv(DATA_PATH)
    ions = list(frame.columns)
    rows: list[dict[str, object]] = []
    contribution_order = {key: idx for idx, key in enumerate(CONTRIBUTION_MAP)}
    ion_order = {ion: idx for idx, ion in enumerate(ions)}

    for contribution, meta in CONTRIBUTION_MAP.items():
        for ion in ions:
            terms, idx = overlay._infinite_dilution_terms("advanced", ion, "water")
            suffix = str(meta["suffix"])
            paper_mu = _paper_value(frame, contribution, ion)
            mu_term = float(np.asarray(terms[f"mu_{suffix}"], dtype=float)[idx])
            a_term = float(terms[f"a_{suffix}"])
            z_term = float(terms[f"z_{suffix}"])
            dadx_term = float(np.asarray(terms[f"dadx_{suffix}"], dtype=float)[idx])
            sum_x_dadx_term = float(terms[f"sum_x_dadx_{suffix}"])

            epcsaft_mu = _kjmol(mu_term)
            a_kjmol = _kjmol(a_term)
            z_kjmol = _kjmol(z_term)
            dadx_kjmol = _kjmol(dadx_term)
            sum_x_dadx_kjmol = _kjmol(sum_x_dadx_term)
            rebuilt_mu = a_kjmol + z_kjmol + dadx_kjmol - sum_x_dadx_kjmol

            rows.append(
                {
                    "ion": ion,
                    "contr": contribution,
                    "paper_mu_contr": paper_mu,
                    "epcsaft_mu_contr": epcsaft_mu,
                    "epcsaft_mu_manual_sum": rebuilt_mu,
                    "a_contr": a_kjmol,
                    "z_contr": z_kjmol,
                    "dadx_contr": dadx_kjmol,
                    "sum_xj_dadx_contr": sum_x_dadx_kjmol,
                    "manual_minus_epcsaft": rebuilt_mu - epcsaft_mu,
                    "paper_minus_epcsaft": paper_mu - epcsaft_mu,
                }
            )

    rows.sort(key=lambda row: (contribution_order[str(row["contr"])], ion_order[str(row["ion"])]))
    return rows


def main() -> None:
    rows = _build_rows()
    OUTPUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    with OUTPUT_CSV.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "ion",
                "contr",
                "paper_mu_contr",
                "epcsaft_mu_contr",
                "epcsaft_mu_manual_sum",
                "a_contr",
                "z_contr",
                "dadx_contr",
                "sum_xj_dadx_contr",
                "manual_minus_epcsaft",
                "paper_minus_epcsaft",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {OUTPUT_CSV}", flush=True)


if __name__ == "__main__":
    main()

