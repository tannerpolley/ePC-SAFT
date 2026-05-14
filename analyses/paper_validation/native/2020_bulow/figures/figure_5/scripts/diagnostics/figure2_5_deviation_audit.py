from __future__ import annotations

import csv
import math
import sys


from collections import defaultdict
from pathlib import Path
import sys as _bootstrap_sys
from pathlib import Path as _BootstrapPath

for _candidate in _BootstrapPath(__file__).resolve().parents:
    if (_candidate / "scripts" / "plot_outputs.py").is_file():
        if str(_candidate) not in _bootstrap_sys.path:
            _bootstrap_sys.path.insert(0, str(_candidate))
        break
else:
    raise ModuleNotFoundError("Could not locate repo root containing scripts/plot_outputs.py")
from scripts.plot_outputs import REPO_ROOT

import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
FIGURE_DIR = SCRIPT_DIR.parent
ANALYSIS_ROOT = FIGURE_DIR.parents[2]
OUTPUT_DIR = common.analysis_runs_path(__file__, "_placeholder", category=("figure_5", "diagnostics")).parent

if str(ANALYSIS_ROOT) not in sys.path:
    sys.path.insert(0, str(ANALYSIS_ROOT))
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

import _plot_common as common
import _model_overlay as overlay

FIG3_PAPER_ROW_KEYS = {
    "hc": ("hc avg", "hc"),
    "disp": ("disp avg", "disp"),
    "assoc": ("assoc avg", "assoc"),
    "born": ("born avg", "born"),
}


def _summary(rows: list[dict[str, object]], keys: tuple[str, ...]) -> list[dict[str, object]]:
    grouped: dict[tuple[object, ...], list[float]] = defaultdict(list)
    for row in rows:
        grouped[tuple(row[key] for key in keys)].append(float(row["delta"]))

    out: list[dict[str, object]] = []
    for group, deltas in grouped.items():
        arr = np.asarray(deltas, dtype=float)
        item = {key: value for key, value in zip(keys, group, strict=False)}
        item["count"] = int(arr.size)
        item["mean_delta"] = float(np.mean(arr))
        item["mean_abs_delta"] = float(np.mean(np.abs(arr)))
        item["rmse"] = float(np.sqrt(np.mean(arr**2)))
        item["max_abs_delta"] = float(np.max(np.abs(arr)))
        item["positive_fraction"] = float(np.mean(arr > 0.0))
        out.append(item)
    return sorted(out, key=lambda entry: tuple(entry[key] for key in keys))


def _figure3_paper_values(frame: common.IndexedCsv, term_key: str, ions: list[str]) -> np.ndarray:
    for row_key in FIG3_PAPER_ROW_KEYS[term_key]:
        if row_key in frame.index:
            return frame.values(row_key, ions)
    raise KeyError(f"Missing figure 3 paper row for {term_key}")


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    detail_rows: list[dict[str, object]] = []

    # Figure 2
    fig2 = common.load_indexed_csv(
        common.analysis_data_path(ANALYSIS_ROOT / "figure_2", "water_comparisons.csv", kind="processed", category="figure_2")
    )
    ions = list(fig2.columns)
    for variant in ("advanced", "revised"):
        paper = fig2.values(variant, ions)
        model = np.asarray([overlay.gsolv_ion(variant, ion, "water") for ion in ions], dtype=float)
        for ion, paper_val, model_val in zip(ions, paper, model, strict=False):
            detail_rows.append(
                {
                    "figure": "figure_2",
                    "group": variant,
                    "subgroup": "water",
                    "item": ion,
                    "paper": float(paper_val),
                    "model": float(model_val),
                    "delta": float(model_val - paper_val),
                }
            )

    # Figure 3
    fig3 = common.load_indexed_csv(
        common.analysis_data_path(ANALYSIS_ROOT / "figure_3", "water_contributions.csv", kind="processed", category="figure_3")
    )
    ions = list(fig3.columns)
    for term in ("hc", "disp", "assoc", "dh", "born"):
        if term == "dh":
            paper = np.zeros(len(ions), dtype=float)
        else:
            paper = _figure3_paper_values(fig3, term, ions)
        model = np.asarray(
            [overlay.contribution_breakdown("advanced", ion, "water")[term] for ion in ions], dtype=float
        )
        for ion, paper_val, model_val in zip(ions, paper, model, strict=False):
            detail_rows.append(
                {
                    "figure": "figure_3",
                    "group": term,
                    "subgroup": "water",
                    "item": ion,
                    "paper": float(paper_val),
                    "model": float(model_val),
                    "delta": float(model_val - paper_val),
                }
            )

    # Figure 4
    for solvent in ("methanol", "ethanol"):
        fig4 = common.load_indexed_csv(
            common.analysis_data_path(
                ANALYSIS_ROOT / "figure_4", f"water-{solvent}-comparison.csv", kind="processed", category="figure_4"
            )
        )
        ions = list(fig4.columns)
        for variant in ("advanced", "revised"):
            paper = fig4.values(variant, ions)
            model = np.asarray([overlay.transfer_total(variant, ion, solvent) for ion in ions], dtype=float)
            for ion, paper_val, model_val in zip(ions, paper, model, strict=False):
                detail_rows.append(
                    {
                        "figure": "figure_4",
                        "group": variant,
                        "subgroup": solvent,
                        "item": ion,
                        "paper": float(paper_val),
                        "model": float(model_val),
                        "delta": float(model_val - paper_val),
                    }
                )

    # Figure 5
    for solvent in ("methanol", "ethanol"):
        fig5 = common.load_indexed_csv(
            common.analysis_data_path(
                ANALYSIS_ROOT / "figure_5", f"water-{solvent}-contributions.csv", kind="processed", category="figure_5"
            )
        )
        for term in ("hc", "disp", "assoc", "dh", "Born"):
            for ion in ("Na+", "Cl-", "I-"):
                if term == "dh":
                    paper_val = 0.0
                    term_key = "dh"
                else:
                    paper_val = fig5.scalar(term, ion)
                    term_key = term.lower() if term != "Born" else "born"
                model_val = overlay.transfer_breakdown("advanced", ion, solvent)[term_key]
                detail_rows.append(
                    {
                        "figure": "figure_5",
                        "group": term_key,
                        "subgroup": solvent,
                        "item": ion,
                        "paper": float(paper_val),
                        "model": float(model_val),
                        "delta": float(model_val - paper_val),
                    }
                )

    detail_path = OUTPUT_DIR / "figure2_5_deviation_details.csv"
    with detail_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["figure", "group", "subgroup", "item", "paper", "model", "delta"])
        writer.writeheader()
        writer.writerows(detail_rows)

    summary_specs = [
        ("Figure 2", ("figure", "group")),
        ("Figure 3", ("figure", "group")),
        ("Figure 4", ("figure", "subgroup", "group")),
        ("Figure 5", ("figure", "subgroup", "group")),
    ]

    summary_path = OUTPUT_DIR / "figure2_5_deviation_summary.md"
    with summary_path.open("w", encoding="utf-8") as handle:
        handle.write("# Figure 2-5 Deviation Audit\n\n")
        handle.write("All deltas are `epcsaft model - paper value` in kJ/mol.\n\n")
        for title, keys in summary_specs:
            handle.write(f"## {title}\n\n")
            rows = _summary(detail_rows, keys)
            header = " | ".join(
                [*keys, "count", "mean_delta", "mean_abs_delta", "rmse", "max_abs_delta", "positive_fraction"]
            )
            sep = " | ".join(["---"] * (len(keys) + 5))
            handle.write(f"| {header} |\n")
            handle.write(f"| {sep} |\n")
            for row in rows:
                values = [row[key] for key in keys]
                values.extend(
                    [
                        row["count"],
                        f"{row['mean_delta']:.3f}",
                        f"{row['mean_abs_delta']:.3f}",
                        f"{row['rmse']:.3f}",
                        f"{row['max_abs_delta']:.3f}",
                        f"{row['positive_fraction']:.3f}",
                    ]
                )
                handle.write("| " + " | ".join(str(v) for v in values) + " |\n")
            handle.write("\n")

        handle.write("## Largest Absolute Deviations\n\n")
        worst = sorted(detail_rows, key=lambda row: abs(float(row["delta"])), reverse=True)[:15]
        handle.write("| figure | group | subgroup | item | paper | model | delta |\n")
        handle.write("| --- | --- | --- | --- | --- | --- | --- |\n")
        for row in worst:
            handle.write(
                f"| {row['figure']} | {row['group']} | {row['subgroup']} | {row['item']} | "
                f"{float(row['paper']):.3f} | {float(row['model']):.3f} | {float(row['delta']):.3f} |\n"
            )

    print(f"Wrote {detail_path}")
    print(f"Wrote {summary_path}")


if __name__ == "__main__":
    main()

