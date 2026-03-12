from __future__ import annotations

import csv
import math
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
ANALYSIS_ROOT = SCRIPT_DIR.parent
REPO_ROOT = ANALYSIS_ROOT.parents[2]
OUTPUT_DIR = SCRIPT_DIR / "output"

if str(ANALYSIS_ROOT) not in sys.path:
    sys.path.insert(0, str(ANALYSIS_ROOT))
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

import _plot_common as common
from data.epcsaft_properties import get_prop_dict
from pcsaft import pcsaft_den, pcsaft_gsolv, pcsaft_lnfugcoef_terms


R_GAS = 8.31446261815324
T_REF = 298.15
P_REF = 1.0e5
EPS = 1.0e-8
SOLVENT_SPECIES = {"water": "Water", "methanol": "Methanol", "ethanol": "Ethanol"}
VARIANT_DATASET = {"advanced": "bulow_2020", "revised": "held_2014"}


def _species_for_ion(ion: str, solvent: str) -> list[str]:
    solvent_species = SOLVENT_SPECIES[solvent]
    if ion in {"Li+", "Na+", "K+"}:
        return [ion, "Cl-", solvent_species]
    if ion == "F-":
        return ["Na+", "F-", solvent_species]
    return ["Na+", ion, solvent_species]


def _gsolv_ion(variant: str, ion: str, solvent: str) -> float:
    species = _species_for_ion(ion, solvent)
    x = np.asarray([EPS, EPS, 1.0 - 2.0 * EPS], dtype=float)
    params = get_prop_dict(VARIANT_DATASET[variant], species, x, T_REF, user_options={})
    rho = pcsaft_den(T_REF, P_REF, x, params, phase="liq")
    values = pcsaft_gsolv(T_REF, rho, x, params, species=species)
    return float(values[ion]) / 1000.0


def _term_breakdown_advanced(ion: str, solvent: str) -> dict[str, float]:
    species = _species_for_ion(ion, solvent)
    x = np.asarray([EPS, EPS, 1.0 - 2.0 * EPS], dtype=float)
    params = get_prop_dict("bulow_2020", species, x, T_REF, user_options={})
    rho = pcsaft_den(T_REF, P_REF, x, params, phase="liq")
    terms = pcsaft_lnfugcoef_terms(T_REF, rho, x, params)
    idx = species.index(ion)
    return {
        "hc": float(R_GAS * T_REF * terms["mu_hc"][idx] / 1000.0),
        "disp": float(R_GAS * T_REF * terms["mu_disp"][idx] / 1000.0),
        "assoc": float(R_GAS * T_REF * terms["mu_assoc"][idx] / 1000.0),
        "born": float(R_GAS * T_REF * terms["mu_born"][idx] / 1000.0),
        "total": float(R_GAS * T_REF * terms["lnfugcoef_total"][idx] / 1000.0),
    }


def _transfer_total(variant: str, ion: str, solvent: str) -> float:
    return _gsolv_ion(variant, ion, solvent) - _gsolv_ion(variant, ion, "water")


def _transfer_breakdown_advanced(ion: str, solvent: str) -> dict[str, float]:
    organic = _term_breakdown_advanced(ion, solvent)
    water = _term_breakdown_advanced(ion, "water")
    return {key: organic[key] - water[key] for key in organic}


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


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    detail_rows: list[dict[str, object]] = []

    # Figure 2
    fig2 = common.load_indexed_csv(ANALYSIS_ROOT / "figure_2" / "data" / "water_comparisons.csv")
    ions = list(fig2.columns)
    for variant in ("advanced", "revised"):
        paper = fig2.values(variant, ions)
        model = np.asarray([_gsolv_ion(variant, ion, "water") for ion in ions], dtype=float)
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
    fig3 = common.load_indexed_csv(ANALYSIS_ROOT / "figure_3" / "data" / "water_contributions.csv")
    ions = list(fig3.columns)
    for term in ("hc", "disp", "assoc", "born"):
        paper = fig3.values(term, ions)
        model = np.asarray([_term_breakdown_advanced(ion, "water")[term] for ion in ions], dtype=float)
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
        fig4 = common.load_indexed_csv(ANALYSIS_ROOT / "figure_4" / "data" / f"water-{solvent}-comparison.csv")
        ions = list(fig4.columns)
        for variant in ("advanced", "revised"):
            paper = fig4.values(variant, ions)
            model = np.asarray([_transfer_total(variant, ion, solvent) for ion in ions], dtype=float)
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
        fig5 = common.load_indexed_csv(ANALYSIS_ROOT / "figure_5" / "data" / f"water-{solvent}-contributions.csv")
        for term in ("hc", "disp", "assoc", "Born"):
            for ion in ("Na+", "Cl-", "I-"):
                paper_val = fig5.scalar(term, ion)
                model_val = _transfer_breakdown_advanced(ion, solvent)[term.lower() if term != "Born" else "born"]
                detail_rows.append(
                    {
                        "figure": "figure_5",
                        "group": term.lower(),
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
        handle.write("All deltas are `pcsaft model - paper value` in kJ/mol.\n\n")
        for title, keys in summary_specs:
            handle.write(f"## {title}\n\n")
            rows = _summary(detail_rows, keys)
            header = " | ".join([*keys, "count", "mean_delta", "mean_abs_delta", "rmse", "max_abs_delta", "positive_fraction"])
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
