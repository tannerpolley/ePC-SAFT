from __future__ import annotations

import csv
import json
import subprocess
import sys
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch


SCRIPT_DIR = Path(__file__).resolve().parent
FIGURE3_DIR = SCRIPT_DIR.parent
ANALYSIS_ROOT = FIGURE3_DIR.parent
REPO_ROOT = ANALYSIS_ROOT.parents[2]

if str(ANALYSIS_ROOT) not in sys.path:
    sys.path.insert(0, str(ANALYSIS_ROOT))
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import _model_overlay as overlay
import _plot_common as common
import feos_extractor


DATA_PATH = FIGURE3_DIR / "data" / "water_contributions.csv"
FIGURE2_TOTALS_PATH = ANALYSIS_ROOT / "figure_2" / "data" / "water_comparisons.csv"
FEOS_RAW_PATH = SCRIPT_DIR / "feos_raw.json"
CLAPEYRON_RAW_PATH = SCRIPT_DIR / "clapeyron_raw.json"
RAW_LONG_PATH = SCRIPT_DIR / "raw_comparison_long.csv"
TOTAL_CHECK_PATH = SCRIPT_DIR / "total_check.csv"
SUMMARY_PATH = SCRIPT_DIR / "comparison_summary.csv"

CLAPEYRON_ROOT = Path(r"C:\Users\Tanner\Documents\git\Clapeyron.jl")
CLAPEYRON_EXTRACTOR = SCRIPT_DIR / "extract_clapeyron_figure3.jl"

TERMS = [
    ("hc", "Hard chain", "#9F9F9F"),
    ("disp", "Dispersion", "#5F5F5F"),
    ("assoc", "Association", "#111111"),
    ("dh", "Debye-Huckel", "#8C564B"),
    ("born", "Born", "#D8891C"),
]
TERM_KEYS = tuple(term for term, _, _ in TERMS)
METRICS = (*TERM_KEYS, "mu_sum", "adjusted_sum", "total")

SOURCE_ORDER = ("paper", "epcsaft", "feos", "clapeyron")
SOURCE_LABELS = {
    "paper": "Paper",
    "epcsaft": "Current repo ePC-SAFT",
    "feos": "feos",
    "clapeyron": "Clapeyron.jl",
}
SOURCE_HATCH = {
    "paper": None,
    "epcsaft": "////",
    "feos": "\\\\\\\\",
    "clapeyron": "xx",
}
SOURCE_ALPHA = {
    "paper": 1.0,
    "epcsaft": 0.82,
    "feos": 0.72,
    "clapeyron": 0.72,
}
TOTAL_SOURCE_COLORS = {
    "paper": "#7F7F7F",
    "epcsaft": "#2CA02C",
    "feos": "#1F77B4",
    "clapeyron": "#FF7F0E",
}
METRIC_BASIS = {
    "hc": "mu",
    "disp": "mu",
    "assoc": "mu",
    "dh": "mu",
    "born": "mu",
    "mu_sum": "mu",
    "adjusted_sum": "lnfug",
    "total": "lnfug",
}
PAPER_ROW_KEYS = {
    "hc": ("hc avg", "hc"),
    "disp": ("disp avg", "disp"),
    "assoc": ("assoc avg", "assoc"),
    "born": ("born avg", "born"),
}
PAPER_RANGE_ROW_KEYS = {
    "hc": ("hc low", "hc avg", "hc hi"),
    "disp": ("disp low", "disp avg", "disp hi"),
    "assoc": ("assoc low", "assoc avg", "assoc hi"),
    "born": ("born low", "born avg", "born hi"),
}


def _paper_values(frame: common.Table, term_key: str, ions: list[str]) -> np.ndarray:
    for row_key in PAPER_ROW_KEYS.get(term_key, (term_key,)):
        if row_key in frame.index:
            return frame.values(row_key, ions)
    raise KeyError(f"Missing paper row for Figure 3 term {term_key!r}.")


def _paper_range(frame: common.Table, term_key: str, ions: list[str]) -> np.ndarray | None:
    row_keys = PAPER_RANGE_ROW_KEYS.get(term_key)
    if row_keys is None or not all(row_key in frame.index for row_key in row_keys):
        return None
    low = frame.values(row_keys[0], ions)
    avg = frame.values(row_keys[1], ions)
    high = frame.values(row_keys[2], ions)
    return np.vstack([np.abs(avg - low), np.abs(high - avg)])


def _load_paper_data(frame: common.Table, totals: common.Table, ions: list[str]) -> tuple[dict[str, np.ndarray], dict[str, np.ndarray | None]]:
    values: dict[str, np.ndarray] = {}
    ranges: dict[str, np.ndarray | None] = {}
    for term_key in TERM_KEYS:
        if term_key == "dh":
            values[term_key] = np.zeros(len(ions), dtype=float)
            ranges[term_key] = None
        else:
            values[term_key] = _paper_values(frame, term_key, ions)
            ranges[term_key] = _paper_range(frame, term_key, ions)
    values["mu_sum"] = sum(values[term_key] for term_key in TERM_KEYS)
    values["adjusted_sum"] = values["mu_sum"].copy()
    values["total"] = totals.values("advanced", ions)
    return values, ranges


def _load_epcsaft_data(ions: list[str]) -> dict[str, np.ndarray]:
    values = {metric: np.empty(len(ions), dtype=float) for metric in METRICS}
    for idx, ion in enumerate(ions):
        mu = overlay.contribution_breakdown("advanced", ion, "water", basis="mu")
        lnf = overlay.contribution_breakdown("advanced", ion, "water", basis="lnfug")
        total = lnf["total"]
        for term_key in TERM_KEYS:
            values[term_key][idx] = mu[term_key]
        values["mu_sum"][idx] = sum(mu[term_key] for term_key in TERM_KEYS)
        values["adjusted_sum"][idx] = sum(lnf[term_key] for term_key in TERM_KEYS)
        values["total"][idx] = total
    return values


def _load_package_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def _package_payload_to_arrays(payload: dict, ions: list[str]) -> dict[str, np.ndarray]:
    results = payload["results"]
    values = {metric: np.empty(len(ions), dtype=float) for metric in METRICS}
    for idx, ion in enumerate(ions):
        row = results[ion]
        for term_key in TERM_KEYS:
            values[term_key][idx] = float(row[term_key])
        values["mu_sum"][idx] = float(row["mu_sum_kj_mol"])
        values["adjusted_sum"][idx] = float(row["lnfug_sum_kj_mol"])
        values["total"][idx] = float(row["total_kj_mol"])
    return values


def _run_clapeyron_extractor() -> dict:
    cmd = [
        "julia",
        f"--project={CLAPEYRON_ROOT}",
        str(CLAPEYRON_EXTRACTOR),
        str(CLAPEYRON_RAW_PATH),
        str(CLAPEYRON_ROOT),
    ]
    completed = subprocess.run(
        cmd,
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
        check=True,
    )
    if completed.stdout.strip():
        print(completed.stdout.strip(), flush=True)
    if completed.stderr.strip():
        print(completed.stderr.strip(), flush=True)
    return _load_package_json(CLAPEYRON_RAW_PATH)


def _validate_values(source_values: dict[str, dict[str, np.ndarray]], ions: list[str]) -> None:
    for source, metrics in source_values.items():
        for metric, values in metrics.items():
            if values.shape != (len(ions),):
                raise ValueError(f"{source} {metric} has shape {values.shape}, expected {(len(ions),)}")
            if not np.all(np.isfinite(values)):
                raise ValueError(f"{source} {metric} contains non-finite values.")


def _write_raw_long_csv(source_values: dict[str, dict[str, np.ndarray]], ions: list[str]) -> Path:
    with RAW_LONG_PATH.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["ion", "term", "source", "basis", "value_kj_mol"])
        writer.writeheader()
        for ion_idx, ion in enumerate(ions):
            for source in SOURCE_ORDER:
                for metric in METRICS:
                    writer.writerow(
                        {
                            "ion": ion,
                            "term": metric,
                            "source": source,
                            "basis": METRIC_BASIS[metric],
                            "value_kj_mol": f"{float(source_values[source][metric][ion_idx]):.12f}",
                        }
                    )
    return RAW_LONG_PATH


def _write_total_check_csv(source_values: dict[str, dict[str, np.ndarray]], ions: list[str]) -> Path:
    with TOTAL_CHECK_PATH.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "ion",
                "source",
                "total_kj_mol",
                "mu_sum_kj_mol",
                "adjusted_sum_kj_mol",
                "delta_total_minus_mu_sum",
                "delta_total_minus_adjusted_sum",
            ],
        )
        writer.writeheader()
        for ion_idx, ion in enumerate(ions):
            for source in SOURCE_ORDER:
                total = float(source_values[source]["total"][ion_idx])
                mu_sum = float(source_values[source]["mu_sum"][ion_idx])
                adjusted_sum = float(source_values[source]["adjusted_sum"][ion_idx])
                writer.writerow(
                    {
                        "ion": ion,
                        "source": source,
                        "total_kj_mol": f"{total:.12f}",
                        "mu_sum_kj_mol": f"{mu_sum:.12f}",
                        "adjusted_sum_kj_mol": f"{adjusted_sum:.12f}",
                        "delta_total_minus_mu_sum": f"{(total - mu_sum):.12f}",
                        "delta_total_minus_adjusted_sum": f"{(total - adjusted_sum):.12f}",
                    }
                )
    return TOTAL_CHECK_PATH


def _pct_delta(value: float, reference: float) -> float:
    denom = abs(reference)
    if denom <= 1.0e-12:
        return float("nan")
    return 100.0 * (value - reference) / denom


def _write_summary_csv(source_values: dict[str, dict[str, np.ndarray]], ions: list[str]) -> Path:
    with SUMMARY_PATH.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "ion",
                "term",
                "basis",
                "source",
                "value_kj_mol",
                "delta_vs_epcsaft_kj_mol",
                "pct_vs_epcsaft",
                "delta_vs_paper_kj_mol",
                "pct_vs_paper",
            ],
        )
        writer.writeheader()
        for ion_idx, ion in enumerate(ions):
            for metric in METRICS:
                paper_value = float(source_values["paper"][metric][ion_idx])
                epcsaft_value = float(source_values["epcsaft"][metric][ion_idx])
                for source in SOURCE_ORDER:
                    value = float(source_values[source][metric][ion_idx])
                    writer.writerow(
                        {
                            "ion": ion,
                            "term": metric,
                            "basis": METRIC_BASIS[metric],
                            "source": source,
                            "value_kj_mol": f"{value:.12f}",
                            "delta_vs_epcsaft_kj_mol": f"{(value - epcsaft_value):.12f}",
                            "pct_vs_epcsaft": f"{_pct_delta(value, epcsaft_value):.12f}",
                            "delta_vs_paper_kj_mol": f"{(value - paper_value):.12f}",
                            "pct_vs_paper": f"{_pct_delta(value, paper_value):.12f}",
                        }
                    )
    return SUMMARY_PATH


def _plot_term(
    term_key: str,
    term_label: str,
    color: str,
    ions: list[str],
    source_values: dict[str, dict[str, np.ndarray]],
    paper_yerr: np.ndarray | None,
) -> None:
    x = np.arange(len(ions), dtype=float)
    width = 0.18
    offsets = np.array([-1.5, -0.5, 0.5, 1.5]) * width

    fig, ax = plt.subplots(figsize=(11.8, 6.1))
    all_values: list[np.ndarray] = []
    for offset, source in zip(offsets, SOURCE_ORDER, strict=True):
        values = source_values[source][term_key]
        all_values.append(values)
        bars = ax.bar(
            x + offset,
            values,
            width=width,
            color=color,
            edgecolor="black",
            linewidth=0.45,
            hatch=SOURCE_HATCH[source],
            alpha=SOURCE_ALPHA[source],
            label=SOURCE_LABELS[source],
        )
        if source == "paper" and paper_yerr is not None:
            ax.errorbar(
                x + offset,
                values,
                yerr=paper_yerr,
                fmt="none",
                ecolor="0.15",
                elinewidth=1.1,
                capsize=5.0,
                capthick=1.1,
                zorder=6,
            )
        common.annotate_bar_values(ax, bars, fontsize=6)

    ax.axhline(0.0, color="black", linewidth=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(ions)
    ax.set_ylabel(r"Contribution to $\Delta G_{\mathrm{hyd},i}^{\infty}$ / kJ mol$^{-1}$ ($\mu$ basis)")
    ax.set_title(f"Bulow 2020 Figure 3 {term_label}: paper vs ePC-SAFT vs feos vs Clapeyron ($\\mu$ basis)")
    common.set_strict_bar_ylim(
        ax,
        np.asarray(all_values, dtype=float),
        step=5.0,
        top_pad_frac=0.20,
        bottom_pad_frac=0.12,
    )
    ax.grid(axis="y", alpha=0.22)
    ax.legend(ncol=2, frameon=True)

    output_path = SCRIPT_DIR / f"figure_3_{term_key}.png"
    common.save_figure(fig, output_path)
    plt.close(fig)
    print(f"Wrote {output_path}", flush=True)


def _plot_comprehensive(ions: list[str], source_values: dict[str, dict[str, np.ndarray]]) -> None:
    width = 0.09
    source_gap = 0.01
    term_gap = 0.11
    ion_gap = 0.80
    group_centers: list[float] = []
    tick_positions: list[float] = []
    tick_labels: list[str] = []
    separator_positions: list[float] = []
    x_cursor = 0.0

    x_positions: dict[tuple[str, str], list[float]] = {(source, term): [] for source in SOURCE_ORDER for term in TERM_KEYS}
    for ion_idx, _ion in enumerate(ions):
        ion_start = x_cursor
        for term_key in TERM_KEYS:
            term_start = x_cursor
            for source_idx, source in enumerate(SOURCE_ORDER):
                x_positions[(source, term_key)].append(x_cursor + source_idx * (width + source_gap))
            term_end = term_start + len(SOURCE_ORDER) * width + (len(SOURCE_ORDER) - 1) * source_gap
            tick_positions.append(0.5 * (term_start + term_end - width))
            tick_labels.append(term_key.upper())
            x_cursor = term_end + term_gap
        ion_end = x_cursor - term_gap
        group_centers.append(0.5 * (ion_start + ion_end))
        if ion_idx < len(ions) - 1:
            separator_positions.append(x_cursor + 0.5 * ion_gap)
        x_cursor += ion_gap

    fig, ax = plt.subplots(figsize=(20.0, 7.8))
    stacked_values = []
    for term_key, _term_label, color in TERMS:
        for source in SOURCE_ORDER:
            bars = ax.bar(
                x_positions[(source, term_key)],
                source_values[source][term_key],
                width=width,
                color=color,
                edgecolor="black",
                linewidth=0.4,
                hatch=SOURCE_HATCH[source],
                alpha=SOURCE_ALPHA[source],
            )
            stacked_values.append(source_values[source][term_key])
            if source == "paper":
                common.annotate_bar_values(ax, bars, fontsize=4)

    for xpos in separator_positions:
        ax.axvline(xpos, color="0.75", linewidth=0.8, linestyle="--", zorder=0)

    values = np.asarray(stacked_values, dtype=float)
    y_min = float(np.nanmin(values))
    y_max = float(np.nanmax(values))
    pad = max(15.0, 0.08 * max(abs(y_min), abs(y_max), 1.0))

    ax.axhline(0.0, color="black", linewidth=0.8)
    ax.set_ylim(y_min - pad, y_max + pad)
    ax.set_ylabel(r"Contribution to $\Delta G_{\mathrm{hyd},i}^{\infty}$ / kJ mol$^{-1}$ ($\mu$ basis)")
    ax.set_title("Bulow 2020 Figure 3 Comprehensive: paper vs ePC-SAFT vs feos vs Clapeyron")
    ax.set_xticks(tick_positions)
    ax.set_xticklabels(tick_labels)
    ax.grid(axis="y", alpha=0.22)

    for center, ion in zip(group_centers, ions, strict=True):
        ax.text(center, y_min - 0.70 * pad, ion, ha="center", va="bottom", fontsize=9, fontweight="bold")

    contribution_handles = [Patch(facecolor=color, edgecolor="black", label=label) for _key, label, color in TERMS]
    source_handles = [
        Patch(facecolor="white", edgecolor="black", hatch=SOURCE_HATCH[source], label=SOURCE_LABELS[source])
        for source in SOURCE_ORDER
    ]
    legend1 = ax.legend(handles=contribution_handles, loc="upper left", frameon=True, title="Contribution")
    ax.add_artist(legend1)
    ax.legend(handles=source_handles, loc="upper right", frameon=True, title="Source")

    output_path = SCRIPT_DIR / "figure_3_comprehensive.png"
    common.save_figure(fig, output_path)
    plt.close(fig)
    print(f"Wrote {output_path}", flush=True)


def _plot_total_check(ions: list[str], source_values: dict[str, dict[str, np.ndarray]]) -> None:
    x = np.arange(len(ions), dtype=float)
    order = [
        ("paper", "total"),
        ("paper", "mu_sum"),
        ("epcsaft", "total"),
        ("epcsaft", "adjusted_sum"),
        ("feos", "total"),
        ("feos", "adjusted_sum"),
        ("clapeyron", "total"),
        ("clapeyron", "adjusted_sum"),
    ]
    width = 0.095
    offsets = (np.arange(len(order), dtype=float) - (len(order) - 1) / 2.0) * width

    fig, ax = plt.subplots(figsize=(14.8, 6.7))
    stacked_values = []
    for offset, (source, metric) in zip(offsets, order, strict=True):
        values = source_values[source][metric]
        stacked_values.append(values)
        if metric == "total":
            facecolor = TOTAL_SOURCE_COLORS[source]
            edgecolor = "black"
            hatch = None
            alpha = 0.85
            linewidth = 0.45
        else:
            facecolor = "white"
            edgecolor = TOTAL_SOURCE_COLORS[source]
            hatch = SOURCE_HATCH[source]
            alpha = 1.0
            linewidth = 1.1
        bars = ax.bar(
            x + offset,
            values,
            width=width,
            color=facecolor,
            edgecolor=edgecolor,
            linewidth=linewidth,
            hatch=hatch,
            alpha=alpha,
        )
        common.annotate_bar_values(ax, bars, fontsize=5)

    values = np.asarray(stacked_values, dtype=float)
    y_min = float(np.nanmin(values))
    y_max = float(np.nanmax(values))
    pad = max(8.0, 0.09 * max(abs(y_min), abs(y_max), 1.0))

    ax.axhline(0.0, color="black", linewidth=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(ions)
    ax.set_ylabel(r"$\Delta G_{\mathrm{hyd},i}^{\infty}$ / kJ mol$^{-1}$ (total vs summed contributions)")
    ax.set_title("Bulow 2020 Figure 3 Total Check: paper-model totals and package fugacity-basis contribution sums")
    ax.set_ylim(y_min - pad, y_max + 2.0 * pad)
    ax.grid(axis="y", alpha=0.22)

    source_handles = [
        Patch(facecolor=TOTAL_SOURCE_COLORS[source], edgecolor="black", label=SOURCE_LABELS[source])
        for source in SOURCE_ORDER
    ]
    style_handles = [
        Patch(facecolor="0.75", edgecolor="black", label="Total"),
        Patch(facecolor="white", edgecolor="black", hatch="////", label="Contribution sum"),
    ]
    legend1 = ax.legend(handles=source_handles, loc="upper left", frameon=True, title="Source")
    ax.add_artist(legend1)
    ax.legend(handles=style_handles, loc="upper right", frameon=True, title="Bar style")

    output_path = SCRIPT_DIR / "figure_3_total.png"
    common.save_figure(fig, output_path)
    plt.close(fig)
    print(f"Wrote {output_path}", flush=True)


def main() -> None:
    common.configure_style()
    frame = common.load_indexed_csv(DATA_PATH)
    totals = common.load_indexed_csv(FIGURE2_TOTALS_PATH)
    ions = list(frame.columns)

    paper_values, paper_ranges = _load_paper_data(frame, totals, ions)
    epcsaft_values = _load_epcsaft_data(ions)

    feos_extractor.write_results(FEOS_RAW_PATH)
    feos_payload = _load_package_json(FEOS_RAW_PATH)
    feos_values = _package_payload_to_arrays(feos_payload, ions)

    clapeyron_payload = _run_clapeyron_extractor()
    clapeyron_values = _package_payload_to_arrays(clapeyron_payload, ions)

    source_values = {
        "paper": paper_values,
        "epcsaft": epcsaft_values,
        "feos": feos_values,
        "clapeyron": clapeyron_values,
    }
    _validate_values(source_values, ions)

    _write_raw_long_csv(source_values, ions)
    _write_total_check_csv(source_values, ions)
    _write_summary_csv(source_values, ions)

    for term_key, term_label, color in TERMS:
        _plot_term(term_key, term_label, color, ions, source_values, paper_ranges[term_key])
    _plot_comprehensive(ions, source_values)
    _plot_total_check(ions, source_values)


if __name__ == "__main__":
    main()

