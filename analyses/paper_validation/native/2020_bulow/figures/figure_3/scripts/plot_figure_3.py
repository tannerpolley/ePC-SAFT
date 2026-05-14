from __future__ import annotations

import sys


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

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

SCRIPT_DIR = Path(__file__).resolve().parent
ANALYSIS_ROOT = SCRIPT_DIR.parents[2]

if str(ANALYSIS_ROOT) not in sys.path:
    sys.path.insert(0, str(ANALYSIS_ROOT))
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

import _plot_common as common
import _model_overlay as overlay

DATA_PATH = common.analysis_data_path(__file__, "water_contributions.csv", kind="input")
FIGURE2_TOTALS_PATH = common.analysis_data_path(__file__, "water_comparisons.csv", kind="input")
CONTRIBUTIONS = [
    ("hc", "Hard chain", "#9f9f9f"),
    ("disp", "Dispersion", "#5f5f5f"),
    ("assoc", "Association", "#111111"),
    ("dh", "Debye-Huckel", "#8c564b"),
    ("born", "Born", "#d8891c"),
]
COMPREHENSIVE_KEYS = ("hc", "disp", "assoc", "born")
TOTAL_GREEN = "#2ca02c"
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


def _model_values(ions: list[str], term_key: str, basis: str = "mu", d_born_mode: int | None = None) -> np.ndarray:
    out = np.empty(len(ions), dtype=float)
    for idx, ion in enumerate(ions):
        out[idx] = overlay.contribution_breakdown("advanced", ion, "water", basis=basis, d_born_mode=d_born_mode)[
            term_key
        ]
    return out


def _paper_values(frame: common.Table, term_key: str, ions: list[str]) -> np.ndarray:
    for row_key in PAPER_ROW_KEYS.get(term_key, (term_key,)):
        if row_key in frame.index:
            return frame.values(row_key, ions)
    raise KeyError(
        f"Missing paper row for Figure 3 term '{term_key}'. Tried: {PAPER_ROW_KEYS.get(term_key, (term_key,))}"
    )


def _paper_range(frame: common.Table, term_key: str, ions: list[str]) -> tuple[np.ndarray, np.ndarray] | None:
    row_keys = PAPER_RANGE_ROW_KEYS.get(term_key)
    if row_keys is None:
        return None
    if not all(row_key in frame.index for row_key in row_keys):
        return None
    low = frame.values(row_keys[0], ions)
    avg = frame.values(row_keys[1], ions)
    high = frame.values(row_keys[2], ions)
    yerr = np.vstack([np.abs(avg - low), np.abs(high - avg)])
    return avg, yerr


def _model_zquotient_sum(ions: list[str]) -> np.ndarray:
    out = np.empty(len(ions), dtype=float)
    for idx, ion in enumerate(ions):
        mu = overlay.contribution_breakdown("advanced", ion, "water", basis="mu")
        lnf = overlay.contribution_breakdown("advanced", ion, "water", basis="lnfug")
        out[idx] = sum(lnf[key] - mu[key] for key in ("hc", "disp", "assoc", "dh", "born"))
    return out


def _plot_total(
    ions: list[str],
    figure2_paper: np.ndarray,
    figure2_model: np.ndarray,
    figure3_sum_paper: np.ndarray,
    figure3_sum_model: np.ndarray,
) -> None:
    x = np.arange(len(ions), dtype=float)
    width = 0.17
    offsets = np.array([-1.5, -0.5, 0.5, 1.5]) * width
    fig, ax = plt.subplots(figsize=(11.4, 6.0))

    bar_specs = [
        (offsets[0], figure2_paper, TOTAL_GREEN, "black", None, 0.9, "Figure 2 advanced total (paper)"),
        (offsets[1], figure2_model, TOTAL_GREEN, "black", "////", 0.75, "Figure 2 advanced total (epcsaft)"),
        (offsets[2], figure3_sum_paper, "white", TOTAL_GREEN, None, 1.0, "Figure 3 contribution sum (paper)"),
        (offsets[3], figure3_sum_model, "white", TOTAL_GREEN, "////", 1.0, "Figure 3 contribution sum (epcsaft)"),
    ]

    legend_added: set[str] = set()
    for offset, values, facecolor, edgecolor, hatch, alpha, label in bar_specs:
        show_label = label if label not in legend_added else None
        bars = ax.bar(
            x + offset,
            values,
            width=width,
            color=facecolor,
            edgecolor=edgecolor,
            linewidth=1.0 if facecolor == "white" else 0.45,
            hatch=hatch,
            alpha=alpha,
            label=show_label,
        )
        common.annotate_bar_values(ax, bars, fontsize=6)
        if show_label is not None:
            legend_added.add(label)

    stacked = np.vstack([figure2_paper, figure2_model, figure3_sum_paper, figure3_sum_model])
    y_min = float(np.nanmin(stacked))
    y_max = float(np.nanmax(stacked))
    pad = max(3.0, 0.10 * max(abs(y_min), abs(y_max), 1.0))

    ax.axhline(0.0, color="black", linewidth=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(ions)
    ax.set_ylabel(r"$\Delta G_{\mathrm{hyd},i}^{\infty}$ / kJ mol$^{-1}$ (EOS total / EOS-sum check)")
    ax.set_title(
        "Bulow 2020 Part I Figure 3 Total Check: Figure 2 advanced totals vs summed Figure 3 contribution bars"
    )
    ax.set_ylim(y_min - pad, y_max + 2.7 * pad)
    ax.grid(axis="y", alpha=0.22)
    common.annotate_percent_deltas(ax, x + offsets[1], figure2_paper, figure2_model, xs_ref=x + offsets[0])
    common.annotate_percent_deltas(ax, x + offsets[3], figure3_sum_paper, figure3_sum_model, xs_ref=x + offsets[2])
    ax.text(
        0.99,
        0.12,
        "Filled green bars: Figure 2 advanced totals\nWhite/green bars: Figure 3 sum\nPaper sum uses published bars; epcsaft sum uses EOS-consistent lnfug contributions",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=8,
        bbox={"facecolor": "white", "edgecolor": "0.5", "alpha": 0.9, "boxstyle": "round,pad=0.2"},
    )
    common.add_percent_note(ax)
    ax.legend(ncol=2, frameon=True)

    output_path = SCRIPT_DIR / "figure_3_total.png"
    common.save_figure(fig, output_path)
    plt.close(fig)
    print(f"Wrote {output_path}", flush=True)


def _plot_z_adjusted_total(
    ions: list[str],
    figure2_paper: np.ndarray,
    figure2_model: np.ndarray,
    figure3_sum_paper: np.ndarray,
    figure3_mu_sum_model: np.ndarray,
    zquotient_sum_model: np.ndarray,
) -> None:
    adjusted_paper = figure3_sum_paper + zquotient_sum_model
    adjusted_model = figure3_mu_sum_model + zquotient_sum_model

    x = np.arange(len(ions), dtype=float)
    width = 0.13
    offsets = np.array([-2.5, -1.5, -0.5, 0.5, 1.5, 2.5]) * width
    fig, ax = plt.subplots(figsize=(13.2, 6.3))

    bars1 = ax.bar(
        x + offsets[0],
        figure2_paper,
        width=width,
        color=TOTAL_GREEN,
        edgecolor="black",
        linewidth=0.45,
        label="Figure 2 advanced total (paper)",
    )
    bars2 = ax.bar(
        x + offsets[1],
        figure3_sum_paper,
        width=width,
        color="#bdbdbd",
        edgecolor="black",
        linewidth=0.45,
        label="Figure 3 paper sum",
    )
    bars3 = ax.bar(
        x + offsets[2],
        adjusted_paper,
        width=width,
        color="#4c78a8",
        edgecolor="black",
        linewidth=0.45,
        hatch="////",
        alpha=0.8,
        label=r"Figure 3 paper sum + model $\sum_{\alpha}\left[-\frac{Z^\alpha}{Z-1}\ln Z\right]$",
    )
    bars4 = ax.bar(
        x + offsets[3],
        figure2_model,
        width=width,
        color=TOTAL_GREEN,
        edgecolor="black",
        linewidth=0.45,
        hatch="////",
        alpha=0.75,
        label="Figure 2 advanced total (epcsaft)",
    )
    bars5 = ax.bar(
        x + offsets[4],
        figure3_mu_sum_model,
        width=width,
        color="white",
        edgecolor="black",
        linewidth=1.0,
        label=r"Figure 3 epcsaft $\mu$-sum",
    )
    bars6 = ax.bar(
        x + offsets[5],
        adjusted_model,
        width=width,
        color="#4c78a8",
        edgecolor="black",
        linewidth=0.45,
        hatch="xx",
        alpha=0.85,
        label=r"Figure 3 epcsaft $\mu$-sum + model $\sum_{\alpha}\left[-\frac{Z^\alpha}{Z-1}\ln Z\right]$",
    )
    common.annotate_bar_values(ax, bars1, fontsize=6)
    common.annotate_bar_values(ax, bars2, fontsize=6)
    common.annotate_bar_values(ax, bars3, fontsize=6)
    common.annotate_bar_values(ax, bars4, fontsize=6)
    common.annotate_bar_values(ax, bars5, fontsize=6)
    common.annotate_bar_values(ax, bars6, fontsize=6)

    stacked = np.vstack(
        [figure2_paper, figure3_sum_paper, adjusted_paper, figure2_model, figure3_mu_sum_model, adjusted_model]
    )
    y_min = float(np.nanmin(stacked))
    y_max = float(np.nanmax(stacked))
    pad = max(3.0, 0.10 * max(abs(y_min), abs(y_max), 1.0))

    before_gap = figure2_paper - figure3_sum_paper
    after_gap = figure2_paper - adjusted_paper
    mae_before = float(np.mean(np.abs(before_gap)))
    mae_after = float(np.mean(np.abs(after_gap)))
    model_identity_gap = float(np.max(np.abs(adjusted_model - figure2_model)))

    ax.axhline(0.0, color="black", linewidth=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(ions)
    ax.set_ylabel(r"$\Delta G_{\mathrm{hyd},i}^{\infty}$ / kJ mol$^{-1}$")
    ax.set_title(r"Figure 3 bookkeeping check: paper and epcsaft sums plus summed $Z$-quotient correction")
    ax.set_ylim(y_min - pad, y_max + 2.5 * pad)
    ax.grid(axis="y", alpha=0.22)
    common.annotate_percent_deltas(ax, x + offsets[2], figure2_paper, adjusted_paper, xs_ref=x + offsets[0], fontsize=7)
    common.annotate_percent_deltas(ax, x + offsets[5], figure2_model, adjusted_model, xs_ref=x + offsets[3], fontsize=7)
    ax.text(
        0.99,
        0.12,
        "Paper side:\n"
        f"raw Figure 3 sum gap = {mae_before:.2f} kJ/mol\n"
        f"+ summed model Z quotient gap = {mae_after:.2f} kJ/mol\n"
        "Model side:\n"
        f"max |(mu-sum + Z quotient) - Figure 2 total| = {model_identity_gap:.2e} kJ/mol",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=8,
        bbox={"facecolor": "white", "edgecolor": "0.5", "alpha": 0.9, "boxstyle": "round,pad=0.2"},
    )
    common.add_percent_note(ax)
    ax.legend(frameon=True)

    output_path = SCRIPT_DIR / "figure_3_zquotient_total.png"
    common.save_figure(fig, output_path)
    plt.close(fig)
    print(f"Wrote {output_path}", flush=True)


def _plot_one(
    term_key: str,
    term_label: str,
    color: str,
    ions: list[str],
    paper: np.ndarray,
    model: np.ndarray,
    paper_yerr: np.ndarray | None = None,
) -> None:
    x = np.arange(len(ions), dtype=float)
    width = 0.34

    fig, ax = plt.subplots(figsize=(10.8, 5.8))
    paper_offset = -width / 2.0
    model_offset = width / 2.0
    paper_bars = ax.bar(
        x + paper_offset,
        paper,
        width=width,
        color=color,
        edgecolor="black",
        linewidth=0.45,
        label=f"{term_label} (paper)",
    )
    if paper_yerr is not None:
        ax.errorbar(
            x + paper_offset,
            paper,
            yerr=paper_yerr,
            fmt="none",
            ecolor="0.15",
            elinewidth=1.1,
            capsize=5.0,
            capthick=1.1,
            zorder=6,
        )
    model_bars = ax.bar(
        x + model_offset,
        model,
        width=width,
        color=color,
        edgecolor="black",
        linewidth=0.45,
        hatch="////",
        alpha=0.75,
        label=f"{term_label} (epcsaft)",
    )
    common.annotate_bar_values(ax, paper_bars, fontsize=6)
    common.annotate_bar_values(ax, model_bars, fontsize=6)

    values = np.vstack([paper, model])

    ax.axhline(0.0, color="black", linewidth=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(ions)
    ax.set_ylabel(r"Contribution to $\Delta G_{\mathrm{hyd},i}^{\infty}$ / kJ mol$^{-1}$ ($\mu$ basis)")
    ax.set_title(f"Bulow 2020 Part I Figure 3: {term_label} ($\\mu$ basis)")
    common.set_strict_bar_ylim(ax, values, step=5.0, top_pad_frac=0.18, bottom_pad_frac=0.12)
    ax.grid(axis="y", alpha=0.22)
    common.annotate_percent_deltas(ax, x + model_offset, paper, model, xs_ref=x + paper_offset)
    common.add_percent_note(ax)
    ax.legend(frameon=True)

    output_path = SCRIPT_DIR / f"figure_3_{term_key}.png"
    common.save_figure(fig, output_path)
    plt.close(fig)
    print(f"Wrote {output_path}", flush=True)


def _plot_comprehensive(
    ions: list[str],
    paper_contributions: dict[str, np.ndarray],
    paper_ranges: dict[str, np.ndarray | None],
    model_contributions: dict[str, np.ndarray],
) -> None:
    term_meta = {key: (label, color) for key, label, color in CONTRIBUTIONS if key in COMPREHENSIVE_KEYS}
    width = 0.16
    term_gap = 0.10
    ion_gap = 0.95
    group_centers: list[float] = []
    paper_x: list[float] = []
    model_x: list[float] = []
    tick_labels: list[str] = []
    separator_positions: list[float] = []
    x_cursor = 0.0

    for ion_idx, ion in enumerate(ions):
        start_x = x_cursor
        for term_key in COMPREHENSIVE_KEYS:
            paper_x.append(x_cursor)
            model_x.append(x_cursor + width)
            tick_labels.append(term_key.upper())
            x_cursor += 2.0 * width + term_gap
        end_x = x_cursor - term_gap
        group_centers.append(0.5 * (start_x + end_x))
        if ion_idx < len(ions) - 1:
            separator_positions.append(x_cursor + 0.5 * ion_gap)
        x_cursor += ion_gap

    fig, ax = plt.subplots(figsize=(18.0, 7.2))
    values_for_limits = []

    bar_index = 0
    for ion_idx, ion in enumerate(ions):
        for term_key in COMPREHENSIVE_KEYS:
            _, color = term_meta[term_key]
            paper_value = float(paper_contributions[term_key][ion_idx])
            model_value = float(model_contributions[term_key][ion_idx])
            values_for_limits.extend([paper_value, model_value])
            paper_bars = ax.bar(
                paper_x[bar_index],
                paper_value,
                width=width,
                color=color,
                edgecolor="black",
                linewidth=0.45,
            )
            if paper_ranges.get(term_key) is not None:
                ax.errorbar(
                    paper_x[bar_index],
                    paper_value,
                    yerr=np.asarray(paper_ranges[term_key])[:, ion_idx : ion_idx + 1],
                    fmt="none",
                    ecolor="0.15",
                    elinewidth=0.9,
                    capsize=3.5,
                    capthick=0.9,
                    zorder=6,
                )
            model_bars = ax.bar(
                model_x[bar_index],
                model_value,
                width=width,
                color=color,
                edgecolor="black",
                linewidth=0.45,
                hatch="////",
                alpha=0.75,
            )
            common.annotate_bar_values(ax, paper_bars, fontsize=5)
            common.annotate_bar_values(ax, model_bars, fontsize=5)
            bar_index += 1

    for xpos in separator_positions:
        ax.axvline(xpos, color="0.75", linewidth=0.8, linestyle="--", zorder=0)

    ax.axhline(0.0, color="black", linewidth=0.8)
    ax.set_ylim(-1000.0, 50.0)
    ax.set_ylabel(r"Contribution to $\Delta G_{\mathrm{hyd},i}^{\infty}$ / kJ mol$^{-1}$ ($\mu$ basis)")
    ax.set_title("Bulow 2020 Part I Figure 3 Comprehensive: paper vs epcsaft contribution bars ($\\mu$ basis)")
    ax.set_xticks([(px + mx) / 2.0 for px, mx in zip(paper_x, model_x)])
    ax.set_xticklabels(tick_labels, rotation=0)
    ax.grid(axis="y", alpha=0.22)

    for center, ion in zip(group_centers, ions):
        ax.text(center, -985.0, ion, ha="center", va="bottom", fontsize=9, fontweight="bold")

    contribution_handles = [
        Patch(facecolor=term_meta[key][1], edgecolor="black", label=term_meta[key][0]) for key in COMPREHENSIVE_KEYS
    ]
    basis_handles = [
        Patch(facecolor="white", edgecolor="black", label="Paper"),
        Patch(facecolor="white", edgecolor="black", hatch="////", label="epcsaft"),
    ]
    legend1 = ax.legend(handles=contribution_handles, loc="upper left", frameon=True, title="Contribution")
    ax.add_artist(legend1)
    ax.legend(handles=basis_handles, loc="upper right", frameon=True, title="Bar style")

    output_path = SCRIPT_DIR / "figure_3_comprehensive.png"
    common.save_figure(fig, output_path)
    plt.close(fig)
    print(f"Wrote {output_path}", flush=True)


def main() -> None:
    common.configure_style()
    frame = common.load_indexed_csv(DATA_PATH)
    totals = common.load_indexed_csv(FIGURE2_TOTALS_PATH)
    ions = list(frame.columns)
    paper_contributions: dict[str, np.ndarray] = {}
    paper_ranges: dict[str, np.ndarray | None] = {}
    model_contributions: dict[str, np.ndarray] = {}
    for term_key, term_label, color in CONTRIBUTIONS:
        print(f"Computing Figure 3 {term_key}", flush=True)
        if term_key == "dh":
            paper = np.zeros(len(ions), dtype=float)
            paper_yerr = None
        else:
            paper = _paper_values(frame, term_key, ions)
            range_info = _paper_range(frame, term_key, ions)
            paper_yerr = None if range_info is None else range_info[1]
        model = _model_values(ions, term_key, basis="mu")
        paper_contributions[term_key] = paper
        paper_ranges[term_key] = paper_yerr
        model_contributions[term_key] = model
        _plot_one(term_key, term_label, color, ions, paper, model, paper_yerr=paper_yerr)

    _plot_comprehensive(ions, paper_contributions, paper_ranges, model_contributions)

    figure2_paper = totals.values("advanced", ions)
    figure2_model = np.asarray([overlay.gsolv_ion("advanced", ion, "water") for ion in ions], dtype=float)
    figure3_sum_paper = sum(paper_contributions[key] for key, _, _ in CONTRIBUTIONS)
    figure3_mu_sum_model = sum(_model_values(ions, key, basis="mu") for key, _, _ in CONTRIBUTIONS)
    figure3_sum_model = sum(_model_values(ions, key, basis="lnfug") for key, _, _ in CONTRIBUTIONS)
    zquotient_sum_model = _model_zquotient_sum(ions)
    _plot_total(ions, figure2_paper, figure2_model, figure3_sum_paper, figure3_sum_model)
    _plot_z_adjusted_total(
        ions, figure2_paper, figure2_model, figure3_sum_paper, figure3_mu_sum_model, zquotient_sum_model
    )


if __name__ == "__main__":
    main()

