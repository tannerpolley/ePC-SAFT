from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import run_package_comparisons as rpc


REDUCED_SUFFIX = '_reduced'
REDUCED_SOURCES = ('paper', 'pcsaft', 'feos')
REDUCED_LABELS = {
    'paper': 'Paper',
    'pcsaft': 'Current repo PC-SAFT',
    'feos': 'feos (patched)',
}
REDUCED_HATCH = {
    'paper': None,
    'pcsaft': '////',
    'feos': '\\\\\\\\',
}
REDUCED_ALPHA = {
    'paper': 1.0,
    'pcsaft': 0.82,
    'feos': 0.72,
}
REDUCED_TOTAL_COLORS = {
    'paper': '#7F7F7F',
    'pcsaft': '#2CA02C',
    'feos': '#1F77B4',
}


def _plot_term_reduced(
    term_key: str,
    term_label: str,
    color: str,
    ions: list[str],
    source_values: dict[str, dict[str, np.ndarray]],
    paper_yerr: np.ndarray | None,
) -> Path:
    x = np.arange(len(ions), dtype=float)
    width = 0.22
    offsets = (np.arange(len(REDUCED_SOURCES), dtype=float) - (len(REDUCED_SOURCES) - 1) / 2.0) * width

    fig, ax = plt.subplots(figsize=(11.8, 6.1))
    all_values: list[np.ndarray] = []
    for offset, source in zip(offsets, REDUCED_SOURCES, strict=True):
        values = source_values[source][term_key]
        all_values.append(values)
        bars = ax.bar(
            x + offset,
            values,
            width=width,
            color=color,
            edgecolor='black',
            linewidth=0.45,
            hatch=REDUCED_HATCH[source],
            alpha=REDUCED_ALPHA[source],
            label=REDUCED_LABELS[source],
        )
        if source == 'paper' and paper_yerr is not None:
            ax.errorbar(
                x + offset,
                values,
                yerr=paper_yerr,
                fmt='none',
                ecolor='0.15',
                elinewidth=1.1,
                capsize=5.0,
                capthick=1.1,
                zorder=6,
            )
        rpc.common.annotate_bar_values(ax, bars, fontsize=6)

    ax.axhline(0.0, color='black', linewidth=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(ions)
    ax.set_ylabel(r'Contribution to $\Delta G_{\mathrm{hyd},i}^{\infty}$ / kJ mol$^{-1}$ ($\mu$ basis)')
    ax.set_title(f'Bulow 2020 Figure 3 {term_label}: paper vs PC-SAFT vs feos ($\\mu$ basis)')
    rpc.common.set_strict_bar_ylim(
        ax,
        np.asarray(all_values, dtype=float),
        step=5.0,
        top_pad_frac=0.20,
        bottom_pad_frac=0.12,
    )
    ax.grid(axis='y', alpha=0.22)
    ax.legend(ncol=2, frameon=True)

    output_path = SCRIPT_DIR / f'figure_3_{term_key}{REDUCED_SUFFIX}.png'
    rpc.common.save_figure(fig, output_path)
    plt.close(fig)
    print(f'Wrote {output_path}', flush=True)
    return output_path


def _plot_comprehensive_reduced(ions: list[str], source_values: dict[str, dict[str, np.ndarray]]) -> Path:
    width = 0.11
    source_gap = 0.015
    term_gap = 0.13
    ion_gap = 0.85
    group_centers: list[float] = []
    tick_positions: list[float] = []
    tick_labels: list[str] = []
    separator_positions: list[float] = []
    x_cursor = 0.0

    x_positions: dict[tuple[str, str], list[float]] = {(source, term): [] for source in REDUCED_SOURCES for term in rpc.TERM_KEYS}
    for ion_idx, _ion in enumerate(ions):
        ion_start = x_cursor
        for term_key in rpc.TERM_KEYS:
            term_start = x_cursor
            for source_idx, source in enumerate(REDUCED_SOURCES):
                x_positions[(source, term_key)].append(x_cursor + source_idx * (width + source_gap))
            term_end = term_start + len(REDUCED_SOURCES) * width + (len(REDUCED_SOURCES) - 1) * source_gap
            tick_positions.append(0.5 * (term_start + term_end - width))
            tick_labels.append(term_key.upper())
            x_cursor = term_end + term_gap
        ion_end = x_cursor - term_gap
        group_centers.append(0.5 * (ion_start + ion_end))
        if ion_idx < len(ions) - 1:
            separator_positions.append(x_cursor + 0.5 * ion_gap)
        x_cursor += ion_gap

    fig, ax = plt.subplots(figsize=(18.0, 7.2))
    stacked_values = []
    for term_key, _term_label, color in rpc.TERMS:
        for source in REDUCED_SOURCES:
            bars = ax.bar(
                x_positions[(source, term_key)],
                source_values[source][term_key],
                width=width,
                color=color,
                edgecolor='black',
                linewidth=0.4,
                hatch=REDUCED_HATCH[source],
                alpha=REDUCED_ALPHA[source],
            )
            stacked_values.append(source_values[source][term_key])
            if source == 'paper':
                rpc.common.annotate_bar_values(ax, bars, fontsize=4)

    for xpos in separator_positions:
        ax.axvline(xpos, color='0.75', linewidth=0.8, linestyle='--', zorder=0)

    values = np.asarray(stacked_values, dtype=float)
    y_min = float(np.nanmin(values))
    y_max = float(np.nanmax(values))
    pad = max(15.0, 0.08 * max(abs(y_min), abs(y_max), 1.0))

    ax.axhline(0.0, color='black', linewidth=0.8)
    ax.set_ylim(y_min - pad, y_max + pad)
    ax.set_ylabel(r'Contribution to $\Delta G_{\mathrm{hyd},i}^{\infty}$ / kJ mol$^{-1}$ ($\mu$ basis)')
    ax.set_title('Bulow 2020 Figure 3 Comprehensive: paper vs PC-SAFT vs feos')
    ax.set_xticks(tick_positions)
    ax.set_xticklabels(tick_labels)
    ax.grid(axis='y', alpha=0.22)

    for center, ion in zip(group_centers, ions, strict=True):
        ax.text(center, y_min - 0.70 * pad, ion, ha='center', va='bottom', fontsize=9, fontweight='bold')

    contribution_handles = [Patch(facecolor=color, edgecolor='black', label=label) for _key, label, color in rpc.TERMS]
    source_handles = [Patch(facecolor='white', edgecolor='black', hatch=REDUCED_HATCH[source], label=REDUCED_LABELS[source]) for source in REDUCED_SOURCES]
    legend1 = ax.legend(handles=contribution_handles, loc='upper left', frameon=True, title='Contribution')
    ax.add_artist(legend1)
    ax.legend(handles=source_handles, loc='upper right', frameon=True, title='Source')

    output_path = SCRIPT_DIR / f'figure_3_comprehensive{REDUCED_SUFFIX}.png'
    rpc.common.save_figure(fig, output_path)
    plt.close(fig)
    print(f'Wrote {output_path}', flush=True)
    return output_path


def _plot_total_check_reduced(ions: list[str], source_values: dict[str, dict[str, np.ndarray]]) -> Path:
    x = np.arange(len(ions), dtype=float)
    order = [
        ('paper', 'total'),
        ('paper', 'mu_sum'),
        ('pcsaft', 'total'),
        ('pcsaft', 'adjusted_sum'),
        ('feos', 'total'),
        ('feos', 'adjusted_sum'),
    ]
    width = 0.12
    offsets = (np.arange(len(order), dtype=float) - (len(order) - 1) / 2.0) * width

    fig, ax = plt.subplots(figsize=(13.6, 6.4))
    stacked_values = []
    for offset, (source, metric) in zip(offsets, order, strict=True):
        values = source_values[source][metric]
        stacked_values.append(values)
        if metric == 'total':
            facecolor = REDUCED_TOTAL_COLORS[source]
            edgecolor = 'black'
            hatch = None
            alpha = 0.85
            linewidth = 0.45
        else:
            facecolor = 'white'
            edgecolor = REDUCED_TOTAL_COLORS[source]
            hatch = REDUCED_HATCH[source]
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
        rpc.common.annotate_bar_values(ax, bars, fontsize=5)

    values = np.asarray(stacked_values, dtype=float)
    y_min = float(np.nanmin(values))
    y_max = float(np.nanmax(values))
    pad = max(8.0, 0.09 * max(abs(y_min), abs(y_max), 1.0))

    ax.axhline(0.0, color='black', linewidth=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(ions)
    ax.set_ylabel(r'$\Delta G_{\mathrm{hyd},i}^{\infty}$ / kJ mol$^{-1}$ (total vs summed contributions)')
    ax.set_title('Bulow 2020 Figure 3 Total Check: paper-model totals and package fugacity-basis contribution sums')
    ax.set_ylim(y_min - pad, y_max + 2.0 * pad)
    ax.grid(axis='y', alpha=0.22)

    source_handles = [Patch(facecolor=REDUCED_TOTAL_COLORS[source], edgecolor='black', label=REDUCED_LABELS[source]) for source in REDUCED_SOURCES]
    style_handles = [
        Patch(facecolor='0.75', edgecolor='black', label='Total'),
        Patch(facecolor='white', edgecolor='black', hatch='////', label='Contribution sum'),
    ]
    legend1 = ax.legend(handles=source_handles, loc='upper left', frameon=True, title='Source')
    ax.add_artist(legend1)
    ax.legend(handles=style_handles, loc='upper right', frameon=True, title='Bar style')

    output_path = SCRIPT_DIR / f'figure_3_total{REDUCED_SUFFIX}.png'
    rpc.common.save_figure(fig, output_path)
    plt.close(fig)
    print(f'Wrote {output_path}', flush=True)
    return output_path


def main() -> None:
    rpc.common.configure_style()
    rpc.SOURCE_ORDER = REDUCED_SOURCES
    rpc.RAW_LONG_PATH = SCRIPT_DIR / f'raw_comparison_long{REDUCED_SUFFIX}.csv'
    rpc.TOTAL_CHECK_PATH = SCRIPT_DIR / f'total_check{REDUCED_SUFFIX}.csv'
    rpc.SUMMARY_PATH = SCRIPT_DIR / f'comparison_summary{REDUCED_SUFFIX}.csv'
    reduced_feos_path = SCRIPT_DIR / f'feos_raw{REDUCED_SUFFIX}.json'

    frame = rpc.common.load_indexed_csv(rpc.DATA_PATH)
    totals = rpc.common.load_indexed_csv(rpc.FIGURE2_TOTALS_PATH)
    ions = list(frame.columns)

    paper_values, paper_ranges = rpc._load_paper_data(frame, totals, ions)
    pcsaft_values = rpc._load_pcsaft_data(ions)

    rpc.feos_extractor.write_results(reduced_feos_path)
    feos_payload = rpc._load_package_json(reduced_feos_path)
    feos_values = rpc._package_payload_to_arrays(feos_payload, ions)

    source_values = {
        'paper': paper_values,
        'pcsaft': pcsaft_values,
        'feos': feos_values,
    }
    rpc._validate_values(source_values, ions)

    rpc._write_raw_long_csv(source_values, ions)
    rpc._write_total_check_csv(source_values, ions)
    rpc._write_summary_csv(source_values, ions)

    for term_key, term_label, color in rpc.TERMS:
        _plot_term_reduced(term_key, term_label, color, ions, source_values, paper_ranges[term_key])
    _plot_comprehensive_reduced(ions, source_values)
    _plot_total_check_reduced(ions, source_values)


if __name__ == '__main__':
    main()
