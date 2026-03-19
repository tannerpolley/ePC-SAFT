from __future__ import annotations

from pathlib import Path
import sys

import numpy as np
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import _common as common

OUTPUT = Path(__file__).with_name('figure_4.png')
OUTPUT_COMBINED = Path(__file__).with_name('figure_4_combined.png')
OUTPUT_A = Path(__file__).with_name('figure_4a.png')
OUTPUT_B = Path(__file__).with_name('figure_4b.png')
DATA = Path(__file__).with_name('data') / 'water.csv'


def _load_literature() -> dict[str, float]:
    _, rows = common.read_csv_rows(DATA)
    out: dict[str, float] = {}
    for row in rows:
        ion = str(row.get('Ion', '')).strip()
        value = common.parse_float(row.get('Gsolv (kJ/mol)'))
        if ion and value is not None:
            out[f'{ion}+' if ion in {'H', 'Li', 'Na', 'K'} else f'{ion}-' if ion in {'Cl', 'Br', 'I'} else ion] = value
    return out


def _safe_model(dataset: str, ion: str) -> float:
    try:
        return -common.gsolv_ion(dataset, ion, 'water', {'water': 1.0})
    except Exception as exc:
        print(f'[figure_4] skipping {dataset} {ion}: {exc}')
        return float('nan')


def _bar_values(ions):
    literature = _load_literature()
    lit_vals = np.array([-literature.get(ion, np.nan) for ion in ions], dtype=float)
    figiel_vals = np.array([_safe_model('figiel_2025', ion) for ion in ions], dtype=float)
    bulow_vals = np.array([_safe_model('bulow_2020', ion) for ion in ions], dtype=float)
    return lit_vals, figiel_vals, bulow_vals


def _ion_labels(ions):
    mapping = {
        'Li+': r'$Li^+$',
        'Na+': r'$Na^+$',
        'K+': r'$K^+$',
        'Cl-': r'$Cl^-$',
        'Br-': r'$Br^-$',
        'I-': r'$I^-$',
    }
    return [mapping[ion] for ion in ions]


def _draw_bars(ax, ions, title: str | None = None, ylim=(0.0, 800.0)):
    x = np.arange(len(ions), dtype=float)
    width = 0.22
    lit_vals, figiel_vals, bulow_vals = _bar_values(ions)

    mask_lit = np.isfinite(lit_vals)
    mask_figiel = np.isfinite(figiel_vals)
    mask_bulow = np.isfinite(bulow_vals)

    ax.bar(x[mask_lit] - width, lit_vals[mask_lit], width=width, color=common.LIGHT_GRAY, edgecolor='black', linewidth=0.8, label='Literature')
    ax.bar(x[mask_figiel], figiel_vals[mask_figiel], width=width, color=common.BLUE_COLOR, edgecolor='black', linewidth=0.8, label='ePC-SAFT 2025')
    ax.bar(x[mask_bulow] + width, bulow_vals[mask_bulow], width=width, color=common.BROWN_COLOR, edgecolor='black', linewidth=0.8, label='ePC-SAFT 2020')

    ax.set_xticks(x)
    ax.set_xticklabels(_ion_labels(ions))
    if title:
        ax.set_title(title, fontsize=10)
    ax.set_ylim(*ylim)
    ax.set_ylabel(r'$-\Delta G_i^{solv,\infty,x}$ / kJ mol$^{-1}$')


def _plot_cations(ax) -> None:
    _draw_bars(ax, ['Li+', 'Na+', 'K+'], 'Cations in water')


def _plot_anions(ax) -> None:
    _draw_bars(ax, ['Cl-', 'Br-', 'I-'], 'Anions in water')


def _plot_combined(ax) -> None:
    _draw_bars(ax, ['Li+', 'Na+', 'K+', 'Cl-', 'Br-', 'I-'], 'Gibbs energy of solvation at infinite dilution in water of different ions at 298.15 K and 1 bar.')


def main() -> None:
    common.configure_style()

    fig_combined, ax = plt.subplots(1, 1, figsize=(8.0, 4.2))
    _plot_combined(ax)
    handles, labels = ax.get_legend_handles_labels()
    fig_combined.legend(handles, labels, loc='upper center', ncol=3, bbox_to_anchor=(0.5, 0.98), fontsize=9, frameon=False)
    fig_combined.subplots_adjust(left=0.09, right=0.98, bottom=0.14, top=0.82)
    common.save_figure(fig_combined, OUTPUT)
    common.save_figure(fig_combined, OUTPUT_COMBINED)
    plt.close(fig_combined)

    common.save_panel_figure(_plot_cations, OUTPUT_A, figsize=(4.1, 3.6), legend_handles=handles, legend_labels=labels, legend_kwargs={'loc': 'upper right', 'fontsize': 7})
    common.save_panel_figure(_plot_anions, OUTPUT_B, figsize=(4.1, 3.6), legend_handles=handles, legend_labels=labels, legend_kwargs={'loc': 'upper right', 'fontsize': 7})


if __name__ == '__main__':
    main()

