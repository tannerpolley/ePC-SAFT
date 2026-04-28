from __future__ import annotations

from functools import lru_cache
from pathlib import Path
import sys

import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import _common as common
import figure_data

OUTPUT = Path(__file__).with_name('figure_8.png')
OUTPUTS = [
    Path(__file__).with_name('figure_8a.png'),
    Path(__file__).with_name('figure_8b.png'),
    Path(__file__).with_name('figure_8c.png'),
]
PANELS = figure_data.FIG8_PANELS


@lru_cache(maxsize=1)
def _payload():
    return tuple(figure_data.read_payload('figure_8'))


def _plot_panel(ax, label, salt, m_max, y_max, include_legend: bool = False):
    rows = list(_payload())
    x_data_meoh, y_data_meoh = figure_data.xy(figure_data.select_rows(rows, panel_id=label, series_id=f'data_{salt}_methanol'))
    x_data_etoh, y_data_etoh = figure_data.xy(figure_data.select_rows(rows, panel_id=label, series_id=f'data_{salt}_ethanol'))
    m_grid_meoh, y_meoh = figure_data.xy(figure_data.select_rows(rows, panel_id=label, series_id=f'model_{salt}_methanol'))
    m_grid_etoh, y_etoh = figure_data.xy(figure_data.select_rows(rows, panel_id=label, series_id=f'model_{salt}_ethanol'))
    ax.scatter(x_data_meoh, y_data_meoh, s=24, facecolor='none', edgecolor=common.GRAY_COLOR, linewidth=0.9, label='Methanol data')
    ax.scatter(x_data_etoh, y_data_etoh, s=24, marker='s', facecolor=common.GREEN_COLOR, edgecolor=common.GREEN_COLOR, linewidth=0.8, label='Ethanol data')
    ax.plot(m_grid_meoh, y_meoh, color=common.GRAY_COLOR, linewidth=1.5, label='Methanol model')
    ax.plot(m_grid_etoh, y_etoh, color='black', linewidth=1.5, label='Ethanol model')
    ax.set_xlim(0.0, m_max)
    ax.set_ylim(0.0, y_max)
    ax.set_title(salt, fontsize=10)
    ax.set_xlabel(r'$\bar{m}_{salt}$ / mol kg$^{-1}$')
    ax.set_ylabel(r'$\gamma_{\pm}^{m,*}$ / -')
    if include_legend:
        ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.20), ncol=2, fontsize=7.5, frameon=False, columnspacing=1.0, handletextpad=0.5)


def main() -> None:
    common.configure_style()
    fig, axes = plt.subplots(1, 3, figsize=(11.2, 3.8))
    for ax, cfg in zip(axes, PANELS):
        _plot_panel(ax, *cfg)
    handles = [
        plt.Line2D([0], [0], marker='o', linestyle='None', markerfacecolor='none', markeredgecolor=common.GRAY_COLOR, color=common.GRAY_COLOR, label='Methanol data'),
        plt.Line2D([0], [0], marker='s', linestyle='None', markerfacecolor=common.GREEN_COLOR, markeredgecolor=common.GREEN_COLOR, color=common.GREEN_COLOR, label='Ethanol data'),
        plt.Line2D([0], [0], color=common.GRAY_COLOR, linewidth=1.5, label='Methanol model'),
        plt.Line2D([0], [0], color='black', linewidth=1.5, label='Ethanol model'),
    ]
    axes[1].legend(handles=handles, loc='upper center', bbox_to_anchor=(0.5, 0.98), ncol=2, fontsize=8, frameon=False, columnspacing=1.0, handletextpad=0.5)
    fig.suptitle('Molality-based salt mean ionic activity coefficients $\\gamma_{\\pm}^{m,*}$ in methanol and ethanol\nat 298.15 K and 1 bar.', fontsize=11, y=0.995)
    fig.subplots_adjust(left=0.07, right=0.99, bottom=0.17, top=0.80, wspace=0.26)
    common.save_figure(fig, OUTPUT)
    plt.close(fig)

    for cfg, out in zip(PANELS, OUTPUTS):
        common.save_panel_figure(lambda ax, cfg=cfg: _plot_panel(ax, *cfg, include_legend=True), out, figsize=(4.1, 3.9))


if __name__ == '__main__':
    main()




