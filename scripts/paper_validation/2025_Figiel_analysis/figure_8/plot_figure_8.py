from __future__ import annotations

from pathlib import Path
import sys

import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import _common as common

OUTPUT = Path(__file__).with_name('figure_8.png')
OUTPUTS = [
    Path(__file__).with_name('figure_8a.png'),
    Path(__file__).with_name('figure_8b.png'),
    Path(__file__).with_name('figure_8c.png'),
]
DATA_ROOT = Path(__file__).with_name('data')
PANELS = [
    ('a)', 'LiBr', 5.0, 3.0),
    ('b)', 'NaI', 1.5, 1.125),
    ('c)', 'NaBr', 1.5, 1.125),
]
CURVE_MAX_OVERRIDES = {
    ('NaI', 'methanol'): 0.77,
    ('NaBr', 'ethanol'): 0.195,
}


def _plot_panel(ax, label, salt, m_max, y_max, include_legend: bool = False):
    methanol_data = common.read_miac_dataset(DATA_ROOT / 'methanol' / f'methanol-{salt}.csv', 'methanol')
    ethanol_data = common.read_miac_dataset(DATA_ROOT / 'ethanol' / f'ethanol-{salt}.csv', 'ethanol')
    meoh_curve_max = CURVE_MAX_OVERRIDES.get((salt, 'methanol'), m_max)
    etoh_curve_max = CURVE_MAX_OVERRIDES.get((salt, 'ethanol'), m_max)
    m_grid_meoh, y_meoh = common.mean_ionic_activity_curve('2025_Figiel', salt, 'methanol', {'methanol': 1.0}, meoh_curve_max, points=600)
    m_grid_etoh, y_etoh = common.mean_ionic_activity_curve('2025_Figiel', salt, 'ethanol', {'ethanol': 1.0}, etoh_curve_max, points=600)
    ax.scatter([r['molality'] for r in methanol_data], [r['miac_m'] for r in methanol_data], s=24, facecolor='none', edgecolor=common.GRAY_COLOR, linewidth=0.9, label='Methanol data')
    ax.scatter([r['molality'] for r in ethanol_data], [r['miac_m'] for r in ethanol_data], s=24, marker='s', facecolor=common.GREEN_COLOR, edgecolor=common.GREEN_COLOR, linewidth=0.8, label='Ethanol data')
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




