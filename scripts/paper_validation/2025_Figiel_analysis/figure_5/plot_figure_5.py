from __future__ import annotations

from pathlib import Path
import sys

import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import _common as common

OUTPUT = Path(__file__).with_name('figure_5.png')
OUTPUT_A = Path(__file__).with_name('figure_5a.png')
OUTPUT_B = Path(__file__).with_name('figure_5b.png')
DATA_ROOT = Path(__file__).with_name('data') / 'water'

SERIES = {
    'Li': {'color': common.ORGANIC_COLOR, 'marker': '^'},
    'Na': {'color': common.GREEN_COLOR, 'marker': 's'},
    'K': {'color': common.GRAY_COLOR, 'marker': 'o'},
}
PANELS = [
    ('a)', ['LiCl', 'NaCl', 'KCl'], 'Cl$^-$ salts in water'),
    ('b)', ['LiBr', 'NaBr', 'KBr'], 'Br$^-$ salts in water'),
]


def _plot_panel(ax, label, salts, title):
    for salt in salts:
        cation = salt.split('C')[0] if 'Cl' in salt else salt.split('B')[0]
        style = SERIES[cation]
        data = common.read_miac_dataset(DATA_ROOT / f'water-{salt}.csv', 'water')
        x_data = [row['molality'] for row in data]
        y_data = [row['miac_m'] for row in data]
        m_grid, y_model = common.mean_ionic_activity_curve('figiel_2025', salt, 'water', {'water': 1.0}, 6.0, points=600)
        ax.plot(m_grid, y_model, color=style['color'], linewidth=1.8)
        ax.scatter(x_data, y_data, marker=style['marker'], s=26, facecolor='none', edgecolor=style['color'], linewidth=1.0, label=salt)
    ax.set_title(title, fontsize=10)
    ax.set_xlim(0.0, 6.0)
    ax.set_ylim(0.4, 4.6)
    ax.set_xlabel(r'$\bar{m}_{salt}$ / mol kg$^{-1}$')
    ax.set_ylabel(r'$\gamma_{\pm}^{m,*}$ / -')


def _plot_panel_a(ax) -> None:
    _plot_panel(ax, *PANELS[0])


def _plot_panel_b(ax) -> None:
    _plot_panel(ax, *PANELS[1])


def main() -> None:
    common.configure_style()
    fig, axes = plt.subplots(1, 2, figsize=(8.2, 3.9), sharey=True)
    _plot_panel_a(axes[0])
    _plot_panel_b(axes[1])
    axes[0].legend(loc='upper left', fontsize=8, frameon=False)
    axes[1].legend(loc='upper left', fontsize=8, frameon=False)
    fig.suptitle('Molality-based mean ionic activity coefficients $\\gamma_{\\pm}^{m,*}$ of alkali halides in water\nat 298.15 K and 1 bar.', fontsize=11, y=0.99)
    fig.subplots_adjust(left=0.08, right=0.98, bottom=0.16, top=0.84, wspace=0.16)
    common.save_figure(fig, OUTPUT)
    plt.close(fig)

    handles_a, labels_a = [], []
    tmp_fig, tmp_ax = plt.subplots()
    _plot_panel_a(tmp_ax)
    handles_a, labels_a = tmp_ax.get_legend_handles_labels()
    plt.close(tmp_fig)
    tmp_fig, tmp_ax = plt.subplots()
    _plot_panel_b(tmp_ax)
    handles_b, labels_b = tmp_ax.get_legend_handles_labels()
    plt.close(tmp_fig)

    common.save_panel_figure(_plot_panel_a, OUTPUT_A, figsize=(4.1, 3.6), legend_handles=handles_a, legend_labels=labels_a, legend_kwargs={'loc': 'upper left', 'fontsize': 8})
    common.save_panel_figure(_plot_panel_b, OUTPUT_B, figsize=(4.1, 3.6), legend_handles=handles_b, legend_labels=labels_b, legend_kwargs={'loc': 'upper left', 'fontsize': 8})


if __name__ == '__main__':
    main()

