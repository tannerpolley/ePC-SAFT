from __future__ import annotations

from pathlib import Path
import sys

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import _common as common

OUTPUT = Path(__file__).with_name('figure_6.png')
OUTPUT_METHANOL = Path(__file__).with_name('figure_6_methanol.png')
OUTPUT_ETHANOL = Path(__file__).with_name('figure_6_ethanol.png')
OUTPUTS = [
    Path(__file__).with_name('figure_6a.png'),
    Path(__file__).with_name('figure_6b.png'),
    Path(__file__).with_name('figure_6c.png'),
    Path(__file__).with_name('figure_6d.png'),
]
DATA_ROOT = Path(__file__).with_name('data') / 'G_trans' / 'water'
PANELS = [
    ('a)', 'K+', 'methanol', DATA_ROOT / 'methanol' / 'K.csv', r'$x_{MeOH}$ / -', (0.0, 12.5)),
    ('b)', 'Br-', 'methanol', DATA_ROOT / 'methanol' / 'Br.csv', r'$x_{MeOH}$ / -', (0.0, 12.0)),
    ('c)', 'Na+', 'ethanol', DATA_ROOT / 'ethanol' / 'Na.csv', r'$x_{EtOH}$ / -', (0.0, 20.0)),
    ('d)', 'Cl-', 'ethanol', DATA_ROOT / 'ethanol' / 'Cl.csv', r'$x_{EtOH}$ / -', (0.0, 30.0)),
]


TITLE_ALL = 'Gibbs energies of transfer at infinite dilution from water to water + organic solvent systems\n$\\Delta G_i^{trans,\\infty}$ of different ions at 298.15 K and 1 bar.'
TITLE_METHANOL = 'Gibbs energies of transfer at infinite dilution from water to water + $\\mathbf{methanol}$ systems\n$\\Delta G_i^{trans,\\infty}$ of different ions at 298.15 K and 1 bar.'
TITLE_ETHANOL = 'Gibbs energies of transfer at infinite dilution from water to water + $\\mathbf{ethanol}$ systems\n$\\Delta G_i^{trans,\\infty}$ of different ions at 298.15 K and 1 bar.'


def _load_xy(path: Path):
    fields, rows = common.read_csv_rows(path)
    x_key = fields[0]
    y_key = fields[1]
    xs, ys = [], []
    for row in rows:
        x = common.parse_float(row.get(x_key))
        y = common.parse_float(row.get(y_key))
        if x is None or y is None:
            continue
        xs.append(x)
        ys.append(y)
    return np.asarray(xs, dtype=float), np.asarray(ys, dtype=float)


def _plot_panel(ax, label, ion, organic, csv_path, xlabel, ylim) -> None:
    x_data, y_data = _load_xy(csv_path)
    x_grid = np.linspace(0.0, 1.0, 401)
    y_model = common.transfer_curve('2025_Figiel', ion, organic, x_grid)
    ax.plot(x_grid, y_model, color='black', linewidth=1.5)
    ax.scatter(x_data, y_data, s=24, facecolor=common.LIGHT_GRAY, edgecolor=common.GRAY_COLOR, linewidth=0.8)
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(*ylim)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(r'$\Delta G_i^{trans,\infty,x}$ / kJ mol$^{-1}$')


def _legend_handles():
    return [
        Line2D([0], [0], color='black', linewidth=1.5, label='ePC-SAFT fit'),
        Line2D([0], [0], marker='o', linestyle='None', markersize=5.5, markerfacecolor=common.LIGHT_GRAY, markeredgecolor=common.GRAY_COLOR, label='Literature data'),
    ]


def _save_combined(panels, output: Path, title: str, shape: tuple[int, int], figsize: tuple[float, float], top: float, legend_y: float):
    common.configure_style()
    fig, axes = plt.subplots(*shape, figsize=figsize)
    axes_flat = np.atleast_1d(axes).ravel()
    for ax, cfg in zip(axes_flat, panels):
        _plot_panel(ax, *cfg)
    handles = _legend_handles()
    fig.legend(handles=handles, loc='upper center', ncol=2, bbox_to_anchor=(0.5, legend_y), fontsize=9, frameon=False)
    fig.suptitle(title, fontsize=11, y=0.995)
    fig.subplots_adjust(left=0.10, right=0.98, bottom=0.14, top=top, wspace=0.24, hspace=0.28)
    common.save_figure(fig, output)
    plt.close(fig)


def main() -> None:
    _save_combined(PANELS, OUTPUT, TITLE_ALL, (2, 2), (7.6, 6.3), 0.87, 0.945)
    _save_combined(PANELS[:2], OUTPUT_METHANOL, TITLE_METHANOL, (1, 2), (7.6, 3.8), 0.72, 0.88)
    _save_combined(PANELS[2:], OUTPUT_ETHANOL, TITLE_ETHANOL, (1, 2), (7.6, 3.8), 0.72, 0.88)

    for cfg, out in zip(PANELS, OUTPUTS):
        common.save_panel_figure(lambda ax, cfg=cfg: _plot_panel(ax, *cfg), out, figsize=(3.8, 3.3))


if __name__ == '__main__':
    main()

