from __future__ import annotations

import csv
import math
import sys
from pathlib import Path

import matplotlib
import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[6]
DIAG_DIR = Path(__file__).resolve().parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
if str(DIAG_DIR) not in sys.path:
    sys.path.insert(0, str(DIAG_DIR))
FIG6_DIR = DIAG_DIR.parent
ANALYSIS_ROOT = FIG6_DIR.parent.parent
if str(ANALYSIS_ROOT) not in sys.path:
    sys.path.insert(0, str(ANALYSIS_ROOT))

from scripts.plot_outputs import paper_validation_path
from figure6b_libr_ethanol_contributions import _calc_ln_miac_contributions, _build_params, _salt_mole_fraction_from_molality
from _plot_common import configure_style, save_figure

OUTPUT_ROOT = paper_validation_path(DIAG_DIR, 'output')
OUTPUT_DATA = OUTPUT_ROOT / 'data'
OUTPUT_PLOTS = OUTPUT_ROOT / 'plots'

matplotlib.use('Agg')
import matplotlib.pyplot as plt

CURVE_PATH = REPO_ROOT / 'scripts' / 'paper_validation' / '2020_Bulow_analysis' / 'figure_6' / 'figure_6b' / 'data' / 'Figure6b_curves.csv'
SERIES_MAP = {
    'total': 'data',
    'hc': 'hc',
    'disp': 'disp',
    'assoc': 'assoc',
    'dh': 'DH',
    'born': 'Born',
}
PLOT_LABELS = {
    'total': 'Total',
    'dh': 'DH',
}
COLORS = {
    'total': 'green',
    'dh': 'tab:blue',
}


def _read_digitized_series(path: Path) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    with path.open('r', newline='', encoding='utf-8-sig') as handle:
        reader = csv.reader(handle)
        header = next(reader)
        rows = list(reader)
    out: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    ncols = len(header)
    for i in range(0, ncols, 2):
        if i + 1 >= ncols:
            continue
        y_name = str(header[i + 1]).strip()
        if not y_name:
            continue
        xs: list[float] = []
        ys: list[float] = []
        for row in rows:
            if i >= len(row) or i + 1 >= len(row):
                continue
            try:
                x_val = float(str(row[i]).strip())
                y_val = float(str(row[i + 1]).strip())
            except ValueError:
                continue
            if math.isfinite(x_val) and math.isfinite(y_val):
                xs.append(x_val)
                ys.append(y_val)
        if xs:
            x_arr = np.asarray(xs, dtype=float)
            y_arr = np.asarray(ys, dtype=float)
            order = np.argsort(x_arr)
            out[y_name] = (x_arr[order], y_arr[order])
    return out


def _metrics(y_model: np.ndarray, y_ref: np.ndarray) -> dict[str, float]:
    delta = np.asarray(y_model, dtype=float) - np.asarray(y_ref, dtype=float)
    return {
        'rmse': float(np.sqrt(np.mean(delta * delta))),
        'mae': float(np.mean(np.abs(delta))),
        'bias': float(np.mean(delta)),
        'max_abs': float(np.max(np.abs(delta))),
    }


def _write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', newline='', encoding='utf-8') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    configure_style()
    OUTPUT_DATA.mkdir(parents=True, exist_ok=True)
    OUTPUT_PLOTS.mkdir(parents=True, exist_ok=True)

    digitized = _read_digitized_series(CURVE_PATH)
    x_total, y_total = digitized['data']
    max_x = float(np.max(x_total))
    max_m = max_x * (1.0 / 46.068e-3) / max(1.0 - max_x, 1.0e-12)
    m_grid = np.linspace(0.0, max_m, 1201)
    x_grid = _salt_mole_fraction_from_molality(m_grid)

    configs = {
        'baseline_rule1': {},
        'salt_basis_rule7': {'elec_model': {'rel_perm': {'rule': 7}}},
    }
    curves_by_cfg: dict[str, dict[str, np.ndarray]] = {}
    metric_rows: list[dict[str, object]] = []

    for cfg_name, user_opts in configs.items():
        params = _build_params(user_options=user_opts)
        curves = _calc_ln_miac_contributions(m_grid, params, method='mu')
        curves_by_cfg[cfg_name] = curves
        for model_key, paper_key in SERIES_MAP.items():
            x_ref, y_ref = digitized[paper_key]
            y_model = np.interp(x_ref, x_grid, curves[model_key])
            stats = _metrics(y_model, y_ref)
            metric_rows.append({
                'config': cfg_name,
                'series': model_key,
                'paper_series': paper_key,
                'rmse': stats['rmse'],
                'mae': stats['mae'],
                'bias': stats['bias'],
                'max_abs': stats['max_abs'],
            })

    _write_csv(OUTPUT_DATA / 'figure6b_relperm_rule7_metrics.csv', metric_rows)

    fig, axes = plt.subplots(1, 2, figsize=(10.4, 4.4), sharex=True)
    for ax, series in zip(axes, ('total', 'dh')):
        paper_key = SERIES_MAP[series]
        x_ref, y_ref = digitized[paper_key]
        ax.scatter(x_ref, y_ref, s=24, facecolor='white', edgecolor='black', linewidth=0.8, label='Paper points', zorder=6)
        for cfg_name, style in (
            ('baseline_rule1', {'linestyle': '-', 'label': 'rule 1'}),
            ('salt_basis_rule7', {'linestyle': '--', 'label': 'rule 7 (salt basis)'}),
        ):
            ax.plot(x_grid, curves_by_cfg[cfg_name][series], color=COLORS[series], linewidth=2.0, linestyle=style['linestyle'], label=f"{PLOT_LABELS[series]} {style['label']}", zorder=4)
        ax.set_title(f"Figure 6b {PLOT_LABELS[series]}")
        ax.set_xlabel(r"salt mole fraction, $x_{salt}$")
        ax.grid(True, alpha=0.25)
    axes[0].set_ylabel(r"$\ln(\gamma_{\pm}^{*})$")
    axes[0].set_xlim(0.0, 0.2)
    axes[0].set_ylim(-3.0, 4.0)
    axes[1].set_xlim(0.0, 0.2)
    axes[1].set_ylim(-3.0, 1.0)
    axes[1].legend(loc='best', fontsize=8)
    save_figure(fig, OUTPUT_PLOTS / 'figure6b_relperm_rule7_comparison.png')
    plt.close(fig)


if __name__ == '__main__':
    main()





