from __future__ import annotations

import csv
import math
import sys
from pathlib import Path

import matplotlib
import numpy as np

THIS_DIR = Path(__file__).resolve().parent
FIGURE_DIR = THIS_DIR.parent
ANALYSIS_ROOT = FIGURE_DIR.parent
REPO_ROOT = FIGURE_DIR.parents[3]
if str(FIGURE_DIR) not in sys.path:
    sys.path.insert(0, str(FIGURE_DIR))
if str(ANALYSIS_ROOT) not in sys.path:
    sys.path.insert(0, str(ANALYSIS_ROOT))
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts._env import require_epcsaft_install
from scripts.plot_outputs import paper_validation_dir

require_epcsaft_install()

from plot_figure_7 import (  # type: ignore
    PANELS,
    T_REF,
    P_REF,
    _build_params,
    _molality_for_salt_mole_fraction,
    _molality_to_species_molefraction,
    _read_dataset,
    _resolve_pair_key,
    _species_for_combo,
)
from _plot_common import configure_style, save_figure  # type: ignore
from scripts._epcsaft_oop import epcsaft_activity_coefficient, epcsaft_density

matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT_DIR = paper_validation_dir(__file__)

PAIR_SCANS = [
    {
        'panel_id': '7a',
        'pair_label': 'Li-Br',
        'runtime_pair_label': 'Li-Br',
        'pair_indices': (0, 1),
        'k_values': np.linspace(0.300, 1.100, 41),
        'slug': 'figure7a_li_br',
        'note': 'Direct salt-ion unlike pair in LiBr/ethanol.',
    },
    {
        'panel_id': '7a',
        'pair_label': 'Li-H2O (requested)',
        'runtime_pair_label': 'Li-Ethanol',
        'pair_indices': (0, 2),
        'k_values': np.linspace(-0.500, 0.300, 41),
        'slug': 'figure7a_li_solvent',
        'note': 'No H2O is present in Figure 7a; Li-Ethanol is the actual active solvent pair.',
    },
    {
        'panel_id': '7c',
        'pair_label': 'Li-H2O (requested)',
        'runtime_pair_label': 'Li-Methanol',
        'pair_indices': (0, 2),
        'k_values': np.linspace(-0.500, 0.300, 41),
        'slug': 'figure7c_li_solvent',
        'note': 'No H2O is present in Figure 7c; Li-Methanol is the actual active solvent pair.',
    },
]


def _panel_map() -> dict[str, dict[str, object]]:
    return {str(panel['id']): panel for panel in PANELS}


def _read_paper_points(panel: dict[str, object]) -> tuple[np.ndarray, np.ndarray]:
    path = Path(panel['data'])
    solvent = str(panel['solvent'])
    xs: list[float] = []
    ys: list[float] = []
    with path.open('r', newline='', encoding='utf-8-sig') as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)
    use_bulow = any(str(row.get('source', '')).strip() == 'Bulow 2020' for row in rows)
    for row in rows:
        if use_bulow and str(row.get('source', '')).strip() != 'Bulow 2020':
            continue
        try:
            x_val = float(str(row.get('mole_fraction', '')).strip())
        except ValueError:
            x_val = math.nan
        try:
            y_val = float(str(row.get('miac', '')).strip())
        except ValueError:
            y_val = math.nan
        if not math.isfinite(y_val):
            try:
                m_val = float(str(row.get('molality', '')).strip())
                y_m = float(str(row.get('miac_m', row.get('gamma', ''))).strip())
                mw = 32.04e-3 if solvent == 'methanol' else 46.068e-3
                y_val = y_m * (1.0 + mw * m_val * 2.0)
            except ValueError:
                y_val = math.nan
        if not math.isfinite(x_val):
            try:
                m_val = float(str(row.get('molality', '')).strip())
                x_val = float(_molality_for_salt_mole_fraction(np.asarray([m_val], dtype=float), solvent)[0])
            except ValueError:
                x_val = math.nan
        if math.isfinite(x_val) and math.isfinite(y_val):
            xs.append(float(x_val))
            ys.append(float(y_val))
    x_arr = np.asarray(xs, dtype=float)
    y_arr = np.asarray(ys, dtype=float)
    order = np.argsort(x_arr)
    return x_arr[order], y_arr[order]


def _curve_with_override(dataset: str, salt: str, solvent: str, x_points: np.ndarray, pair_indices: tuple[int, int], k_value: float) -> np.ndarray:
    params = _build_params(dataset, salt, solvent)
    i, j = pair_indices
    k_mat = np.asarray(params['k_ij'], dtype=float).copy()
    k_mat[i, j] = float(k_value)
    k_mat[j, i] = float(k_value)
    params['k_ij'] = k_mat
    species = _species_for_combo(salt, solvent)
    pair_key = None
    y = np.empty_like(x_points, dtype=float)
    m_grid = _molality_for_salt_mole_fraction(x_points, solvent)
    for idx, m in enumerate(m_grid):
        x = _molality_to_species_molefraction(float(max(m, 1.0e-12)), salt, solvent)
        rho = epcsaft_density(T_REF, P_REF, x, params, phase='liq')
        vals = epcsaft_activity_coefficient(T_REF, rho, x, params, species=species, mean_ionic_form=True, basis="mole")
        if pair_key is None:
            pair_key = _resolve_pair_key(vals, salt)
        y[idx] = float(vals[pair_key])
    return y


def _metrics(y_model: np.ndarray, y_ref: np.ndarray) -> dict[str, float]:
    delta = np.asarray(y_model, dtype=float) - np.asarray(y_ref, dtype=float)
    return {
        'rmse': float(np.sqrt(np.mean(delta**2))),
        'mae': float(np.mean(np.abs(delta))),
        'bias': float(np.mean(delta)),
        'max_abs': float(np.max(np.abs(delta))),
    }


def _visible_mask(panel: dict[str, object], x_vals: np.ndarray, y_vals: np.ndarray) -> np.ndarray:
    x_min, x_max = panel['xlim']
    _y_min, y_max = panel['ylim']
    return (x_vals >= float(x_min) - 1e-12) & (x_vals <= float(x_max) + 1e-12) & (y_vals <= float(y_max) + 1e-12)


def _write_scan_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        return
    fields = list(rows[0].keys())
    with path.open('w', newline='', encoding='utf-8') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def _plot_scan(panel: dict[str, object], pair_label: str, runtime_pair_label: str, x_data: np.ndarray, y_data: np.ndarray, x_grid: np.ndarray, y_base: np.ndarray, y_best: np.ndarray, best_k: float, base_k: float, out_path: Path) -> None:
    configure_style()
    fig, ax = plt.subplots(figsize=(5.4, 4.2))
    ax.scatter(x_data, y_data, s=28, facecolor='white', edgecolor='black', linewidth=0.8, label='Paper points', zorder=5)
    ax.plot(x_grid, y_base, color='#228b22', linewidth=1.9, label=f'Advanced baseline ({runtime_pair_label}={base_k:.3f})', zorder=3)
    ax.plot(x_grid, y_best, color='#b22222', linewidth=1.9, linestyle='--', label=f'Best scan ({runtime_pair_label}={best_k:.3f})', zorder=4)
    ax.set_xlim(*panel['xlim'])
    ax.set_ylim(*panel['ylim'])
    ax.set_xlabel(r'salt mole fraction, $x_{salt}$')
    ax.set_ylabel(r'$\gamma_{\pm}^{*}$ / -')
    ax.set_title(f"Figure {panel['id']} sensitivity: {pair_label}")
    ax.grid(True, alpha=0.25)
    ax.legend(loc='best')
    save_figure(fig, out_path)
    plt.close(fig)


def main() -> None:
    panel_map = _panel_map()
    summary_rows: list[dict[str, object]] = []
    for cfg in PAIR_SCANS:
        panel = panel_map[cfg['panel_id']]
        salt = str(panel['salt'])
        solvent = str(panel['solvent'])
        x_paper, y_paper = _read_paper_points(panel)
        visible_mask = _visible_mask(panel, x_paper, y_paper)
        x_eval = x_paper[visible_mask]
        y_eval = y_paper[visible_mask]
        x_grid = np.linspace(float(panel['xlim'][0]), float(panel['xlim'][1]), 801)
        base_params = _build_params('2020_Bulow', salt, solvent)
        i, j = cfg['pair_indices']
        base_k = float(np.asarray(base_params['k_ij'], dtype=float)[i, j])

        scan_rows: list[dict[str, object]] = []
        best = None
        y_base_grid = _curve_with_override('2020_Bulow', salt, solvent, x_grid, cfg['pair_indices'], base_k)
        for k_value in np.asarray(cfg['k_values'], dtype=float):
            y_model_full = _curve_with_override('2020_Bulow', salt, solvent, x_paper, cfg['pair_indices'], float(k_value))
            y_model = y_model_full[visible_mask]
            stats = _metrics(y_model, y_eval)
            stats_full = _metrics(y_model_full, y_paper)
            row = {
                'panel': cfg['panel_id'],
                'salt': salt,
                'solvent': solvent,
                'requested_pair': cfg['pair_label'],
                'runtime_pair': cfg['runtime_pair_label'],
                'k_value': float(k_value),
                'rmse': stats['rmse'],
                'mae': stats['mae'],
                'bias': stats['bias'],
                'max_abs': stats['max_abs'],
                'rmse_full': stats_full['rmse'],
                'mae_full': stats_full['mae'],
                'bias_full': stats_full['bias'],
                'max_abs_full': stats_full['max_abs'],
                'delta_k_from_baseline': float(k_value - base_k),
            }
            scan_rows.append(row)
            if best is None or row['rmse'] < best['rmse']:
                best = row
        if best is None:
            continue
        csv_path = OUT_DIR / f"{cfg['slug']}_scan.csv"
        _write_scan_csv(csv_path, scan_rows)
        y_best_grid = _curve_with_override('2020_Bulow', salt, solvent, x_grid, cfg['pair_indices'], float(best['k_value']))
        plot_path = OUT_DIR / f"{cfg['slug']}_scan.png"
        _plot_scan(panel, str(cfg['pair_label']), str(cfg['runtime_pair_label']), x_paper[visible_mask], y_paper[visible_mask], x_grid, y_base_grid, y_best_grid, float(best['k_value']), base_k, plot_path)
        y_base_eval_full = _curve_with_override('2020_Bulow', salt, solvent, x_paper, cfg['pair_indices'], base_k)
        baseline_metrics = _metrics(y_base_eval_full[visible_mask], y_eval)
        summary_rows.append({
            'panel': cfg['panel_id'],
            'salt': salt,
            'solvent': solvent,
            'requested_pair': cfg['pair_label'],
            'runtime_pair': cfg['runtime_pair_label'],
            'baseline_k': base_k,
            'best_k_in_scan': float(best['k_value']),
            'baseline_rmse': baseline_metrics['rmse'],
            'best_rmse': best['rmse'],
            'baseline_bias': baseline_metrics['bias'],
            'best_bias': best['bias'],
            'rmse_improvement': float(baseline_metrics['rmse'] - best['rmse']),
            'note': cfg['note'],
            'scan_csv': csv_path.name,
            'scan_plot': plot_path.name,
        })
    _write_scan_csv(OUT_DIR / 'figure7_kij_sensitivity_summary.csv', summary_rows)


if __name__ == '__main__':
    main()



