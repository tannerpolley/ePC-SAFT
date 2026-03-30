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

from scripts._env import require_pcsaft_install

require_pcsaft_install()

from plot_figure_7 import (  # type: ignore
    PANELS,
    T_REF,
    P_REF,
    _build_params,
    _molality_for_salt_mole_fraction,
    _molality_to_species_molefraction,
    _resolve_pair_key,
    _species_for_combo,
    SOLVENT_MW,
)
from _plot_common import configure_style, save_figure  # type: ignore
from pcsaft import pcsaft_den, pcsaft_miac

matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT_DIR = THIS_DIR

SCAN_BOUNDS = {
    'ethanol': np.linspace(18.0, 32.0, 57),
    'methanol': np.linspace(26.0, 40.0, 57),
}


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
                y_val = y_m * (1.0 + SOLVENT_MW[solvent] * m_val * 2.0)
            except ValueError:
                y_val = math.nan
        if not math.isfinite(x_val):
            try:
                m_val = float(str(row.get('molality', '')).strip())
                x_val = float(m_val / ((1.0 / SOLVENT_MW[solvent]) + m_val))
            except ValueError:
                x_val = math.nan
        if math.isfinite(x_val) and math.isfinite(y_val):
            xs.append(float(x_val))
            ys.append(float(y_val))
    x_arr = np.asarray(xs, dtype=float)
    y_arr = np.asarray(ys, dtype=float)
    order = np.argsort(x_arr)
    return x_arr[order], y_arr[order]


def _curve_with_solvent_dielc(salt: str, solvent: str, x_points: np.ndarray, solvent_eps: float) -> np.ndarray:
    params = _build_params('2020_Bulow', salt, solvent)
    dielc = np.asarray(params['dielc'], dtype=float).copy()
    dielc[2] = float(solvent_eps)
    params['dielc'] = dielc
    species = _species_for_combo(salt, solvent)
    pair_key = None
    y = np.empty_like(x_points, dtype=float)
    m_grid = _molality_for_salt_mole_fraction(x_points, solvent)
    for idx, m in enumerate(m_grid):
        x = _molality_to_species_molefraction(float(max(m, 1.0e-12)), salt, solvent)
        rho = pcsaft_den(T_REF, P_REF, x, params, phase='liq')
        vals = pcsaft_miac(T_REF, rho, x, params, species=species)
        if pair_key is None:
            pair_key = _resolve_pair_key(vals, salt)
        y[idx] = float(vals[pair_key])
    return y


def _metrics(y_model: np.ndarray, y_ref: np.ndarray) -> dict[str, float]:
    delta = np.asarray(y_model, dtype=float) - np.asarray(y_ref, dtype=float)
    return {
        'rmse': float(np.sqrt(np.mean(delta ** 2))),
        'mae': float(np.mean(np.abs(delta))),
        'bias': float(np.mean(delta)),
        'max_abs': float(np.max(np.abs(delta))),
    }


def _write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        return
    with path.open('w', newline='', encoding='utf-8') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def _plot_scan(panel: dict[str, object], x_data: np.ndarray, y_data: np.ndarray, x_grid: np.ndarray, y_base: np.ndarray, y_best: np.ndarray, base_eps: float, best_eps: float, out_path: Path) -> None:
    configure_style()
    fig, ax = plt.subplots(figsize=(5.4, 4.2))
    ax.scatter(x_data, y_data, s=28, facecolor='white', edgecolor='black', linewidth=0.8, label='Paper points', zorder=5)
    ax.plot(x_grid, y_base, color='#228b22', linewidth=1.9, label=f'Advanced baseline ($\\varepsilon_{{solv}}={base_eps:.2f}$)', zorder=3)
    ax.plot(x_grid, y_best, color='#b22222', linestyle='--', linewidth=1.9, label=f'Best scan ($\\varepsilon_{{solv}}={best_eps:.2f}$)', zorder=4)
    ax.set_xlim(*panel['xlim'])
    ax.set_ylim(*panel['ylim'])
    ax.set_xlabel(r'salt mole fraction, $x_{salt}$')
    ax.set_ylabel(r'$\gamma_{\pm}^{*}$ / -')
    ax.set_title(f"Figure {panel['id']} solvent dielectric sensitivity")
    ax.grid(True, alpha=0.25)
    ax.legend(loc='best')
    save_figure(fig, out_path)
    plt.close(fig)


def main() -> None:
    summary_rows: list[dict[str, object]] = []
    for panel in PANELS:
        panel_id = str(panel['id'])
        salt = str(panel['salt'])
        solvent = str(panel['solvent'])
        x_data, y_data = _read_paper_points(panel)
        x_grid = np.linspace(float(panel['xlim'][0]), float(panel['xlim'][1]), 801)
        base_params = _build_params('2020_Bulow', salt, solvent)
        base_eps = float(np.asarray(base_params['dielc'], dtype=float)[2])
        scan_eps = np.asarray(SCAN_BOUNDS[solvent], dtype=float)
        scan_rows: list[dict[str, object]] = []
        best_row = None
        y_base_grid = _curve_with_solvent_dielc(salt, solvent, x_grid, base_eps)
        for eps in scan_eps:
            y_model = _curve_with_solvent_dielc(salt, solvent, x_data, float(eps))
            stats = _metrics(y_model, y_data)
            row = {
                'panel': panel_id,
                'salt': salt,
                'solvent': solvent,
                'base_solvent_dielc': base_eps,
                'trial_solvent_dielc': float(eps),
                'delta_eps': float(eps - base_eps),
                'rmse': stats['rmse'],
                'mae': stats['mae'],
                'bias': stats['bias'],
                'max_abs': stats['max_abs'],
            }
            scan_rows.append(row)
            if best_row is None or row['rmse'] < best_row['rmse']:
                best_row = row
        csv_path = OUT_DIR / f'figure_{panel_id}_relperm_scan.csv'
        _write_csv(csv_path, scan_rows)
        if best_row is None:
            continue
        y_best_grid = _curve_with_solvent_dielc(salt, solvent, x_grid, float(best_row['trial_solvent_dielc']))
        plot_path = OUT_DIR / f'figure_{panel_id}_relperm_scan.png'
        _plot_scan(panel, x_data, y_data, x_grid, y_base_grid, y_best_grid, base_eps, float(best_row['trial_solvent_dielc']), plot_path)
        baseline_metrics = _metrics(_curve_with_solvent_dielc(salt, solvent, x_data, base_eps), y_data)
        summary_rows.append({
            'panel': panel_id,
            'salt': salt,
            'solvent': solvent,
            'base_solvent_dielc': base_eps,
            'best_solvent_dielc': float(best_row['trial_solvent_dielc']),
            'baseline_rmse': baseline_metrics['rmse'],
            'best_rmse': best_row['rmse'],
            'baseline_bias': baseline_metrics['bias'],
            'best_bias': best_row['bias'],
            'rmse_improvement': float(baseline_metrics['rmse'] - best_row['rmse']),
            'scan_csv': csv_path.name,
            'scan_plot': plot_path.name,
        })
    _write_csv(OUT_DIR / 'figure7_relperm_sensitivity_summary.csv', summary_rows)


if __name__ == '__main__':
    main()
