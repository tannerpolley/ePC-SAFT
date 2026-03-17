from __future__ import annotations

from pathlib import Path
from copy import deepcopy
import csv
import math
import sys

import matplotlib
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
REPO_ROOT = ROOT.parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from _common import build_params
from pcsaft import pcsaft_den, pcsaft_lnfugcoef_terms

matplotlib.use('Agg')
import matplotlib.pyplot as plt

T_REF = 298.15
P_REF = 1.0e5
CATION = 'C4mim+'
ANION = 'BF4-'
SPECIES = ['H2O', CATION, ANION]
WATER_INDEX = 0
CATION_INDEX = 1
ANION_INDEX = 2
THIS_DIR = Path(__file__).resolve().parent
FIG2A_CSV = THIS_DIR / '2019_Figure2a_curves.csv'
FIG2B_CSV = THIS_DIR / '2019_Figure2b_curves.csv'
SERIES_META = [
    ('water', r'constant $\varepsilon_{water}$', 'tab:blue'),
    ('mixed', r'ePC-SAFT, $\varepsilon(x_{IL})$', 'green'),
    ('il', r'constant $\varepsilon_{IL}$', 'tab:orange'),
]
KB = 1.380649e-23
E_CHRG = 1.602176634e-19
EPS0 = 8.854187817e-22
N_AV = 6.02214076e23
PI = math.pi


def _composition_from_x_il(x_il: float) -> np.ndarray:
    x_il_clamped = min(max(float(x_il), 1.0e-8), 1.0 - 1.0e-8)
    return np.asarray([1.0 - x_il_clamped, 0.5 * x_il_clamped, 0.5 * x_il_clamped], dtype=float)


def _figure2_params(model_mode: str) -> dict:
    params = build_params(SPECIES, T_REF, use_kij=False, model_mode=model_mode)
    if model_mode == 'epc':
        params = deepcopy(params)
        elec = dict(params['elec_model'])
        rel_perm = dict(elec['rel_perm'])
        rel_perm['rule'] = 7
        rel_perm['differential_mode'] = 'analytical'
        elec['rel_perm'] = rel_perm
        params['elec_model'] = elec
    return params


def _ion_diameter_angstrom(i: int, params: dict) -> float:
    sigma_i = float(np.asarray(params['s'], dtype=float)[i])
    z_i = float(np.asarray(params['z'], dtype=float)[i])
    if abs(z_i) <= 1.0e-12:
        return sigma_i
    mode = int(params['elec_model']['DH_model']['d_ion_mode'])
    epsilon_i = float(np.asarray(params['e'], dtype=float)[i])
    if mode == 0:
        return sigma_i
    if mode == 1:
        return 0.88 * sigma_i
    if mode == 2:
        return sigma_i * (1.0 - 0.12 * math.exp(-3.0 * epsilon_i / T_REF))
    raise ValueError(f'Unsupported d_ion_mode={mode}')


def _manual_eps(model_mode: str, x: np.ndarray, params: dict) -> float:
    dielc = np.asarray(params['dielc'], dtype=float)
    eps_water = float(dielc[WATER_INDEX])
    eps_il = float(dielc[CATION_INDEX])
    if model_mode == 'orig_water':
        return eps_water
    if model_mode == 'orig_il':
        return eps_il
    x_water = float(x[WATER_INDEX])
    return eps_water * x_water + eps_il * (1.0 - x_water)


def _manual_aion_eq8(model_mode: str, x_il: float) -> float:
    params = _figure2_params(model_mode)
    x = _composition_from_x_il(float(x_il))
    rho = float(pcsaft_den(T_REF, P_REF, x, params, phase='liq'))
    z = np.asarray(params['z'], dtype=float)
    eps = _manual_eps(model_mode, x, params)
    den = rho * N_AV / 1.0e30
    qsum = float(np.sum((z * z) * x))
    if eps <= 0.0 or qsum <= 0.0:
        return float('nan')
    kappa = math.sqrt(den * E_CHRG * E_CHRG / (KB * T_REF * eps * EPS0) * qsum)
    if not math.isfinite(kappa) or kappa <= 0.0:
        return float('nan')
    d = np.asarray([_ion_diameter_angstrom(i, params) for i in range(len(z))], dtype=float)
    chi = np.zeros_like(z, dtype=float)
    for i in range(len(z)):
        ka = kappa * d[i]
        if ka == 0.0 or not math.isfinite(ka):
            continue
        chi[i] = 3.0 / (ka ** 3) * (1.5 + math.log(1.0 + ka) - 2.0 * (1.0 + ka) + 0.5 * (1.0 + ka) ** 2)
    s_term = float(np.sum(x * (z * z) * chi))
    k0 = E_CHRG * E_CHRG / (12.0 * PI * KB * T_REF * EPS0)
    return -k0 * kappa / eps * s_term


def _manual_dadx_binary_path(model_mode: str, x_il: float) -> float:
    x0 = min(max(float(x_il), 1.0e-6), 1.0 - 1.0e-6)
    h = min(1.0e-5, 0.25 * min(x0, 1.0 - x0))
    if h <= 0.0:
        return float('nan')
    ap = _manual_aion_eq8(model_mode, x0 + h)
    am = _manual_aion_eq8(model_mode, x0 - h)
    return 2.0 * (ap - am) / (2.0 * h)


def _curve_for_mode(model_mode: str, x_il_grid: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    params = _figure2_params(model_mode)
    dadx_cat = np.full_like(x_il_grid, np.nan, dtype=float)
    dadx_manual = np.full_like(x_il_grid, np.nan, dtype=float)
    lnfug_cat = np.full_like(x_il_grid, np.nan, dtype=float)
    for idx, x_il in enumerate(x_il_grid):
        x = _composition_from_x_il(float(x_il))
        try:
            rho = pcsaft_den(T_REF, P_REF, x, params, phase='liq')
            terms = pcsaft_lnfugcoef_terms(T_REF, rho, x, params)
            dadx_ion = np.asarray(terms['dadx_ion'], dtype=float)
            dadx_cat[idx] = float(dadx_ion[CATION_INDEX])
            dadx_manual[idx] = _manual_dadx_binary_path(model_mode, float(x_il))
            lnfug_cat[idx] = float(np.asarray(terms['lnfugcoef'], dtype=float)[CATION_INDEX])
        except Exception:
            continue
    return dadx_cat, dadx_manual, lnfug_cat


def _read_digitized_wide(path: Path) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    out: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    with path.open('r', newline='', encoding='utf-8-sig') as handle:
        reader = csv.reader(handle)
        rows = [row for row in reader if row]
    for idx, key in enumerate(('water', 'mixed', 'il')):
        x_col = 2 * idx
        y_col = 2 * idx + 1
        xs: list[float] = []
        ys: list[float] = []
        for row in rows[1:]:
            if x_col >= len(row) or y_col >= len(row):
                continue
            x_raw = str(row[x_col]).strip()
            y_raw = str(row[y_col]).strip()
            if not x_raw or not y_raw:
                continue
            try:
                x_val = float(x_raw)
                y_val = float(y_raw)
            except ValueError:
                continue
            if math.isfinite(x_val) and math.isfinite(y_val):
                xs.append(x_val)
                ys.append(y_val)
        out[key] = (np.asarray(xs, dtype=float), np.asarray(ys, dtype=float))
    return out


def _plot_single(
    ax: plt.Axes,
    x_il: np.ndarray,
    water_curve: np.ndarray,
    il_curve: np.ndarray,
    epc_curve: np.ndarray,
    digitized: dict[str, tuple[np.ndarray, np.ndarray]],
    *,
    ylabel: str,
    title: str,
    ylim: tuple[float, float],
    water_manual: np.ndarray | None = None,
    il_manual: np.ndarray | None = None,
    epc_manual: np.ndarray | None = None,
) -> None:
    curve_map = {
        'water': water_curve,
        'mixed': epc_curve,
        'il': il_curve,
    }
    manual_map = {
        'water': water_manual,
        'mixed': epc_manual,
        'il': il_manual,
    }
    for key, label, color in SERIES_META:
        ax.plot(x_il, curve_map[key], color=color, linewidth=2.0 if key != 'mixed' else 2.2, label=label, zorder=3)
        manual = manual_map.get(key)
        if manual is not None:
            ax.plot(x_il, manual, color=color, linewidth=1.6, linestyle='--', zorder=2)
        x_pts, y_pts = digitized.get(key, (np.asarray([], dtype=float), np.asarray([], dtype=float)))
        if x_pts.size:
            ax.scatter(x_pts, y_pts, s=22, facecolor='white', edgecolor=color, linewidth=0.9, zorder=5)
    if any(arr is not None for arr in (water_manual, il_manual, epc_manual)):
        ax.plot([], [], color='0.25', linestyle='--', linewidth=1.6, label=r'Eq. 8 binary-path FD')
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(float(ylim[0]), float(ylim[1]))
    ax.set_xlabel(r'IL mole fraction, $x_{IL}$')
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(True, alpha=0.3)


def _auto_ylim(*series: np.ndarray) -> tuple[float, float]:
    finite = [np.asarray(s, dtype=float)[np.isfinite(np.asarray(s, dtype=float))] for s in series if s is not None]
    finite = [s for s in finite if s.size]
    if not finite:
        return (-1.0, 1.0)
    vals = np.concatenate(finite)
    y_min = float(np.min(vals))
    y_max = float(np.max(vals))
    if math.isclose(y_min, y_max, rel_tol=0.0, abs_tol=1.0e-12):
        pad = max(0.5, abs(y_min) * 0.1)
        return (y_min - pad, y_max + pad)
    span = y_max - y_min
    pad = 0.07 * span
    return (y_min - pad, y_max + pad)


def main() -> None:
    x_il = np.linspace(1.0e-4, 1.0 - 1.0e-4, 601)
    dadx_water, dadx_water_manual, lnfug_water = _curve_for_mode('orig_water', x_il)
    dadx_il, dadx_il_manual, lnfug_il = _curve_for_mode('orig_il', x_il)
    dadx_epc, dadx_epc_manual, lnfug_epc = _curve_for_mode('epc', x_il)
    digitized_a = _read_digitized_wide(FIG2A_CSV)
    digitized_b = _read_digitized_wide(FIG2B_CSV)

    fig_a, ax_a = plt.subplots(figsize=(6.4, 4.8))
    _plot_single(
        ax_a,
        x_il,
        dadx_water,
        dadx_il,
        dadx_epc,
        digitized_a,
        ylabel=rf'$\partial a^{{ion}} / \partial x_{{{CATION.replace("+", "")}}}$ / -',
        title=rf'Bulow 2019 Figure 2a style: water + [{CATION.replace("+", "")}][{ANION.replace("-", "")}]',
        ylim=(-3.5, 0.0),
        water_manual=dadx_water_manual,
        il_manual=dadx_il_manual,
        epc_manual=dadx_epc_manual,
    )
    ax_a.legend(fontsize=8)
    out_a = THIS_DIR / 'figure_2a.png'
    fig_a.tight_layout()
    fig_a.savefig(out_a, dpi=220)
    plt.close(fig_a)

    fig_a_full, ax_a_full = plt.subplots(figsize=(6.4, 4.8))
    ylim_a_full = _auto_ylim(
        dadx_water,
        dadx_il,
        dadx_epc,
        dadx_water_manual,
        dadx_il_manual,
        dadx_epc_manual,
        *(digitized_a[key][1] for key in digitized_a),
    )
    _plot_single(
        ax_a_full,
        x_il,
        dadx_water,
        dadx_il,
        dadx_epc,
        digitized_a,
        ylabel=rf'$\partial a^{{ion}} / \partial x_{{{CATION.replace("+", "")}}}$ / -',
        title=rf'Bulow 2019 Figure 2a style: water + [{CATION.replace("+", "")}][{ANION.replace("-", "")}] (full range)',
        ylim=ylim_a_full,
        water_manual=dadx_water_manual,
        il_manual=dadx_il_manual,
        epc_manual=dadx_epc_manual,
    )
    ax_a_full.legend(fontsize=8)
    out_a_full = THIS_DIR / 'figure_2a_fullrange.png'
    fig_a_full.tight_layout()
    fig_a_full.savefig(out_a_full, dpi=220)
    plt.close(fig_a_full)

    fig_b, ax_b = plt.subplots(figsize=(6.4, 4.8))
    _plot_single(
        ax_b,
        x_il,
        lnfug_water,
        lnfug_il,
        lnfug_epc,
        digitized_b,
        ylabel=rf'$\ln \varphi_{{{CATION.replace("+", "")}}}$ / -',
        title=rf'Bulow 2019 Figure 2b style: water + [{CATION.replace("+", "")}][{ANION.replace("-", "")}]',
        ylim=(-2.0, 14.0),
    )
    ax_b.legend(fontsize=8)
    out_b = THIS_DIR / 'figure_2b.png'
    fig_b.tight_layout()
    fig_b.savefig(out_b, dpi=220)
    plt.close(fig_b)

    fig_b_full, ax_b_full = plt.subplots(figsize=(6.4, 4.8))
    ylim_b_full = _auto_ylim(
        lnfug_water,
        lnfug_il,
        lnfug_epc,
        *(digitized_b[key][1] for key in digitized_b),
    )
    _plot_single(
        ax_b_full,
        x_il,
        lnfug_water,
        lnfug_il,
        lnfug_epc,
        digitized_b,
        ylabel=rf'$\ln \varphi_{{{CATION.replace("+", "")}}}$ / -',
        title=rf'Bulow 2019 Figure 2b style: water + [{CATION.replace("+", "")}][{ANION.replace("-", "")}] (full range)',
        ylim=ylim_b_full,
    )
    ax_b_full.legend(fontsize=8)
    out_b_full = THIS_DIR / 'figure_2b_fullrange.png'
    fig_b_full.tight_layout()
    fig_b_full.savefig(out_b_full, dpi=220)
    plt.close(fig_b_full)

    fig, axes = plt.subplots(1, 2, figsize=(12.4, 4.8), sharex=False)
    _plot_single(
        axes[0],
        x_il,
        dadx_water,
        dadx_il,
        dadx_epc,
        digitized_a,
        ylabel=rf'$\partial a^{{ion}} / \partial x_{{{CATION.replace("+", "")}}}$ / -',
        title='(a)',
        ylim=(-3.5, 0.0),
        water_manual=dadx_water_manual,
        il_manual=dadx_il_manual,
        epc_manual=dadx_epc_manual,
    )
    _plot_single(
        axes[1],
        x_il,
        lnfug_water,
        lnfug_il,
        lnfug_epc,
        digitized_b,
        ylabel=rf'$\ln \varphi_{{{CATION.replace("+", "")}}}$ / -',
        title='(b)',
        ylim=(-2.0, 14.0),
    )
    axes[1].legend(fontsize=8, loc='best')
    out = THIS_DIR / 'figure_2.png'
    fig.tight_layout()
    fig.savefig(out, dpi=220)
    plt.close(fig)
    print(f'Wrote: {out_a}')
    print(f'Wrote: {out_a_full}')
    print(f'Wrote: {out_b}')
    print(f'Wrote: {out_b_full}')
    print(f'Wrote: {out}')


if __name__ == '__main__':
    main()
