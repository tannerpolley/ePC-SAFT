from __future__ import annotations

import csv
import importlib.util
import math
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use('Agg')
import matplotlib.pyplot as plt

FIG2_DIR = Path(__file__).resolve().parents[1]
PLOT_SCRIPT = FIG2_DIR / 'plot_figure_2.py'
spec = importlib.util.spec_from_file_location('fig2mod', PLOT_SCRIPT)
fig2 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(fig2)

THIS_DIR = Path(__file__).resolve().parent
OUT_CSV = THIS_DIR / 'figure2_quantity_fit_summary.csv'
OUT_MD = THIS_DIR / 'figure2_quantity_fit_notes.md'
OUT_PLOT_A = THIS_DIR / 'figure2a_candidate_overlay.png'
OUT_PLOT_B = THIS_DIR / 'figure2b_candidate_overlay.png'
R = 8.31446261815324
SERIES = [('orig_water', 'water', 'tab:blue'), ('epc', 'mixed', 'green'), ('orig_il', 'il', 'tab:orange')]


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


def _interp_metrics(x_model: np.ndarray, y_model: np.ndarray, x_data: np.ndarray, y_data: np.ndarray) -> dict[str, float]:
    mask = np.isfinite(y_model)
    y_interp = np.interp(x_data, x_model[mask], y_model[mask])
    diff = y_interp - y_data
    shift = float(np.mean(y_data - y_interp))
    shifted = y_interp + shift
    return {
        'rmse': float(np.sqrt(np.mean(diff * diff))),
        'mae': float(np.mean(np.abs(diff))),
        'bias': float(np.mean(diff)),
        'shift': shift,
        'shifted_rmse': float(np.sqrt(np.mean((shifted - y_data) ** 2))),
        'y_start': float(y_model[mask][0]),
        'y_end': float(y_model[mask][-1]),
    }


def _build_params(mode: str, *, numerical: bool = False) -> dict:
    params = fig2._figure2_params(mode)
    if numerical:
        params = dict(params)
        elec = dict(params['elec_model'])
        relp = dict(elec['rel_perm'])
        relp['differential_mode'] = 'numerical'
        elec['rel_perm'] = relp
        dh = dict(elec['DH_model'])
        mudh = dict(dh['mu_DH_model'])
        mudh['differential_mode'] = 'numerical'
        dh['mu_DH_model'] = mudh
        elec['DH_model'] = dh
        params['elec_model'] = elec
    return params


def _curve_bundle(mode: str, *, numerical: bool = False) -> tuple[np.ndarray, dict[str, np.ndarray]]:
    x = np.linspace(1.0e-4, 1.0 - 1.0e-4, 601)
    params = _build_params(mode, numerical=numerical)
    out = {name: np.full_like(x, np.nan, dtype=float) for name in [
        'dadx_cat_api', 'dadx_an_api', 'dadx_mean_salt', 'fd_aion_dxil_resolvedrho',
        'lnphi_total', 'lnphi_ion', 'lnphi_ion_stdcat', 'lnphi_ion_stdil', 'lnphi_mean_salt_stdil',
    ]}
    for idx, xil in enumerate(x):
        comp = fig2._composition_from_x_il(float(xil))
        rho = float(fig2.epcsaft_density(fig2.T_REF, fig2.P_REF, comp, params, phase='liq'))
        terms = fig2.epcsaft_fugacity_coefficient_terms(fig2.T_REF, rho, comp, params)
        dadx = np.asarray(terms['dadx_ion'], dtype=float)
        lnphi = np.asarray(terms['lnfugcoef'], dtype=float)
        lnphi_ion = np.asarray(terms['lnfugcoef_ion'], dtype=float)
        out['dadx_cat_api'][idx] = float(dadx[fig2.CATION_INDEX])
        out['dadx_an_api'][idx] = float(dadx[fig2.ANION_INDEX])
        # Mean-ionic salt-basis hypothesis: 0.5 stoichiometric average times dx_i/dx_IL = 0.5.
        out['dadx_mean_salt'][idx] = 0.25 * float(dadx[fig2.CATION_INDEX] + dadx[fig2.ANION_INDEX])
        h = min(1.0e-4, 0.25 * min(float(xil), 1.0 - float(xil)))
        if h > 0.0:
            vals = []
            for xh in (xil - h, xil + h):
                comph = fig2._composition_from_x_il(float(xh))
                rhoh = float(fig2.epcsaft_density(fig2.T_REF, fig2.P_REF, comph, params, phase='liq'))
                termsh = fig2.epcsaft_fugacity_coefficient_terms(fig2.T_REF, rhoh, comph, params)
                vals.append(float(termsh['a_ion']))
            out['fd_aion_dxil_resolvedrho'][idx] = (vals[1] - vals[0]) / (2.0 * h)
        lnxcat = math.log(float(comp[fig2.CATION_INDEX]))
        lnxil = math.log(float(xil))
        std = math.log(rho * R * fig2.T_REF / fig2.P_REF)
        out['lnphi_total'][idx] = float(lnphi[fig2.CATION_INDEX])
        out['lnphi_ion'][idx] = float(lnphi_ion[fig2.CATION_INDEX])
        out['lnphi_ion_stdcat'][idx] = float(lnphi_ion[fig2.CATION_INDEX]) - lnxcat + std
        out['lnphi_ion_stdil'][idx] = float(lnphi_ion[fig2.CATION_INDEX]) - lnxil + std
        out['lnphi_mean_salt_stdil'][idx] = 0.5 * (
            float(lnphi_ion[fig2.CATION_INDEX]) - lnxil + std +
            float(lnphi_ion[fig2.ANION_INDEX]) - lnxil + std
        )
    return x, out


def main() -> None:
    dig_a = _read_digitized_wide(fig2.FIG2A_CSV)
    dig_b = _read_digitized_wide(fig2.FIG2B_CSV)
    rows: list[dict[str, object]] = []
    stored: dict[tuple[str, str], tuple[np.ndarray, dict[str, np.ndarray]]] = {}

    for mode, key, _ in SERIES:
        x, out = _curve_bundle(mode, numerical=False)
        stored[(mode, 'analytical')] = (x, out)
        x_num, out_num = _curve_bundle(mode, numerical=True)
        stored[(mode, 'numerical')] = (x_num, out_num)

        xda, yda = dig_a[key]
        xdb, ydb = dig_b[key]

        for quantity in ['dadx_cat_api', 'dadx_mean_salt', 'fd_aion_dxil_resolvedrho']:
            metrics = _interp_metrics(x, out[quantity], xda, yda)
            rows.append({'panel': '2a', 'mode': mode, 'series': key, 'quantity': quantity, 'variant': 'analytical', **metrics})
        metrics_num = _interp_metrics(x_num, out_num['dadx_cat_api'], xda, yda)
        rows.append({'panel': '2a', 'mode': mode, 'series': key, 'quantity': 'dadx_cat_api', 'variant': 'numerical', **metrics_num})
        metrics_num_mean = _interp_metrics(x_num, out_num['dadx_mean_salt'], xda, yda)
        rows.append({'panel': '2a', 'mode': mode, 'series': key, 'quantity': 'dadx_mean_salt', 'variant': 'numerical', **metrics_num_mean})

        for quantity in ['lnphi_total', 'lnphi_ion', 'lnphi_ion_stdcat', 'lnphi_ion_stdil', 'lnphi_mean_salt_stdil']:
            metrics = _interp_metrics(x, out[quantity], xdb, ydb)
            rows.append({'panel': '2b', 'mode': mode, 'series': key, 'quantity': quantity, 'variant': 'analytical', **metrics})

    rows = sorted(rows, key=lambda r: (r['panel'], r['mode'], r['shifted_rmse'], r['rmse'], r['quantity'], r['variant']))
    fieldnames = ['panel', 'mode', 'series', 'quantity', 'variant', 'rmse', 'mae', 'bias', 'shift', 'shifted_rmse', 'y_start', 'y_end']
    with OUT_CSV.open('w', newline='', encoding='utf-8') as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    best_2a: dict[str, dict[str, object]] = {}
    best_2b: dict[str, dict[str, object]] = {}
    for row in rows:
        if row['panel'] == '2a' and row['mode'] not in best_2a:
            best_2a[row['mode']] = row
        if row['panel'] == '2b' and row['mode'] not in best_2b:
            best_2b[row['mode']] = row

    lines = [
        '# Figure 2 Quantity Audit',
        '',
        '## Key Findings',
        '',
        '- The original Figure 2 replication bug was that it zeroed all $k_{ij}$ values, but the 2019 paper uses the Table 6 water-IL unlike parameters.',
        '- For panel 2a, the raw API cation derivative is not the closest salt-basis match. The strongest current candidate is the mean-ionic salt-basis quantity $0.25\\,(d a^{ion}/d x_{cat} + d a^{ion}/d x_{anion})$.',
        '- For panel 2b, the paper scale around $14$ at low $x_{IL}$ is incompatible with the raw API $\\ln\\varphi_{cat}$. The closest current candidate is the concentration-referenced ionic quantity $\\ln\\varphi^{ion}_{cat} - \\ln x_{IL} + \\ln(\\rho RT/P)$.',
        '- Numerical vs analytical DH derivatives are essentially identical for these curves; the main mismatch is the plotted quantity definition, not numerical derivative error.',
        '',
        '## Best Current Matches',
        '',
    ]
    for row in best_2a.values():
        lines.append(f"- Panel 2a `{row['mode']}`: best current quantity is `{row['quantity']}` (`{row['variant']}`), RMSE `{row['rmse']:.4f}`, shifted RMSE `{row['shifted_rmse']:.4f}`.")
    for row in best_2b.values():
        lines.append(f"- Panel 2b `{row['mode']}`: best current quantity is `{row['quantity']}` (`{row['variant']}`), RMSE `{row['rmse']:.4f}`, shifted RMSE `{row['shifted_rmse']:.4f}`.")
    OUT_MD.write_text('\n'.join(lines) + '\n', encoding='utf-8')

    fig_a, axes_a = plt.subplots(1, 3, figsize=(14.5, 4.2), sharey=False)
    for ax, (mode, key, color) in zip(axes_a, SERIES):
        x, out = stored[(mode, 'analytical')]
        x_num, out_num = stored[(mode, 'numerical')]
        xd, yd = dig_a[key]
        ax.scatter(xd, yd, s=22, facecolor='white', edgecolor=color, linewidth=0.9, zorder=6, label='paper')
        ax.plot(x, out['dadx_cat_api'], color=color, linewidth=1.6, alpha=0.45, label='api cation')
        ax.plot(x, out['dadx_mean_salt'], color='black', linewidth=1.7, linestyle=':', zorder=5, label='mean ionic salt-basis analytical')
        ax.plot(x_num, out_num['dadx_mean_salt'], color='0.35', linewidth=1.5, linestyle='--', zorder=4, label='mean ionic salt-basis numerical')
        ax.set_title(mode)
        ax.set_xlim(0.0, 1.0)
        ax.grid(True, alpha=0.25)
    axes_a[0].set_ylabel(r'candidate value / -')
    axes_a[-1].legend(fontsize=7, loc='best')
    fig_a.tight_layout()
    fig_a.savefig(OUT_PLOT_A, dpi=220)
    plt.close(fig_a)

    fig_b, axes_b = plt.subplots(1, 3, figsize=(14.5, 4.2), sharey=False)
    for ax, (mode, key, color) in zip(axes_b, SERIES):
        x, out = stored[(mode, 'analytical')]
        xd, yd = dig_b[key]
        ax.scatter(xd, yd, s=22, facecolor='white', edgecolor=color, linewidth=0.9, zorder=6, label='paper')
        ax.plot(x, out['lnphi_total'], color=color, linewidth=1.6, alpha=0.55, label='raw lnphi total')
        ax.plot(x, out['lnphi_ion_stdil'], color='black', linewidth=1.8, linestyle=':', zorder=5, label=r'lnphi_ion - ln x_IL + ln(rhoRT/P)')
        ax.plot(x, out['lnphi_mean_salt_stdil'], color='0.35', linewidth=1.4, linestyle='--', zorder=4, label='mean ionic salt-basis')
        ax.set_title(mode)
        ax.set_xlim(0.0, 1.0)
        ax.grid(True, alpha=0.25)
    axes_b[0].set_ylabel(r'candidate value / -')
    axes_b[-1].legend(fontsize=7, loc='best')
    fig_b.tight_layout()
    fig_b.savefig(OUT_PLOT_B, dpi=220)
    plt.close(fig_b)

    print(f'Wrote: {OUT_CSV}')
    print(f'Wrote: {OUT_MD}')
    print(f'Wrote: {OUT_PLOT_A}')
    print(f'Wrote: {OUT_PLOT_B}')


if __name__ == '__main__':
    main()
