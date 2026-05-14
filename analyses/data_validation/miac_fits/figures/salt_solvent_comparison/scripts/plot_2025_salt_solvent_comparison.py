from __future__ import annotations

import statistics
import sys


from pathlib import Path
import sys as _bootstrap_sys
from pathlib import Path as _BootstrapPath

for _candidate in _BootstrapPath(__file__).resolve().parents:
    if (_candidate / "scripts" / "plot_outputs.py").is_file():
        if str(_candidate) not in _bootstrap_sys.path:
            _bootstrap_sys.path.insert(0, str(_candidate))
        break
else:
    raise ModuleNotFoundError("Could not locate repo root containing scripts/plot_outputs.py")
from scripts.plot_outputs import REPO_ROOT

if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

from analyses.data_validation.miac_fits.scripts import validate_miac_fits as vmf
from scripts.plot_outputs import fits_plot_path, save_plot_figure

OUT_DIR = fits_plot_path("miac", "2025_solvent_comparison", "_placeholder").parent
OUT_DIR.mkdir(parents=True, exist_ok=True)

SOLVENT_ORDER = ["water", "methanol", "ethanol"]
SOLVENT_COLORS = {
    "water": "tab:blue",
    "methanol": "tab:orange",
    "ethanol": "tab:green",
}
SALT_ORDER = ["LiCl", "LiBr", "LiI", "NaCl", "NaBr", "NaI", "KCl", "KBr", "KI"]


def _single_solvent_combo_map():
    combos = vmf.discover_combos()
    out = {}
    for combo in combos:
        solvent = str(combo["solvent_system"])
        if solvent not in SOLVENT_ORDER:
            continue
        if combo.get("comp_signature"):
            continue
        out[(str(combo["salt"]), solvent)] = combo
    return out


def _load_rows_for_combo(combo, quantity: str = "miac_m"):
    molal_exp_raw, values_exp_raw, _, source_exp_raw = vmf.load_exp_data(combo, quantity=quantity)
    keep = vmf._high_outlier_mask(values_exp_raw)
    if np.any(keep):
        molal_exp = molal_exp_raw[keep]
        values_exp = values_exp_raw[keep]
        source_exp = source_exp_raw[keep]
    else:
        molal_exp = molal_exp_raw
        values_exp = values_exp_raw
        source_exp = source_exp_raw
    return molal_exp, values_exp, source_exp


def _choose_low_source(molal_exp: np.ndarray, values_exp: np.ndarray, source_exp: np.ndarray) -> str:
    by_source = {}
    for y, src in zip(values_exp.tolist(), source_exp.tolist()):
        by_source.setdefault(str(src), []).append(float(y))
    ranked = sorted(by_source.items(), key=lambda kv: (statistics.median(kv[1]), statistics.mean(kv[1]), kv[0]))
    return ranked[0][0]


def _curve_for_combo(combo, quantity: str = "miac_m"):
    molal_exp, values_exp, source_exp = _load_rows_for_combo(combo, quantity=quantity)
    if str(combo["salt"]) == "LiBr" and str(combo["solvent_system"]) == "methanol":
        chosen = "Safarov"
        mask = source_exp == chosen
        if np.any(mask):
            molal_exp = molal_exp[mask]
            values_exp = values_exp[mask]
        else:
            chosen = None
    else:
        chosen = None

    if len(molal_exp) == 0:
        raise ValueError(f"No experimental points available for {combo['salt']} in {combo['solvent_system']}.")

    m_upper = vmf._molality_axis_upper(float(np.max(molal_exp)))
    molal_grid = np.linspace(0.0, m_upper, 701)
    curve = vmf.calc_curve(combo, "2025_Figiel", molal_grid, quantity=quantity)
    return molal_exp, values_exp, molal_grid, curve, chosen


def make_plot_for_salt(salt: str) -> Path:
    combo_map = _single_solvent_combo_map()
    fig, ax = plt.subplots(figsize=(8.8, 5.8))
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")

    ymax = 1.0
    plotted_solvents = []
    low_source_note = None

    for solvent in SOLVENT_ORDER:
        combo = combo_map.get((salt, solvent))
        if combo is None:
            continue
        color = SOLVENT_COLORS[solvent]
        molal_exp, values_exp, molal_grid, curve, chosen = _curve_for_combo(combo, quantity="miac_m")
        ax.scatter(
            molal_exp,
            values_exp,
            s=34,
            facecolors="none",
            edgecolors=color,
            linewidths=1.1,
            zorder=4,
        )
        ax.plot(
            molal_grid,
            curve,
            color=color,
            linewidth=2.0,
            linestyle="-",
            zorder=3,
        )
        ymax = max(ymax, float(np.max(values_exp)), float(np.nanmax(curve)))
        plotted_solvents.append(solvent)
        if chosen is not None:
            low_source_note = f"Methanol LiBr data uses Safarov source."

    ax.set_title(f"{salt} MIAC_m fit comparison (2025 model)", color="black", fontsize=13)
    ax.set_xlabel(r"molality, $m$ / mol kg$^{-1}$", fontsize=12)
    ax.set_ylabel(r"mean ionic activity coefficient, $\gamma_{\pm}^{m,*}$", fontsize=12)
    ax.set_xlim(left=0.0)
    ax.set_ylim(0.0, max(1.0, np.ceil(ymax)))
    ax.grid(True, alpha=0.28, color="0.75")
    ax.tick_params(colors="black", labelsize=10)
    for spine in ax.spines.values():
        spine.set_color("black")

    legend_handles = []
    for solvent in plotted_solvents:
        color = SOLVENT_COLORS[solvent]
        legend_handles.append(
            Line2D(
                [0],
                [0],
                color=color,
                marker="o",
                markerfacecolor="white",
                markeredgecolor=color,
                linewidth=2.0,
                label=solvent.capitalize(),
            )
        )
    if legend_handles:
        legend = ax.legend(handles=legend_handles, title="Solvent", fontsize=9, title_fontsize=9, loc="best")
        legend.get_frame().set_facecolor("white")
        legend.get_frame().set_edgecolor("black")
        legend.get_frame().set_alpha(1.0)

    if low_source_note:
        ax.text(0.99, 0.02, low_source_note, transform=ax.transAxes, ha="right", va="bottom", fontsize=8, color="black")

    out = OUT_DIR / f"miac_m_2025_{salt}_solvent_comparison.png"
    fig.tight_layout()
    save_plot_figure(fig, out, dpi=220, bbox_inches=None)
    plt.close(fig)
    return out


def main():
    generated = []
    for salt in SALT_ORDER:
        generated.append(make_plot_for_salt(salt))
    for path in generated:
        print(path)


if __name__ == "__main__":
    main()

