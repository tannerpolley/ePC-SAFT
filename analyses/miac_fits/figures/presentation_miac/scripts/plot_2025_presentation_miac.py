from __future__ import annotations

import math
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[5]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

from analyses.miac_fits.scripts import validate_miac_fits as vmf
from scripts.plot_outputs import fits_plot_path, save_plot_figure

OUT_DIR = fits_plot_path("miac", "2025_presentation", "_placeholder").parent
OUT_DIR.mkdir(parents=True, exist_ok=True)

SOLVENT_COLORS = {
    "water": "tab:blue",
    "methanol": "tab:orange",
    "ethanol": "tab:green",
}

SALT_COLORS = {
    "LiCl": "tab:blue",
    "LiBr": "tab:orange",
    "LiI": "tab:green",
    "NaCl": "tab:red",
    "NaBr": "tab:purple",
    "NaI": "tab:brown",
    "KCl": "tab:olive",
    "KBr": "tab:cyan",
    "KI": "tab:pink",
}

SPECIAL_SOURCE = {
    ("LiBr", "methanol"): "Safarov",
}

PRESENTATION_GROUPS = [
    {
        "slug": "nacl_solvents",
        "title": "NaCl in water, methanol, and ethanol",
        "style_by": "solvent",
        "series": [
            {"salt": "NaCl", "solvent": "water", "label": "Water"},
            {"salt": "NaCl", "solvent": "methanol", "label": "Methanol"},
            {"salt": "NaCl", "solvent": "ethanol", "label": "Ethanol"},
        ],
    },
    {
        "slug": "nabr_solvents",
        "title": "NaBr in water, methanol, and ethanol",
        "style_by": "solvent",
        "series": [
            {"salt": "NaBr", "solvent": "water", "label": "Water"},
            {"salt": "NaBr", "solvent": "methanol", "label": "Methanol"},
            {"salt": "NaBr", "solvent": "ethanol", "label": "Ethanol"},
        ],
    },
    {
        "slug": "lii_solvents",
        "title": "LiI in water, methanol, and ethanol",
        "style_by": "solvent",
        "series": [
            {"salt": "LiI", "solvent": "water", "label": "Water"},
            {"salt": "LiI", "solvent": "methanol", "label": "Methanol"},
            {"salt": "LiI", "solvent": "ethanol", "label": "Ethanol"},
        ],
    },
    {
        "slug": "libr_nonaqueous",
        "title": "LiBr in methanol and ethanol",
        "style_by": "solvent",
        "series": [
            {"salt": "LiBr", "solvent": "methanol", "label": "Methanol"},
            {"salt": "LiBr", "solvent": "ethanol", "label": "Ethanol"},
        ],
    },
    {
        "slug": "methanol_chlorides",
        "title": "Chlorides in methanol",
        "style_by": "salt",
        "series": [
            {"salt": "LiCl", "solvent": "methanol", "label": "LiCl"},
            {"salt": "NaCl", "solvent": "methanol", "label": "NaCl"},
            {"salt": "KCl", "solvent": "methanol", "label": "KCl"},
        ],
    },
    {
        "slug": "methanol_iodides",
        "title": "Iodides in methanol",
        "style_by": "salt",
        "series": [
            {"salt": "LiI", "solvent": "methanol", "label": "LiI"},
            {"salt": "NaI", "solvent": "methanol", "label": "NaI"},
            {"salt": "KI", "solvent": "methanol", "label": "KI"},
        ],
    },
    {
        "slug": "water_potassium_halides",
        "title": "Potassium halides in water",
        "style_by": "salt",
        "series": [
            {"salt": "KCl", "solvent": "water", "label": "KCl"},
            {"salt": "KBr", "solvent": "water", "label": "KBr"},
            {"salt": "KI", "solvent": "water", "label": "KI"},
        ],
    },
    {
        "slug": "ethanol_sodium_salts",
        "title": "Sodium salts in ethanol",
        "style_by": "salt",
        "series": [
            {"salt": "NaCl", "solvent": "ethanol", "label": "NaCl"},
            {"salt": "NaBr", "solvent": "ethanol", "label": "NaBr"},
            {"salt": "NaI", "solvent": "ethanol", "label": "NaI"},
        ],
    },
]


def _combo_map():
    out = {}
    for combo in vmf.discover_combos():
        if combo.get("comp_signature"):
            continue
        solvent = str(combo["solvent_system"])
        salt = str(combo["salt"])
        out[(salt, solvent)] = combo
    return out


def _nice_upper(value: float) -> float:
    if value <= 0.0:
        return 1.0
    padded = value * 1.05
    if padded <= 1.0:
        return math.ceil(padded * 10.0) / 10.0
    if padded <= 2.0:
        return math.ceil(padded * 5.0) / 5.0
    if padded <= 10.0:
        return float(math.ceil(padded))
    if padded <= 25.0:
        return float(math.ceil(padded / 2.0) * 2.0)
    return float(math.ceil(padded / 5.0) * 5.0)


def _series_color(style_by: str, salt: str, solvent: str) -> str:
    if style_by == "solvent":
        return SOLVENT_COLORS[solvent]
    return SALT_COLORS[salt]


def _load_series(combo: dict[str, object]):
    salt = str(combo["salt"])
    solvent = str(combo["solvent_system"])
    molal_exp_raw, values_exp_raw, _, source_exp_raw = vmf.load_exp_data(combo, quantity="miac_m")
    keep = vmf._high_outlier_mask(values_exp_raw)
    if np.any(keep):
        molal_exp = molal_exp_raw[keep]
        values_exp = values_exp_raw[keep]
        source_exp = source_exp_raw[keep]
    else:
        molal_exp = molal_exp_raw
        values_exp = values_exp_raw
        source_exp = source_exp_raw

    chosen_note = None
    wanted_source = SPECIAL_SOURCE.get((salt, solvent))
    if wanted_source is not None:
        mask = source_exp == wanted_source
        if np.any(mask):
            molal_exp = molal_exp[mask]
            values_exp = values_exp[mask]
            source_exp = source_exp[mask]
            chosen_note = f"{solvent.capitalize()} {salt}: {wanted_source} only"

    order = np.argsort(molal_exp)
    molal_exp = molal_exp[order]
    values_exp = values_exp[order]
    source_exp = source_exp[order]
    x_end = float(molal_exp[-1]) if len(molal_exp) else 0.0
    molal_grid = np.linspace(0.0, x_end, 401) if x_end > 0.0 else np.array([0.0])
    curve = vmf.calc_curve(combo, "2025_Figiel", molal_grid, quantity="miac_m")
    return molal_exp, values_exp, molal_grid, curve, chosen_note


def _legend_handle(color: str, label: str):
    return Line2D(
        [0],
        [0],
        color=color,
        linewidth=2.1,
        marker="o",
        markerfacecolor="white",
        markeredgecolor=color,
        markersize=5.5,
        label=label,
    )


def make_plot(group: dict[str, object], combo_map: dict[tuple[str, str], dict[str, object]]) -> Path:
    fig, ax = plt.subplots(figsize=(8.4, 5.0))
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")

    max_x = 0.0
    max_y = 0.0
    legend_handles = []
    notes = []

    for spec in group["series"]:
        salt = str(spec["salt"])
        solvent = str(spec["solvent"])
        combo = combo_map[(salt, solvent)]
        color = _series_color(str(group["style_by"]), salt, solvent)
        label = str(spec["label"])
        molal_exp, values_exp, molal_grid, curve, chosen_note = _load_series(combo)
        ax.scatter(molal_exp, values_exp, s=32, facecolors="white", edgecolors=color, linewidths=1.15, zorder=4)
        ax.plot(molal_grid, curve, color=color, linewidth=2.1, zorder=3)
        max_x = max(max_x, float(molal_exp[-1]) if len(molal_exp) else 0.0)
        max_y = max(
            max_y, float(np.max(values_exp)) if len(values_exp) else 0.0, float(np.max(curve)) if len(curve) else 0.0
        )
        legend_handles.append(_legend_handle(color, label))
        if chosen_note:
            notes.append(chosen_note)

    ax.set_title(str(group["title"]), fontsize=13, color="black")
    ax.set_xlabel(r"molality, $m$ / mol kg$^{-1}$", fontsize=12)
    ax.set_ylabel(r"mean ionic activity coefficient, $\gamma_{\pm}^{m,*}$", fontsize=12)
    ax.set_xlim(0.0, _nice_upper(max_x))
    ax.set_ylim(0.0, _nice_upper(max_y))
    ax.grid(True, alpha=0.25, color="0.75")
    ax.tick_params(axis="both", labelsize=10, colors="black")
    for spine in ax.spines.values():
        spine.set_color("black")

    legend = ax.legend(handles=legend_handles, loc="best", fontsize=9, frameon=True)
    legend.get_frame().set_facecolor("white")
    legend.get_frame().set_edgecolor("black")
    legend.get_frame().set_alpha(1.0)

    ax.text(
        0.01,
        0.99,
        "Open circles: data\nSolid lines: 2025 fit",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=8.5,
        color="black",
    )
    if notes:
        ax.text(
            0.99, 0.01, "; ".join(notes), transform=ax.transAxes, ha="right", va="bottom", fontsize=8, color="black"
        )

    out = OUT_DIR / f"miac_m_2025_present_{group['slug']}.png"
    fig.tight_layout()
    save_plot_figure(fig, out, dpi=240, bbox_inches=None)
    plt.close(fig)
    return out


def write_summary(paths: list[Path]) -> Path:
    summary = OUT_DIR / "miac_m_2025_presentation_plot_set.md"
    lines = [
        "# 2025 MIAC presentation plot set",
        "",
        "These plots were curated for slide use rather than exhaustiveness.",
        "",
        "Selection rules:",
        "- group only series with comparable molality and MIAC_m ranges",
        "- stop each fit at the last retained experimental point for that series",
        "- use Safarov only for methanol LiBr",
        "- avoid combinations dominated by the very large aqueous LiBr range",
        "",
        "Generated plots:",
    ]
    for path in paths:
        lines.append(f"- {path.name}")
    summary.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return summary


def main() -> None:
    combo_map = _combo_map()
    generated = [make_plot(group, combo_map) for group in PRESENTATION_GROUPS]
    summary = write_summary(generated)
    for path in generated:
        print(path)
    print(summary)


if __name__ == "__main__":
    main()

