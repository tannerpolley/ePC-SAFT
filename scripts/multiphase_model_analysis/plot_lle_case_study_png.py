from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib
import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from scripts.multiphase_model_analysis import ascani_case2_dataset_comparison as case2

SPECIES = ["H2O", "Butanol", "Na+", "K+", "Cl-"]
SPECIES_LABELS = {
    "H2O": "Water",
    "Butanol": "1-Butanol",
    "Na+": "Na+",
    "K+": "K+",
    "Cl-": "Cl-",
}
SPECIES_COLORS = {
    "H2O": "#2b6cb0",
    "Butanol": "#c05621",
    "Na+": "#2f855a",
    "K+": "#805ad5",
    "Cl-": "#4a5568",
}
PHASE_COLORS = {
    "organic": "#8b1e3f",
    "aqueous": "#2b6cb0",
}
DEFAULT_MODEL_KEY = "ascani2022_params_figiel2025_opts"


def _select_config(model_key: str) -> dict:
    configs = case2._default_model_configs()
    for config in configs:
        if config["key"] == model_key:
            return config
    available = ", ".join(config["key"] for config in configs)
    raise KeyError(f"Unknown model_key '{model_key}'. Available: {available}")


def _fmt_mole_fraction(value: float) -> str:
    value = float(value)
    if value >= 1.0e-2:
        return f"{value:.4f}"
    if value >= 1.0e-4:
        return f"{value:.5f}"
    return f"{value:.2e}"


def build_lle_case_study_png(out_path: Path, model_key: str = DEFAULT_MODEL_KEY) -> Path:
    config = _select_config(model_key)
    result = case2._solve_dataset(config)
    diag = result["paper_compare"]

    feed_vals = np.array([result["feed_z"][sp] for sp in SPECIES], dtype=float)
    org_phase = result["phases"]["organic"]
    aq_phase = result["phases"]["aqueous"]
    org_vals = np.array([org_phase["x"][sp] for sp in SPECIES], dtype=float)
    aq_vals = np.array([aq_phase["x"][sp] for sp in SPECIES], dtype=float)
    org_share = np.array([org_phase["share_of_feed_pct"][sp] for sp in SPECIES], dtype=float)
    aq_share = np.array([aq_phase["share_of_feed_pct"][sp] for sp in SPECIES], dtype=float)

    fig = plt.figure(figsize=(13.6, 7.6), facecolor="white")
    fig.text(0.5, 0.965, "LLE case study: feed splitting into two liquid phases", ha="center", va="top", fontsize=19, fontweight="bold")
    fig.text(0.5, 0.928, "Water + 1-butanol + NaCl + KCl", ha="center", va="top", fontsize=13)

    ax_table = fig.add_axes([0.05, 0.16, 0.43, 0.68])
    ax_table.axis("off")

    row_labels = [SPECIES_LABELS[sp] for sp in SPECIES]
    cell_text = [
        [_fmt_mole_fraction(feed_vals[i]), _fmt_mole_fraction(org_vals[i]), _fmt_mole_fraction(aq_vals[i])]
        for i in range(len(SPECIES))
    ]
    table = ax_table.table(
        cellText=cell_text,
        rowLabels=row_labels,
        colLabels=[r"$z_i^{feed}$", r"$x_i^{org}$", r"$x_i^{aq}$"],
        loc="center",
        cellLoc="center",
        rowLoc="center",
        bbox=[0.05, 0.15, 0.90, 0.70],
    )
    table.auto_set_font_size(False)
    table.set_fontsize(11)
    table.scale(1.2, 1.7)
    for (r, c), cell in table.get_celld().items():
        cell.set_edgecolor("#cbd5e0")
        cell.set_linewidth(1.0)
        if r == 0:
            cell.set_facecolor("#edf2f7")
            cell.get_text().set_fontweight("bold")
        elif c == -1 and r > 0:
            sp = SPECIES[r - 1]
            cell.set_facecolor("white")
            cell.get_text().set_color(SPECIES_COLORS[sp])
            cell.get_text().set_fontweight("bold")
        else:
            cell.set_facecolor("white")

    ax_table.text(0.5, 0.93, "Compositions before and after phase split", ha="center", va="center", fontsize=14, fontweight="bold")
    ax_table.text(0.5, 0.07, "Mole fractions are normalized within the feed or within each equilibrium phase.", ha="center", va="center", fontsize=10)

    ax_beta = fig.add_axes([0.56, 0.72, 0.36, 0.10])
    ax_beta.set_title("Phase fractions at equilibrium", fontsize=13, pad=8)
    ax_beta.barh([0], [100.0 * org_phase["beta"]], color=PHASE_COLORS["organic"], height=0.42)
    ax_beta.barh([0], [100.0 * aq_phase["beta"]], left=[100.0 * org_phase["beta"]], color=PHASE_COLORS["aqueous"], height=0.42)
    ax_beta.text(100.0 * org_phase["beta"] * 0.5, 0.34, rf"$\beta_{{org}}$ = {100.0 * org_phase['beta']:.2f}\%", ha="center", va="bottom", fontsize=10, color=PHASE_COLORS["organic"], fontweight="bold")
    ax_beta.text(100.0 * org_phase["beta"] + 100.0 * aq_phase["beta"] * 0.5, 0.34, rf"$\beta_{{aq}}$ = {100.0 * aq_phase['beta']:.2f}\%", ha="center", va="bottom", fontsize=10, color=PHASE_COLORS["aqueous"], fontweight="bold")
    ax_beta.set_xlim(0.0, 100.0)
    ax_beta.set_xticks([])
    ax_beta.set_yticks([])
    for spine in ax_beta.spines.values():
        spine.set_visible(False)

    ax_thermo = fig.add_axes([0.56, 0.50, 0.36, 0.16])
    ax_thermo.axis("off")
    thermo_lines = [
        r"Thermodynamic driving force",
        rf"$\hat g_{{feed}} = {diag['ghat_feed_j_per_mol']:.1f}\;\mathrm{{J\,mol^{{-1}}}}$",
        rf"$\hat g_{{eq}} = {diag['ghat_eq_j_per_mol']:.1f}\;\mathrm{{J\,mol^{{-1}}}}$",
        rf"$\Delta \hat g = {diag['ghat_delta_j_per_mol']:.2f}\;\mathrm{{J\,mol^{{-1}}}}$",
        rf"$TPDF_{{min}} = {diag['tpdf_min']:.4f}$",
    ]
    ax_thermo.text(
        0.5,
        0.50,
        "\n".join(thermo_lines),
        ha="center",
        va="center",
        fontsize=11,
        bbox={"boxstyle": "round,pad=0.55", "facecolor": "#f7fafc", "edgecolor": "#cbd5e0", "linewidth": 1.0},
    )

    ax_part = fig.add_axes([0.56, 0.14, 0.36, 0.26])
    y_positions = np.arange(len(SPECIES), dtype=float)
    for y, sp in zip(y_positions, SPECIES):
        idx = SPECIES.index(sp)
        x_org = org_share[idx]
        x_aq = aq_share[idx]
        ax_part.hlines(y, min(x_org, x_aq), max(x_org, x_aq), color=SPECIES_COLORS[sp], linewidth=5.0, alpha=0.25)
        ax_part.scatter(x_org, y, s=95, color=PHASE_COLORS["organic"], edgecolor="black", linewidth=0.7, zorder=3)
        ax_part.scatter(x_aq, y, s=95, color=PHASE_COLORS["aqueous"], edgecolor="black", linewidth=0.7, zorder=3)
        ax_part.text(x_org - 1.8, y + 0.17, f"{x_org:.1f}", ha="right", va="bottom", fontsize=9, color=PHASE_COLORS["organic"])
        ax_part.text(x_aq + 1.8, y + 0.17, f"{x_aq:.1f}", ha="left", va="bottom", fontsize=9, color=PHASE_COLORS["aqueous"])

    ax_part.set_xlim(0.0, 100.0)
    ax_part.set_ylim(-0.6, len(SPECIES) - 0.4)
    ax_part.set_yticks(y_positions)
    ax_part.set_yticklabels([SPECIES_LABELS[sp] for sp in SPECIES])
    ax_part.invert_yaxis()
    ax_part.set_xlabel(r"Fraction of each original component reporting to phase (\%)")
    ax_part.set_title("Component partitioning from feed to the two equilibrium phases", fontsize=13, pad=8)
    ax_part.grid(axis="x", color="#d7d7d7", linewidth=0.7, alpha=0.8)
    ax_part.set_axisbelow(True)
    ax_part.spines["top"].set_visible(False)
    ax_part.spines["right"].set_visible(False)
    ax_part.legend(
        handles=[
            Line2D([0], [0], marker="o", color="none", markerfacecolor=PHASE_COLORS["organic"], markeredgecolor="black", markersize=8, label="Organic-rich phase"),
            Line2D([0], [0], marker="o", color="none", markerfacecolor=PHASE_COLORS["aqueous"], markeredgecolor="black", markersize=8, label="Aqueous-rich phase"),
        ],
        loc="lower right",
        frameon=False,
    )

    fig.text(0.74, 0.09, r"Negative $TPDF_{min}$ and lower $\hat g$ indicate that phase splitting is favored.", ha="center", va="center", fontsize=10)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=240, bbox_inches="tight")
    plt.close(fig)
    return out_path


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Create a slide-ready LLE case-study PNG.")
    parser.add_argument(
        "--out",
        type=Path,
        default=REPO_ROOT / "scripts" / "multiphase_model_analysis" / "output" / "ascani_case2_figiel2025_lle_case_study.png",
        help="Output PNG path.",
    )
    parser.add_argument(
        "--model-key",
        default=DEFAULT_MODEL_KEY,
        help="Model key from ascani_case2_dataset_comparison._default_model_configs().",
    )
    return parser.parse_args()


def main() -> None:
    args = _parse_args()
    out = build_lle_case_study_png(args.out, model_key=args.model_key)
    print(f"Saved LLE case-study PNG: {out}")


if __name__ == "__main__":
    main()

