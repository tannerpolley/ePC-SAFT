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


def _fmt_amount_mmol(value_mol: float) -> str:
    value_mmol = 1000.0 * float(value_mol)
    if value_mmol >= 100.0:
        return f"{value_mmol:.1f}"
    if value_mmol >= 10.0:
        return f"{value_mmol:.2f}"
    if value_mmol >= 0.1:
        return f"{value_mmol:.3f}"
    return f"{value_mmol:.4f}"


def _solve_case(model_key: str) -> dict:
    config = _select_config(model_key)
    return case2._solve_dataset(config)


def _save_composition_table(result: dict, out_path: Path) -> None:
    beta_org = float(result["phases"]["organic"]["beta"])
    beta_aq = float(result["phases"]["aqueous"]["beta"])
    feed_vals = np.array([result["feed_z"][sp] for sp in SPECIES], dtype=float)
    org_vals = np.array([beta_org * result["phases"]["organic"]["x"][sp] for sp in SPECIES], dtype=float)
    aq_vals = np.array([beta_aq * result["phases"]["aqueous"]["x"][sp] for sp in SPECIES], dtype=float)

    fig = plt.figure(figsize=(8.0, 4.9), facecolor="white")
    ax = fig.add_axes([0.04, 0.06, 0.92, 0.88])
    ax.axis("off")

    table = ax.table(
        cellText=[
            [_fmt_amount_mmol(feed_vals[i]), _fmt_amount_mmol(org_vals[i]), _fmt_amount_mmol(aq_vals[i])]
            for i in range(len(SPECIES))
        ],
        rowLabels=[SPECIES_LABELS[sp] for sp in SPECIES],
        colLabels=["Feed", "Organic-rich phase", "Aqueous-rich phase"],
        cellLoc="center",
        rowLoc="center",
        loc="center",
        bbox=[0.10, 0.12, 0.82, 0.68],
    )
    table.auto_set_font_size(False)
    table.set_fontsize(11)
    table.scale(1.15, 1.7)
    for (r, c), cell in table.get_celld().items():
        cell.set_edgecolor("#cbd5e0")
        cell.set_linewidth(1.0)
        if r == 0:
            cell.set_facecolor("#edf2f7")
            cell.get_text().set_fontweight("bold")
        elif c == -1 and r > 0:
            sp = SPECIES[r - 1]
            cell.get_text().set_color(SPECIES_COLORS[sp])
            cell.get_text().set_fontweight("bold")
            cell.set_facecolor("white")
        else:
            cell.set_facecolor("white")

    ax.text(0.5, 0.94, "Component amounts by phase", ha="center", va="center", fontsize=16, fontweight="bold")
    ax.text(0.5, 0.885, "mmol per mol feed", ha="center", va="center", fontsize=11)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=240, bbox_inches="tight")
    plt.close(fig)


def _save_phase_split(result: dict, out_path: Path) -> None:
    beta_org = 100.0 * float(result["phases"]["organic"]["beta"])
    beta_aq = 100.0 * float(result["phases"]["aqueous"]["beta"])

    fig = plt.figure(figsize=(7.6, 2.6), facecolor="white")
    ax = fig.add_axes([0.08, 0.28, 0.84, 0.42])
    ax.set_title("Phase fractions at equilibrium", fontsize=15, pad=10)
    ax.barh([0], [beta_org], color=PHASE_COLORS["organic"], height=0.38)
    ax.barh([0], [beta_aq], left=[beta_org], color=PHASE_COLORS["aqueous"], height=0.38)
    ax.text(beta_org * 0.5, 0.34, rf"$\beta_{{org}}$ = {beta_org:.2f}\%", ha="center", va="bottom", fontsize=11, color=PHASE_COLORS["organic"], fontweight="bold")
    ax.text(beta_org + beta_aq * 0.5, 0.34, rf"$\beta_{{aq}}$ = {beta_aq:.2f}\%", ha="center", va="bottom", fontsize=11, color=PHASE_COLORS["aqueous"], fontweight="bold")
    ax.set_xlim(0.0, 100.0)
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=240, bbox_inches="tight")
    plt.close(fig)


def _save_driving_force(result: dict, out_path: Path) -> None:
    diag = result["paper_compare"]
    g_feed = float(diag["ghat_feed_j_per_mol"])
    g_eq = float(diag["ghat_eq_j_per_mol"])
    delta_g = float(diag["ghat_delta_j_per_mol"])
    tpdf = float(diag["tpdf_min"])

    fig = plt.figure(figsize=(7.6, 3.2), facecolor="white")
    ax = fig.add_axes([0.10, 0.24, 0.82, 0.56])
    ax.set_title(r"Thermodynamic driving force for splitting", fontsize=15, pad=10)

    x_min = min(g_feed, g_eq) - 8.0
    x_max = max(g_feed, g_eq) + 8.0
    ax.hlines(0.0, g_eq, g_feed, color="#718096", linewidth=2.5)
    ax.scatter([g_feed], [0.0], s=120, color="#4a5568", edgecolor="black", zorder=3)
    ax.scatter([g_eq], [0.0], s=120, color="#2b6cb0", edgecolor="black", zorder=3)
    ax.text(g_feed, 0.12, rf"$\hat g_{{feed}}$ = {g_feed:.1f}", ha="center", va="bottom", fontsize=11)
    ax.text(g_eq, -0.18, rf"$\hat g_{{eq}}$ = {g_eq:.1f}", ha="center", va="top", fontsize=11)
    ax.text((g_feed + g_eq) * 0.5, 0.23, rf"$\Delta \hat g$ = {delta_g:.2f} J/mol", ha="center", va="bottom", fontsize=11, fontweight="bold")
    ax.text((g_feed + g_eq) * 0.5, -0.33, rf"$TPDF_{{min}}$ = {tpdf:.4f}", ha="center", va="top", fontsize=11)
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(-0.55, 0.45)
    ax.set_yticks([])
    ax.set_xlabel(r"Molar Gibbs energy, $\hat g$ (J/mol)")
    for side in ("left", "right", "top"):
        ax.spines[side].set_visible(False)
    ax.grid(axis="x", color="#d7d7d7", linewidth=0.7, alpha=0.8)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=240, bbox_inches="tight")
    plt.close(fig)


def _save_partitioning(result: dict, out_path: Path) -> None:
    org_share = np.array([result["phases"]["organic"]["share_of_feed_pct"][sp] for sp in SPECIES], dtype=float)
    aq_share = np.array([result["phases"]["aqueous"]["share_of_feed_pct"][sp] for sp in SPECIES], dtype=float)

    fig = plt.figure(figsize=(8.4, 4.8), facecolor="white")
    ax = fig.add_axes([0.12, 0.16, 0.66, 0.70])
    y = np.arange(len(SPECIES), dtype=float)
    bar_h = 0.34

    ax.barh(y - 0.18, org_share, height=bar_h, color=PHASE_COLORS["organic"], alpha=0.92, label="Organic-rich phase")
    ax.barh(y + 0.18, aq_share, height=bar_h, color=PHASE_COLORS["aqueous"], alpha=0.92, label="Aqueous-rich phase")

    x_max = 104.0
    for i, (org_val, aq_val) in enumerate(zip(org_share, aq_share)):
        ax.text(min(org_val + 1.2, x_max - 0.5), y[i] - 0.18, f"{org_val:.1f}%", va="center", ha="left", fontsize=9.5, color=PHASE_COLORS["organic"])
        ax.text(min(aq_val + 1.2, x_max - 0.5), y[i] + 0.18, f"{aq_val:.1f}%", va="center", ha="left", fontsize=9.5, color=PHASE_COLORS["aqueous"])

    ax.set_xlim(0.0, x_max)
    ax.set_xticks([0, 20, 40, 60, 80, 100])
    ax.set_yticks(y)
    ax.set_yticklabels([SPECIES_LABELS[sp] for sp in SPECIES])
    ax.invert_yaxis()
    ax.set_xlabel("Component split between phases (%)")
    ax.set_title("Component partitioning between equilibrium phases", fontsize=15, pad=10)
    ax.grid(axis="x", color="#d7d7d7", linewidth=0.7, alpha=0.8)
    ax.set_axisbelow(True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(loc="center left", bbox_to_anchor=(1.01, 0.08), frameon=False)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=240, bbox_inches="tight")
    plt.close(fig)


def export_panels(out_dir: Path, model_key: str = DEFAULT_MODEL_KEY) -> list[Path]:
    result = _solve_case(model_key)
    paths = [
        out_dir / "ascani_case2_figiel2025_lle_compositions.png",
        out_dir / "ascani_case2_figiel2025_lle_phase_split.png",
        out_dir / "ascani_case2_figiel2025_lle_driving_force.png",
        out_dir / "ascani_case2_figiel2025_lle_partitioning.png",
    ]
    _save_composition_table(result, paths[0])
    _save_phase_split(result, paths[1])
    _save_driving_force(result, paths[2])
    _save_partitioning(result, paths[3])
    return paths


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Export standalone LLE case-study PNG panels.")
    parser.add_argument(
        "--outdir",
        type=Path,
        default=REPO_ROOT / "scripts" / "multiphase_model_analysis" / "output",
        help="Output directory for the PNG panels.",
    )
    parser.add_argument(
        "--model-key",
        default=DEFAULT_MODEL_KEY,
        help="Model key from ascani_case2_dataset_comparison._default_model_configs().",
    )
    return parser.parse_args()


def main() -> None:
    args = _parse_args()
    paths = export_panels(args.outdir, model_key=args.model_key)
    for path in paths:
        print(path)


if __name__ == "__main__":
    main()







