from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

SCRIPT_DIR = Path(__file__).resolve().parent
ANALYSIS_ROOT = SCRIPT_DIR.parent
REPO_ROOT = ANALYSIS_ROOT.parents[2]

if str(ANALYSIS_ROOT) not in sys.path:
    sys.path.insert(0, str(ANALYSIS_ROOT))
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

import _plot_common as common
import _model_overlay as overlay


DATA_PATH = SCRIPT_DIR / "data" / "water_contributions.csv"
FIGURE2_TOTALS_PATH = ANALYSIS_ROOT / "figure_2" / "data" / "water_comparisons.csv"
CONTRIBUTIONS = [
    ("hc", "Hard chain", "#9f9f9f"),
    ("disp", "Dispersion", "#5f5f5f"),
    ("assoc", "Association", "#111111"),
    ("dh", "Debye-Huckel", "#8c564b"),
    ("born", "Born", "#d8891c"),
]


def _model_values(ions: list[str], term_key: str) -> np.ndarray:
    out = np.empty(len(ions), dtype=float)
    for idx, ion in enumerate(ions):
        out[idx] = overlay.contribution_breakdown("advanced", ion, "water")[term_key]
    return out


def _plot_one(term_key: str, term_label: str, color: str, ions: list[str], paper: np.ndarray, model: np.ndarray) -> None:
    x = np.arange(len(ions), dtype=float)
    width = 0.34

    fig, ax = plt.subplots(figsize=(10.8, 5.8))
    ax.bar(
        x - width / 2.0,
        paper,
        width=width,
        color=color,
        edgecolor="black",
        linewidth=0.45,
        label=f"{term_label} (paper)",
    )
    ax.bar(
        x + width / 2.0,
        model,
        width=width,
        color=color,
        edgecolor="black",
        linewidth=0.45,
        hatch="////",
        alpha=0.75,
        label=f"{term_label} (pcsaft)",
    )

    values = np.vstack([paper, model])
    y_min = float(np.nanmin(values))
    y_max = float(np.nanmax(values))

    ax.axhline(0.0, color="black", linewidth=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(ions)
    ax.set_ylabel(r"Contribution to $\Delta G_{\mathrm{hyd},i}^{\infty}$ / kJ mol$^{-1}$")
    ax.set_title(f"Bulow 2020 Part I Figure 3: {term_label}")
    ax.set_ylim(y_min - 30.0, y_max + 30.0)
    ax.grid(axis="y", alpha=0.22)
    common.annotate_percent_deltas(ax, x + width / 2.0, paper, model, fontsize=7)
    common.add_percent_note(ax)
    ax.legend(frameon=True)

    output_path = SCRIPT_DIR / f"figure_3_{term_key}.png"
    common.save_figure(fig, output_path)
    plt.close(fig)
    print(f"Wrote {output_path}", flush=True)


def main() -> None:
    common.configure_style()
    frame = common.load_indexed_csv(DATA_PATH)
    totals = common.load_indexed_csv(FIGURE2_TOTALS_PATH)
    ions = list(frame.columns)
    for term_key, term_label, color in CONTRIBUTIONS:
        print(f"Computing Figure 3 {term_key}", flush=True)
        if term_key == "dh":
            paper = totals.values("advanced", ions)
            for base_key in ("hc", "disp", "assoc", "born"):
                paper = paper - frame.values(base_key, ions)
        else:
            paper = frame.values(term_key, ions)
        model = _model_values(ions, term_key)
        _plot_one(term_key, term_label, color, ions, paper, model)


if __name__ == "__main__":
    main()
