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


DATASETS = [
    ("a", "ethanol", SCRIPT_DIR / "data" / "water-ethanol-contributions.csv"),
    ("b", "methanol", SCRIPT_DIR / "data" / "water-methanol-contributions.csv"),
]
FIGURE4_DATA = {
    "methanol": ANALYSIS_ROOT / "figure_4" / "data" / "water-methanol-comparison.csv",
    "ethanol": ANALYSIS_ROOT / "figure_4" / "data" / "water-ethanol-comparison.csv",
}
IONS = [("Na+", "#2b6cb0"), ("Cl-", "#c44e52"), ("I-", "#3a923a")]
TERMS = [
    ("hc", "Hard chain", "#9f9f9f"),
    ("disp", "Dispersion", "#5f5f5f"),
    ("assoc", "Association", "#111111"),
    ("dh", "Debye-Huckel", "#8c564b"),
    ("Born", "Born", "#d8891c"),
]


def _paper_values(frame, totals_frame, ion: str) -> np.ndarray:
    out = np.empty(len(TERMS), dtype=float)
    total_advanced = totals_frame.scalar("advanced", ion)
    running_sum = 0.0
    for idx, (term_key, _, _) in enumerate(TERMS):
        row_key = "Born" if term_key == "Born" else term_key
        if term_key == "dh":
            out[idx] = total_advanced - running_sum
        else:
            value = frame.scalar(row_key, ion)
            out[idx] = value
            running_sum += value
    return out


def _model_values(solvent: str, ion: str) -> np.ndarray:
    values = overlay.transfer_breakdown("advanced", ion, solvent)
    return np.asarray(
        [
            values["hc"],
            values["disp"],
            values["assoc"],
            values["dh"],
            values["born"],
        ],
        dtype=float,
    )


def _plot_one(panel_tag: str, solvent: str, term_key: str, term_label: str, ions: list[str], paper: np.ndarray, model: np.ndarray) -> None:
    x = np.arange(len(ions), dtype=float)
    width = 0.34
    fig, ax = plt.subplots(figsize=(10.8, 5.8))
    ion_color_map = dict(IONS)
    for idx, ion in enumerate(ions):
        color = ion_color_map[ion]
        ax.bar(
            x[idx] - width / 2.0,
            paper[idx],
            width=width,
            color=color,
            edgecolor="black",
            linewidth=0.45,
            label=f"{ion} (paper)",
        )
        ax.bar(
            x[idx] + width / 2.0,
            model[idx],
            width=width,
            color=color,
            edgecolor="black",
            linewidth=0.45,
            hatch="////",
            alpha=0.75,
            label=f"{ion} (pcsaft)",
        )

    stacked = np.vstack([paper, model])
    y_min = float(np.nanmin(stacked))
    y_max = float(np.nanmax(stacked))
    pad = max(2.0, 0.12 * max(abs(y_min), abs(y_max), 1.0))
    ax.axhline(0.0, color="black", linewidth=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(ions)
    ax.set_ylabel(r"Contribution to $\Delta G_{\mathrm{tr},i}^{\infty}$ / kJ mol$^{-1}$")
    ax.set_title(f"Bulow 2020 Part I Figure 5{panel_tag}: Water to {solvent}, {term_label}")
    ax.set_ylim(y_min - pad, y_max + 2.6 * pad)
    ax.grid(axis="y", alpha=0.22)
    common.annotate_percent_deltas(ax, x + width / 2.0, paper, model, fontsize=7)
    common.add_percent_note(ax)
    ax.legend(frameon=True)

    output_path = SCRIPT_DIR / f"figure_5{panel_tag}_{term_key.lower()}.png"
    common.save_figure(fig, output_path)
    plt.close(fig)
    print(f"Wrote {output_path}", flush=True)


def _plot_total(panel_tag: str, solvent: str, ions: list[str], figure5_paper: np.ndarray, figure5_model: np.ndarray) -> None:
    frame4 = common.load_indexed_csv(FIGURE4_DATA[solvent])
    figure4_paper = frame4.values("advanced", ions)
    figure4_model = np.asarray([overlay.transfer_total("advanced", ion, solvent) for ion in ions], dtype=float)

    x = np.arange(len(ions), dtype=float)
    width = 0.17
    offsets = np.array([-1.5, -0.5, 0.5, 1.5]) * width
    fig, ax = plt.subplots(figsize=(11.4, 6.0))
    ion_color_map = dict(IONS)
    legend_added: set[str] = set()

    for idx, ion in enumerate(ions):
        color = ion_color_map[ion]
        bar_specs = [
            (offsets[0], figure5_paper[idx], color, "black", None, 0.9, "Figure 5 sum (paper)"),
            (offsets[1], figure5_model[idx], color, "black", "////", 0.75, "Figure 5 sum (pcsaft)"),
            (offsets[2], figure4_paper[idx], "white", color, None, 1.0, "Figure 4 total (paper)"),
            (offsets[3], figure4_model[idx], "white", color, "////", 1.0, "Figure 4 total (pcsaft)"),
        ]
        for offset, value, facecolor, edgecolor, hatch, alpha, label in bar_specs:
            show_label = label if label not in legend_added else None
            ax.bar(
                x[idx] + offset,
                value,
                width=width,
                color=facecolor,
                edgecolor=edgecolor,
                linewidth=1.0 if facecolor == "white" else 0.45,
                hatch=hatch,
                alpha=alpha,
                label=show_label,
            )
            if show_label is not None:
                legend_added.add(label)

    stacked = np.vstack([figure5_paper, figure5_model, figure4_paper, figure4_model])
    y_min = float(np.nanmin(stacked))
    y_max = float(np.nanmax(stacked))
    pad = max(2.0, 0.12 * max(abs(y_min), abs(y_max), 1.0))

    ax.axhline(0.0, color="black", linewidth=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(ions)
    ax.set_ylabel(r"$\Delta G_{\mathrm{tr},i}^{\infty}$ / kJ mol$^{-1}$")
    ax.set_title(f"Bulow 2020 Part I Figure 5{panel_tag}: Water to {solvent}, summed contributions vs Figure 4 total")
    ax.set_ylim(y_min - pad, y_max + 2.7 * pad)
    ax.grid(axis="y", alpha=0.22)
    common.annotate_percent_deltas(ax, x + offsets[1], figure5_paper, figure5_model, fontsize=7)
    ax.text(
        0.99,
        0.12,
        "Colored bars: Figure 5 contribution sums\nWhite bars: Figure 4 advanced totals",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=8,
        bbox={"facecolor": "white", "edgecolor": "0.5", "alpha": 0.9, "boxstyle": "round,pad=0.2"},
    )
    common.add_percent_note(ax)
    ax.legend(ncol=2, frameon=True)

    output_path = SCRIPT_DIR / f"figure_5{panel_tag}_total.png"
    common.save_figure(fig, output_path)
    plt.close(fig)
    print(f"Wrote {output_path}", flush=True)


def main() -> None:
    common.configure_style()
    for panel_tag, solvent, data_path in DATASETS:
        print(f"Computing Figure 5{panel_tag} ({solvent})", flush=True)
        frame = common.load_indexed_csv(data_path)
        totals_frame = common.load_indexed_csv(FIGURE4_DATA[solvent])
        ions = [ion for ion, _ in IONS]
        paper_map = {ion: _paper_values(frame, totals_frame, ion) for ion in ions}
        model_map = {ion: _model_values(solvent, ion) for ion in ions}
        for term_idx, (term_key, term_label, _) in enumerate(TERMS):
            paper = np.asarray([paper_map[ion][term_idx] for ion in ions], dtype=float)
            model = np.asarray([model_map[ion][term_idx] for ion in ions], dtype=float)
            _plot_one(panel_tag, solvent, term_key, term_label, ions, paper, model)
        paper_total = np.asarray([np.sum(paper_map[ion]) for ion in ions], dtype=float)
        model_total = np.asarray([np.sum(model_map[ion]) for ion in ions], dtype=float)
        _plot_total(panel_tag, solvent, ions, paper_total, model_total)


if __name__ == "__main__":
    main()
