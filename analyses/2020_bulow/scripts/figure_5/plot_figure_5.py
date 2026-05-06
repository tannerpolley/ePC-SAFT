from __future__ import annotations

import argparse
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


def _paper_values(frame, ion: str) -> np.ndarray:
    out = np.empty(len(TERMS), dtype=float)
    for idx, (term_key, _, _) in enumerate(TERMS):
        if term_key == "dh":
            out[idx] = 0.0
        else:
            row_key = "Born" if term_key == "Born" else term_key
            value = frame.scalar(row_key, ion)
            out[idx] = value
    return out


def _model_values(solvent: str, ion: str, d_born_mode: int | None = None) -> np.ndarray:
    values = overlay.transfer_breakdown("advanced", ion, solvent, basis="mu", d_born_mode=d_born_mode)
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


def _plot_one(
    panel_tag: str,
    solvent: str,
    term_key: str,
    term_label: str,
    ions: list[str],
    paper: np.ndarray,
    model: np.ndarray,
    suffix: str = "",
) -> None:
    x = np.arange(len(ions), dtype=float)
    width = 0.34
    fig, ax = plt.subplots(figsize=(10.8, 5.8))
    ion_color_map = dict(IONS)
    for idx, ion in enumerate(ions):
        color = ion_color_map[ion]
        paper_bars = ax.bar(
            x[idx] - width / 2.0,
            paper[idx],
            width=width,
            color=color,
            edgecolor="black",
            linewidth=0.45,
            label=f"{ion} (paper)",
        )
        model_bars = ax.bar(
            x[idx] + width / 2.0,
            model[idx],
            width=width,
            color=color,
            edgecolor="black",
            linewidth=0.45,
            hatch="////",
            alpha=0.75,
            label=f"{ion} (epcsaft)",
        )
        common.annotate_bar_values(ax, paper_bars, fontsize=6)
        common.annotate_bar_values(ax, model_bars, fontsize=6)

    stacked = np.vstack([paper, model])
    ax.axhline(0.0, color="black", linewidth=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(ions)
    ax.set_ylabel(r"Contribution to $\Delta G_{\mathrm{tr},i}^{\infty}$ / kJ mol$^{-1}$")
    ax.set_title(f"Bulow 2020 Part I Figure 5{panel_tag}: Water to {solvent}, {term_label} ($\\mu$ basis)")
    common.set_strict_bar_ylim(ax, stacked, step=5.0, top_pad_frac=0.18, bottom_pad_frac=0.12)
    ax.grid(axis="y", alpha=0.22)
    common.annotate_percent_deltas(ax, x + width / 2.0, paper, model, xs_ref=x - width / 2.0)
    common.add_percent_note(ax)
    ax.legend(frameon=True)

    output_path = SCRIPT_DIR / f"figure_5{panel_tag}_{term_key.lower()}{suffix}.png"
    common.save_figure(fig, output_path)
    plt.close(fig)
    print(f"Wrote {output_path}", flush=True)


def _plot_total(
    panel_tag: str,
    solvent: str,
    ions: list[str],
    figure5_paper: np.ndarray,
    figure5_model_eos: np.ndarray,
    d_born_mode: int | None = None,
    suffix: str = "",
) -> None:
    frame4 = common.load_indexed_csv(FIGURE4_DATA[solvent])
    figure4_paper = frame4.values("advanced", ions)
    figure4_model = np.asarray(
        [overlay.transfer_total("advanced", ion, solvent, d_born_mode=d_born_mode) for ion in ions], dtype=float
    )

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
            (offsets[1], figure5_model_eos[idx], color, "black", "////", 0.75, "Figure 5 EOS-consistent sum (epcsaft)"),
            (offsets[2], figure4_paper[idx], "white", color, None, 1.0, "Figure 4 total (paper)"),
            (offsets[3], figure4_model[idx], "white", color, "////", 1.0, "Figure 4 total (epcsaft)"),
        ]
        for offset, value, facecolor, edgecolor, hatch, alpha, label in bar_specs:
            show_label = label if label not in legend_added else None
            bars = ax.bar(
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
            common.annotate_bar_values(ax, bars, fontsize=6)
            if show_label is not None:
                legend_added.add(label)

    stacked = np.vstack([figure5_paper, figure5_model_eos, figure4_paper, figure4_model])
    y_min = float(np.nanmin(stacked))
    y_max = float(np.nanmax(stacked))
    pad = max(2.0, 0.12 * max(abs(y_min), abs(y_max), 1.0))

    ax.axhline(0.0, color="black", linewidth=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(ions)
    ax.set_ylabel(r"$\Delta G_{\mathrm{tr},i}^{\infty}$ / kJ mol$^{-1}$ (EOS total / EOS-sum check)")
    ax.set_title(
        f"Bulow 2020 Part I Figure 5{panel_tag}: Water to {solvent}, summed contribution bars vs Figure 4 total"
    )
    ax.set_ylim(y_min - pad, y_max + 2.7 * pad)
    ax.grid(axis="y", alpha=0.22)
    common.annotate_percent_deltas(ax, x + offsets[1], figure5_paper, figure5_model_eos, xs_ref=x + offsets[0])
    ax.text(
        0.99,
        0.12,
        "Colored bars: Figure 5 sums\nPaper sum uses published bars; epcsaft sum uses EOS-consistent lnfug contributions\nWhite bars: Figure 4 advanced totals",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=8,
        bbox={"facecolor": "white", "edgecolor": "0.5", "alpha": 0.9, "boxstyle": "round,pad=0.2"},
    )
    common.add_percent_note(ax)
    ax.legend(ncol=2, frameon=True)

    output_path = SCRIPT_DIR / f"figure_5{panel_tag}_total{suffix}.png"
    common.save_figure(fig, output_path)
    plt.close(fig)
    print(f"Wrote {output_path}", flush=True)


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate Bulow 2020 Figure 5 comparison plots.")
    parser.add_argument(
        "--d-born-mode",
        type=int,
        choices=[0, 1],
        default=None,
        help="Optional override for 2020_Bulow born_model.d_Born_mode.",
    )
    parser.add_argument(
        "--suffix",
        type=str,
        default="",
        help="Optional filename suffix to avoid overwriting the default Figure 5 outputs.",
    )
    return parser.parse_args()


def main() -> None:
    args = _parse_args()
    common.configure_style()
    for panel_tag, solvent, data_path in DATASETS:
        print(f"Computing Figure 5{panel_tag} ({solvent})", flush=True)
        frame = common.load_indexed_csv(data_path)
        ions = [ion for ion, _ in IONS]
        paper_map = {ion: _paper_values(frame, ion) for ion in ions}
        model_map = {ion: _model_values(solvent, ion, d_born_mode=args.d_born_mode) for ion in ions}
        for term_idx, (term_key, term_label, _) in enumerate(TERMS):
            paper = np.asarray([paper_map[ion][term_idx] for ion in ions], dtype=float)
            model = np.asarray([model_map[ion][term_idx] for ion in ions], dtype=float)
            _plot_one(panel_tag, solvent, term_key, term_label, ions, paper, model, suffix=args.suffix)
        paper_total = np.asarray([np.sum(paper_map[ion]) for ion in ions], dtype=float)
        model_total_eos = np.asarray(
            [
                sum(
                    overlay.transfer_breakdown("advanced", ion, solvent, basis="lnfug", d_born_mode=args.d_born_mode)[
                        key
                    ]
                    for key in ("hc", "disp", "assoc", "dh", "born")
                )
                for ion in ions
            ],
            dtype=float,
        )
        _plot_total(
            panel_tag,
            solvent,
            ions,
            paper_total,
            model_total_eos,
            d_born_mode=args.d_born_mode,
            suffix=args.suffix,
        )


if __name__ == "__main__":
    main()
