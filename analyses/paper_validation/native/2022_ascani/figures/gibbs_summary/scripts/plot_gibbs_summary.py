from __future__ import annotations

import csv
import sys
from pathlib import Path

for _candidate in Path(__file__).resolve().parents:
    if (_candidate / "scripts" / "plot_outputs.py").is_file():
        if str(_candidate) not in sys.path:
            sys.path.insert(0, str(_candidate))
        break
else:
    raise ModuleNotFoundError("Could not locate repo root containing scripts/plot_outputs.py")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from scripts.plot_outputs import figure_output_path, save_plot_figure


def _read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _style() -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Serif",
            "font.size": 10,
            "axes.linewidth": 1.0,
            "axes.grid": True,
            "grid.color": "0.88",
            "grid.linewidth": 0.6,
            "legend.frameon": False,
        }
    )


def main() -> None:
    _style()
    rows = _read_rows(figure_output_path(__file__, "data/gibbs_comparison.csv"))
    abs_rows = [row for row in rows if row["quantity"] in {"g_hat_feed_j_per_mol", "g_hat_eq_j_per_mol"}]
    delta_row = next(row for row in rows if row["quantity"] == "delta_g_hat_j_per_mol")

    fig, (ax_abs, ax_delta) = plt.subplots(1, 2, figsize=(8.8, 4.0), gridspec_kw={"width_ratios": [1.4, 1.0]})
    x = np.arange(len(abs_rows))
    width = 0.34
    paper_abs = [float(row["paper"]) for row in abs_rows]
    current_abs = [float(row["current_native_ln_fugacity_basis_j_per_mol"]) for row in abs_rows]
    ax_abs.bar(x - width / 2, paper_abs, width, color="white", edgecolor="0.25", linewidth=1.1, label="paper")
    ax_abs.bar(x + width / 2, current_abs, width, color="black", label="current")
    ax_abs.set_xticks(x, [r"$\hat{g}_{feed}$", r"$\hat{g}_{eq}$"])
    ax_abs.set_ylabel(r"$\hat{g}$ / J mol$^{-1}$")
    ax_abs.set_title("Gibbs objective level", fontsize=12)
    ax_abs.legend(loc="lower left", fontsize=8)

    ax_delta.bar([0 - width / 2], [float(delta_row["paper"])], width, color="white", edgecolor="0.25", linewidth=1.1)
    ax_delta.bar([0 + width / 2], [float(delta_row["current_native_ln_fugacity_basis_j_per_mol"])], width, color="black")
    ax_delta.axhline(0.0, color="0.25", linewidth=0.8)
    ax_delta.set_xticks([0], [r"$\Delta\hat{g}$"])
    ax_delta.set_ylabel(r"J mol$^{-1}$")
    ax_delta.set_title("Split driving force", fontsize=12)
    fig.suptitle("Ascani 2022 Gibbs comparison", fontsize=13)
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.94))
    output = figure_output_path(__file__, "gibbs_summary.png")
    save_plot_figure(fig, output, dpi=300, svg_companion=True)
    plt.close(fig)
    print(f"[write] {output}")


if __name__ == "__main__":
    main()
