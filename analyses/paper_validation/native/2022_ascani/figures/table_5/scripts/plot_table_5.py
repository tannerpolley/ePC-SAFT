from __future__ import annotations

import csv
import math
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

from scripts.plot_outputs import figure_output_path, save_plot_figure


def _read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _float_or_nan(value: str) -> float:
    if value == "":
        return math.nan
    return float(value)


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
    rows = _read_rows(figure_output_path(__file__, "data/table_5_fugacity_comparison.csv"))
    components = ["H2O", "Butanol", "NaCl", "KCl"]
    y_base = {component: idx for idx, component in enumerate(reversed(components))}

    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    for component in components:
        component_rows = [row for row in rows if row["component"] == component]
        paper = float(component_rows[0]["paper_ln_f_bar"])
        y = y_base[component]
        ax.scatter(paper, y, marker="o", s=64, facecolor="white", edgecolor="0.25", linewidth=1.2, label="paper" if component == "H2O" else None)
        for phase, marker, color, offset, label in (
            ("org", "s", "black", -0.12, "current org"),
            ("aq", "^", "#1f77b4", 0.12, "current aq"),
        ):
            row = next(item for item in component_rows if item["phase"] == phase)
            current = _float_or_nan(row["current_ln_f_bar"])
            if math.isfinite(current):
                ax.scatter(current, y + offset, marker=marker, s=48, color=color, label=label if component == "H2O" else None)
    ax.set_yticks([y_base[component] for component in components], [r"$H_2O$", "1-butanol", r"$NaCl$ mean", r"$KCl$ mean"])
    ax.set_xlabel(r"$\ln(f/\mathrm{bar})$")
    ax.set_title("Ascani 2022 Table 5 fugacity comparison", fontsize=13)
    ax.legend(loc="lower right", fontsize=8)
    output = figure_output_path(__file__, "table_5_fugacity.png")
    save_plot_figure(fig, output, dpi=300, svg_companion=True)
    plt.close(fig)
    print(f"[write] {output}")


if __name__ == "__main__":
    main()
