from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import matplotlib
import numpy as np
import pandas as pd

from _paths import ANALYSIS_DIR
from generate_rezaee_2026_figure_comparison_data import FIGURE_SPECS

matplotlib.use("Agg")
import matplotlib.pyplot as plt

RESULTS_DIR = ANALYSIS_DIR / "results" / "reaction_equilibrium"
FIGURES_DIR = ANALYSIS_DIR / "figures"
SUMMARY_JSON = RESULTS_DIR / "rezaee_2026_package_figure_comparison_summary.json"
REPORT_MD = RESULTS_DIR / "rezaee_2026_package_figure_comparison.md"
SVG_METADATA = {"Date": "2026-05-17T00:00:00"}
SVG_HASHSALT = "rezaee_2026_package_figure_comparison"


def _jsonable(value: Any) -> Any:
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, dict):
        return {str(key): _jsonable(item) for key, item in value.items()}
    if isinstance(value, list):
        return [_jsonable(item) for item in value]
    return value


def _aard_pct(rows: pd.DataFrame) -> float:
    x = rows["experimental_value"].to_numpy(dtype=float)
    y = rows["calculated_value"].to_numpy(dtype=float)
    return float((100.0 * np.abs(y - x) / np.clip(x, 1.0e-300, None)).mean())


def _configure_axes(ax: plt.Axes, axis_max: float, x_label: str, y_label: str) -> None:
    ax.set_xlim(-0.05 * axis_max, 1.05 * axis_max)
    ax.set_ylim(-0.05 * axis_max, 1.05 * axis_max)
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.grid(True, which="major", color="#d8d8d8", linewidth=0.8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_linewidth(1.1)
    ax.spines["bottom"].set_linewidth(1.1)
    ax.tick_params(direction="out", length=5, width=1.0)


def _normalize_svg_whitespace(path: Path) -> None:
    text = path.read_text(encoding="utf-8")
    path.write_text("\n".join(line.rstrip() for line in text.splitlines()) + "\n", encoding="utf-8")


def _render_figure(figure_id: str, spec: dict[str, Any]) -> dict[str, Any]:
    figure_folder = FIGURES_DIR / str(spec["folder"])
    data_path = figure_folder / "output" / "data" / f"{spec['folder']}_package_vs_paper_points.csv"
    rows = pd.read_csv(data_path)
    paper = rows.loc[rows["series"] == "digitized_paper_model"].copy()
    package = rows.loc[rows["series"] == "current_epcsaft_package_model"].copy()
    if paper.empty or package.empty:
        raise ValueError(f"{data_path} must contain both digitized paper and current package model series.")

    axis_max = float(spec["axis_max"])
    output_dir = figure_folder / "output"
    output_dir.mkdir(parents=True, exist_ok=True)

    plt.rcParams.update(
        {
            "font.family": "DejaVu Serif",
            "font.size": 12,
            "axes.labelsize": 13,
            "axes.titlesize": 13,
            "xtick.labelsize": 11,
            "ytick.labelsize": 11,
            "legend.fontsize": 10,
            "svg.hashsalt": SVG_HASHSALT,
        }
    )
    fig, ax = plt.subplots(figsize=(7.4, 5.6), constrained_layout=True)
    ax.plot([0.0, axis_max], [0.0, axis_max], color="#8f1d1d", linewidth=1.8, label="1:1 line", zorder=1)
    ax.scatter(
        paper["experimental_value"],
        paper["calculated_value"],
        s=56,
        marker="o",
        facecolor="#1f1f1f",
        edgecolor="white",
        linewidth=0.5,
        label=f"Digitized paper model (AARD {_aard_pct(paper):.1f}%)",
        zorder=3,
    )
    ax.scatter(
        package["experimental_value"],
        package["calculated_value"],
        s=74,
        marker="X",
        facecolor="#2f6fb4",
        edgecolor="#12395f",
        linewidth=0.7,
        label=f"Current epcsaft package (AARD {_aard_pct(package):.1f}%)",
        zorder=4,
    )
    _configure_axes(ax, axis_max, str(spec["x_label"]), str(spec["y_label"]))
    ax.set_title(f"Rezaee 2026 {spec['paper_label']}: package model vs digitized paper")
    ax.legend(loc="best", frameon=True, framealpha=0.94)

    png_path = output_dir / f"{spec['folder']}_package_vs_paper.png"
    svg_path = output_dir / f"{spec['folder']}_package_vs_paper.svg"
    fig.savefig(png_path, dpi=220)
    fig.savefig(svg_path, metadata=SVG_METADATA)
    _normalize_svg_whitespace(svg_path)
    plt.close(fig)

    return {
        "figure_id": figure_id,
        "paper_label": str(spec["paper_label"]),
        "caption": str(spec["caption"]),
        "data": str(data_path.relative_to(ANALYSIS_DIR)),
        "png": str(png_path.relative_to(ANALYSIS_DIR)),
        "svg": str(svg_path.relative_to(ANALYSIS_DIR)),
        "package_case_id": str(spec["package_case_id"]),
        "digitized_paper_points": int(len(paper)),
        "package_model_points": int(len(package)),
        "digitized_paper_aard_pct": _aard_pct(paper),
        "current_epcsaft_package_aard_pct": _aard_pct(package),
    }


def _write_report(summary: dict[str, Any]) -> None:
    lines = [
        "# Rezaee 2026 Package Figure Comparisons",
        "",
        "## Scope",
        "",
        "These figures overlay the digitized Rezaee 2026 paper-model points with the current in-worktree `epcsaft` package-model results. The package series is generated from `rezaee_2026_section32_equilibrium_replication_rows.csv`; downstream scripts are not executed.",
        "",
        "## Outputs",
        "",
    ]
    for entry in summary["figures"]:
        lines.extend(
            [
                f"- {entry['paper_label']}: digitized paper AARD `{entry['digitized_paper_aard_pct']:.3g}%`; current package AARD `{entry['current_epcsaft_package_aard_pct']:.3g}%`.",
                f"  Data: `{entry['data']}`",
                f"  PNG: `{entry['png']}`",
                f"  SVG: `{entry['svg']}`",
            ]
        )
    REPORT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    figures = [_render_figure(figure_id, spec) for figure_id, spec in FIGURE_SPECS.items()]
    summary = {
        "status": "package_figure_comparison_rendered",
        "figure_count": len(figures),
        "figures": figures,
    }
    SUMMARY_JSON.parent.mkdir(parents=True, exist_ok=True)
    SUMMARY_JSON.write_text(json.dumps(_jsonable(summary), indent=2) + "\n", encoding="utf-8")
    _write_report(summary)
    print(json.dumps(_jsonable(summary), indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
