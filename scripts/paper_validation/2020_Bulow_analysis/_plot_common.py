from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


FIG_DPI = 300


def configure_style() -> None:
    plt.rcParams.update(
        {
            "font.size": 10,
            "axes.titlesize": 12,
            "axes.labelsize": 11,
            "legend.fontsize": 9,
            "xtick.labelsize": 10,
            "ytick.labelsize": 10,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "savefig.facecolor": "white",
        }
    )


class Table:
    def __init__(self, columns: list[str], rows: dict[str, dict[str, float]]) -> None:
        self.columns = columns
        self.rows = rows

    def values(self, row_key: str, columns: list[str] | None = None) -> np.ndarray:
        keys = self.columns if columns is None else columns
        row = self.rows[row_key]
        return np.asarray([row[key] for key in keys], dtype=float)

    def scalar(self, row_key: str, column: str) -> float:
        return float(self.rows[row_key][column])

    @property
    def index(self) -> list[str]:
        return list(self.rows.keys())


def load_indexed_csv(path: Path) -> Table:
    with path.open("r", newline="", encoding="utf-8-sig") as handle:
        reader = csv.reader(handle)
        raw_rows = [row for row in reader if row]

    header = [str(value).strip() for value in raw_rows[0][1:]]
    rows: dict[str, dict[str, float]] = {}
    for raw_row in raw_rows[1:]:
        row_key = str(raw_row[0]).strip()
        rows[row_key] = {
            column: float(value)
            for column, value in zip(header, raw_row[1:], strict=False)
        }

    return Table(header, rows)


def save_figure(fig: plt.Figure, output_path: Path) -> None:
    fig.tight_layout()
    fig.savefig(output_path, dpi=FIG_DPI, bbox_inches="tight")


def percent_delta(model_value: float, paper_value: float) -> float:
    denom = abs(float(paper_value))
    if denom <= 1.0e-12:
        return float("nan")
    return 100.0 * (float(model_value) - float(paper_value)) / denom


def annotate_percent_deltas(
    ax: plt.Axes,
    xs: np.ndarray,
    paper: np.ndarray,
    model: np.ndarray,
    *,
    fontsize: int = 7,
    rotation: float = 90.0,
) -> None:
    y_min, y_max = ax.get_ylim()
    span = max(y_max - y_min, 1.0)
    y_pad = 0.018 * span

    for x_val, paper_val, model_val in zip(xs, paper, model, strict=False):
        if not (np.isfinite(paper_val) and np.isfinite(model_val)):
            continue
        pct = percent_delta(float(model_val), float(paper_val))
        if not np.isfinite(pct):
            continue
        y_text = float(model_val) + (y_pad if model_val >= 0.0 else -y_pad)
        va = "bottom" if model_val >= 0.0 else "top"
        ax.text(
            float(x_val),
            y_text,
            f"{pct:+.1f}%",
            ha="center",
            va=va,
            fontsize=fontsize,
            rotation=rotation,
            color="black",
            clip_on=False,
        )


def add_percent_note(ax: plt.Axes, *, xpos: float = 0.99, ypos: float = 0.01) -> None:
    ax.text(
        xpos,
        ypos,
        r"% labels: $(\mathrm{pcsaft} - \mathrm{paper})/|\mathrm{paper}|$",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=8,
        bbox={"facecolor": "white", "edgecolor": "0.5", "alpha": 0.9, "boxstyle": "round,pad=0.2"},
    )
