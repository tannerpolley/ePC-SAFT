from __future__ import annotations

import csv
import math
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[3]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts.plot_outputs import paper_validation_output_path


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
    output_path = paper_validation_output_path(output_path)
    fig.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
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
    xs_ref: np.ndarray | None = None,
    fontsize: int = 8,
    rotation: float = 0.0,
) -> None:
    y_min, y_max = ax.get_ylim()
    span = max(y_max - y_min, 1.0)
    y_pad = 0.05 * span

    if xs_ref is None:
        xs_ref = xs

    for x_val, x_ref_val, paper_val, model_val in zip(xs, xs_ref, paper, model, strict=False):
        if not (np.isfinite(paper_val) and np.isfinite(model_val)):
            continue
        pct = percent_delta(float(model_val), float(paper_val))
        if not np.isfinite(pct):
            continue
        x_mid = 0.5 * (float(x_val) + float(x_ref_val))
        top_value = max(float(paper_val), float(model_val))
        bottom_value = min(float(paper_val), float(model_val))
        if paper_val < 0.0 and model_val < 0.0:
            y_text = bottom_value - y_pad
            va = "top"
        else:
            y_text = top_value + y_pad
            va = "bottom"
        ax.text(
            x_mid,
            y_text,
            f"{pct:+.1f}%",
            ha="center",
            va=va,
            fontsize=fontsize,
            rotation=rotation,
            color="black",
            fontweight="semibold",
            clip_on=False,
        )


def annotate_bar_values(
    ax: plt.Axes,
    bars,
    *,
    fmt: str = "{:.1f}",
    fontsize: int = 8,
    rotation: float = 0.0,
    color: str = "0.2",
) -> None:
    y_min, y_max = ax.get_ylim()
    span = max(y_max - y_min, 1.0)
    y_pad = 0.012 * span

    for bar in bars:
        height = float(bar.get_height())
        if not np.isfinite(height):
            continue
        x_val = float(bar.get_x() + bar.get_width() / 2.0)
        y_text = height + (y_pad if height >= 0.0 else -y_pad)
        va = "bottom" if height >= 0.0 else "top"
        ax.text(
            x_val,
            y_text,
            fmt.format(height),
            ha="center",
            va=va,
            fontsize=fontsize,
            rotation=rotation,
            color=color,
            fontweight="semibold",
            alpha=0.85,
            clip_on=False,
        )


def _round_up_to_multiple(value: float, step: float) -> float:
    return float(step * math.ceil(value / step))


def _round_down_to_multiple(value: float, step: float) -> float:
    return float(step * math.floor(value / step))


def set_strict_bar_ylim(
    ax: plt.Axes,
    values: np.ndarray,
    *,
    step: float = 5.0,
    top_pad_frac: float = 0.12,
    bottom_pad_frac: float = 0.08,
) -> tuple[float, float]:
    finite = np.asarray(values, dtype=float)
    finite = finite[np.isfinite(finite)]
    if finite.size == 0:
        lo, hi = -step, step
        ax.set_ylim(lo, hi)
        return lo, hi

    raw_min = float(np.min(finite))
    raw_max = float(np.max(finite))
    span = max(raw_max - raw_min, step)

    padded_min = raw_min - bottom_pad_frac * span
    if raw_max <= 0.0:
        padded_max = 0.0
    else:
        padded_max = raw_max + top_pad_frac * span

    lo = _round_down_to_multiple(padded_min, step)
    hi = 0.0 if padded_max == 0.0 else _round_up_to_multiple(padded_max, step)
    if lo == hi:
        hi = lo + step

    ax.set_ylim(lo, hi)
    return lo, hi


def add_percent_note(ax: plt.Axes, *, xpos: float = 0.99, ypos: float = 0.01) -> None:
    ax.text(
        xpos,
        ypos,
        r"% labels: $(\mathrm{epcsaft} - \mathrm{paper})/|\mathrm{paper}|$",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=8,
        bbox={"facecolor": "white", "edgecolor": "0.5", "alpha": 0.9, "boxstyle": "round,pad=0.2"},
    )

