from __future__ import annotations

from pathlib import Path
from typing import Iterable

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from epcsaft import ePCSAFTMixture
from epcsaft.regression import _debug_native_pure_neutral_objective
from scripts import plot_outputs
from tests.helpers.regression_cases import _methane_like_records
from tests.helpers.regression_cases import _minimal_neutral_metadata


def hydrocarbon_basis_mixture() -> ePCSAFTMixture:
    params = {
        "m": np.asarray([1.0, 1.6069, 2.0020]),
        "s": np.asarray([3.7039, 3.5206, 3.6184]),
        "e": np.asarray([150.03, 191.42, 208.11]),
        "k_ij": np.asarray(
            [
                [0.0, 3.0e-4, 1.15e-2],
                [3.0e-4, 0.0, 5.10e-3],
                [1.15e-2, 5.10e-3, 0.0],
            ]
        ),
    }
    return ePCSAFTMixture.from_params(params, species=["Methane", "Ethane", "Propane"])


def methanol_cyclohexane_mixture() -> ePCSAFTMixture:
    params = {
        "MW": np.asarray([32.042e-3, 84.147e-3]),
        "m": np.asarray([1.5255, 2.5303]),
        "s": np.asarray([3.2300, 3.8499]),
        "e": np.asarray([188.90, 278.11]),
        "e_assoc": np.asarray([2899.5, 0.0]),
        "vol_a": np.asarray([0.035176, 0.0]),
        "assoc_scheme": ["2B", None],
        "k_ij": np.asarray([[0.0, 0.051], [0.051, 0.0]]),
        "z": np.asarray([0.0, 0.0]),
        "dielc": np.asarray([33.05, 2.02]),
    }
    return ePCSAFTMixture.from_params(params, species=["Methanol", "Cyclohexane"])


def assert_plot_with_data(path: Path) -> None:
    csv_path = path.parent / "data" / f"{path.stem}_plot_data.csv"
    assert path.exists()
    assert csv_path.exists()
    assert csv_path.stat().st_size > 0


def save_comparison_plot(
    filename: str,
    title: str,
    labels: list[str],
    actual: np.ndarray,
    expected: np.ndarray,
    *,
    category: Iterable[str],
    ylabel: str = "Value",
    relative_error: bool = True,
) -> Path:
    actual = np.asarray(actual, dtype=float)
    expected = np.asarray(expected, dtype=float)
    x = np.arange(len(labels), dtype=float)
    fig_height = 5.6 if relative_error else 4.4
    if relative_error:
        fig, axes = plt.subplots(2, 1, figsize=(max(7.0, len(labels) * 0.68), fig_height), height_ratios=[3, 1.35])
        ax, err_ax = axes
    else:
        fig, ax = plt.subplots(figsize=(max(7.0, len(labels) * 0.68), fig_height))
        err_ax = None

    ax.bar(x - 0.18, actual, width=0.36, label="Actual")
    ax.bar(x + 0.18, expected, width=0.36, label="Expected")
    ax.set_title(title)
    ax.set_ylabel(ylabel)
    ax.set_xticks(x, labels, rotation=35, ha="right")
    ax.legend()
    ax.axhline(0.0, color="0.82", linewidth=0.8)

    if err_ax is not None:
        scale = np.maximum(np.abs(expected), 1.0e-30)
        err_ax.bar(x, (actual - expected) / scale, width=0.5, label="Relative error")
        err_ax.axhline(0.0, color="0.25", linewidth=0.8)
        err_ax.set_ylabel("Rel. err.")
        err_ax.set_xticks(x, labels, rotation=35, ha="right")

    output_path = plot_outputs.test_plot_path(__file__, filename, category=category)
    try:
        plot_outputs.save_plot_figure(fig, output_path, dpi=120)
    finally:
        plt.close(fig)
    assert_plot_with_data(output_path)
    return output_path


def save_parity_plot(
    filename: str,
    title: str,
    labels: list[str],
    actual: np.ndarray,
    expected: np.ndarray,
    *,
    category: Iterable[str],
    xlabel: str = "Expected",
    ylabel: str = "Actual",
) -> Path:
    actual = np.asarray(actual, dtype=float)
    expected = np.asarray(expected, dtype=float)
    finite = np.isfinite(actual) & np.isfinite(expected)
    plot_actual = actual[finite]
    plot_expected = expected[finite]

    fig, axes = plt.subplots(1, 2, figsize=(12.0, 4.8), width_ratios=[2.2, 1.5])
    ax, err_ax = axes
    ax.scatter(plot_expected, plot_actual, label="Comparison points")
    if plot_expected.size:
        lo = float(min(np.min(plot_expected), np.min(plot_actual)))
        hi = float(max(np.max(plot_expected), np.max(plot_actual)))
        pad = max((hi - lo) * 0.08, 1.0e-12)
        ax.plot([lo - pad, hi + pad], [lo - pad, hi + pad], color="0.25", linewidth=1.0, label="Parity")
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.legend()

    x = np.arange(len(labels), dtype=float)
    scale = np.maximum(np.abs(expected), 1.0e-30)
    err_ax.bar(x, (actual - expected) / scale, width=0.55)
    err_ax.axhline(0.0, color="0.25", linewidth=0.8)
    err_ax.set_title("Relative error")
    err_ax.set_xticks(x, labels, rotation=45, ha="right")

    output_path = plot_outputs.test_plot_path(__file__, filename, category=category)
    try:
        plot_outputs.save_plot_figure(fig, output_path, dpi=120)
    finally:
        plt.close(fig)
    assert_plot_with_data(output_path)
    return output_path


def append_payload_rows(rows: list[dict[str, object]], prefix: str, payload: dict[str, object]) -> None:
    total = np.asarray(payload["total"], dtype=float)
    term_arrays = {key: np.asarray(value, dtype=float) for key, value in payload["terms"].items()}
    if total.ndim == 0:
        rows.append(
            {
                "label": prefix,
                "terms": {key: float(value) for key, value in payload["terms"].items()},
                "total": float(total),
            }
        )
        return
    for index, total_value in enumerate(total.tolist()):
        rows.append(
            {
                "label": f"{prefix}[{index}]",
                "terms": {key: float(value[index]) for key, value in term_arrays.items()},
                "total": float(total_value),
            }
        )


def save_contribution_closure_plot(
    filename: str,
    title: str,
    rows: list[dict[str, object]],
    *,
    category: Iterable[str],
) -> Path:
    labels = [str(row["label"]) for row in rows]
    term_names = sorted({name for row in rows for name in row["terms"]})
    x = np.arange(len(rows), dtype=float)

    fig, axes = plt.subplots(2, 1, figsize=(max(8.0, len(rows) * 0.72), 7.0), height_ratios=[3.0, 1.4])
    ax, err_ax = axes
    positive_bottom = np.zeros(len(rows), dtype=float)
    negative_bottom = np.zeros(len(rows), dtype=float)
    for term_name in term_names:
        values = np.asarray([float(row["terms"].get(term_name, 0.0)) for row in rows], dtype=float)
        bottoms = np.where(values >= 0.0, positive_bottom, negative_bottom)
        ax.bar(x, values, bottom=bottoms, width=0.58, label=term_name)
        positive_bottom += np.where(values >= 0.0, values, 0.0)
        negative_bottom += np.where(values < 0.0, values, 0.0)

    totals = np.asarray([float(row["total"]) for row in rows], dtype=float)
    term_sums = np.asarray([sum(float(value) for value in row["terms"].values()) for row in rows], dtype=float)
    ax.scatter(x, totals, marker="x", color="black", label="Reported total", zorder=5)
    ax.set_title(title)
    ax.set_ylabel("Contribution value")
    ax.set_xticks(x, labels, rotation=35, ha="right")
    ax.axhline(0.0, color="0.82", linewidth=0.8)
    ax.legend(ncol=min(3, max(1, len(term_names))), fontsize="small")

    scale = np.maximum(np.abs(totals), 1.0e-30)
    err_ax.bar(x, (term_sums - totals) / scale, width=0.55, label="Closure error")
    err_ax.axhline(0.0, color="0.25", linewidth=0.8)
    err_ax.set_ylabel("Rel. closure")
    err_ax.set_xticks(x, labels, rotation=35, ha="right")

    output_path = plot_outputs.test_plot_path(__file__, filename, category=category)
    try:
        plot_outputs.save_plot_figure(fig, output_path, dpi=120)
    finally:
        plt.close(fig)
    assert_plot_with_data(output_path)
    return output_path


def save_contribution_term_breakdown_plot(
    filename: str,
    title: str,
    rows: list[dict[str, object]],
    *,
    category: Iterable[str],
) -> Path:
    labels = [str(row["label"]) for row in rows]
    term_names = sorted({name for row in rows for name in row["terms"]})
    x = np.arange(len(rows), dtype=float)
    width = min(0.75 / max(len(term_names), 1), 0.18)
    offsets = (np.arange(len(term_names), dtype=float) - (len(term_names) - 1) / 2.0) * width

    fig, ax = plt.subplots(figsize=(max(9.0, len(rows) * 0.8), 4.8))
    for offset, term_name in zip(offsets, term_names, strict=False):
        values = np.asarray([float(row["terms"].get(term_name, 0.0)) for row in rows], dtype=float)
        ax.bar(x + offset, values, width=width, label=term_name)

    ax.scatter(x, [float(row["total"]) for row in rows], marker="x", color="black", label="Reported total", zorder=5)
    ax.set_title(title)
    ax.set_ylabel("Term contribution")
    ax.set_xticks(x, labels, rotation=35, ha="right")
    ax.axhline(0.0, color="0.82", linewidth=0.8)
    ax.legend(ncol=min(4, max(1, len(term_names))), fontsize="small")

    output_path = plot_outputs.test_plot_path(__file__, filename, category=category)
    try:
        plot_outputs.save_plot_figure(fig, output_path, dpi=120)
    finally:
        plt.close(fig)
    assert_plot_with_data(output_path)
    return output_path


def finite_difference_gradient_values() -> tuple[np.ndarray, np.ndarray]:
    theta = {"m": 1.05, "s": 3.68, "e": 151.0}

    def objective_at(m: float, s: float, e: float) -> float:
        debug = _debug_native_pure_neutral_objective(
            _methane_like_records(),
            "Methane",
            assoc_scheme="",
            fixed_parameters=_minimal_neutral_metadata(16.043e-3),
            initial_guess=theta,
            x={"m": m, "s": s, "e": e},
        )
        return float(debug["objective"])

    debug = _debug_native_pure_neutral_objective(
        _methane_like_records(),
        "Methane",
        assoc_scheme="",
        fixed_parameters=_minimal_neutral_metadata(16.043e-3),
        initial_guess=theta,
        x=theta,
    )
    exact = np.asarray(debug["gradient"], dtype=float)
    eps = np.asarray([1.0e-6, 1.0e-6, 1.0e-5], dtype=float)
    fd = np.empty(3, dtype=float)
    base = np.asarray([theta["m"], theta["s"], theta["e"]], dtype=float)
    for i in range(3):
        forward = base.copy()
        backward = base.copy()
        forward[i] += eps[i]
        backward[i] -= eps[i]
        fd[i] = (objective_at(*forward) - objective_at(*backward)) / (2.0 * eps[i])
    return exact, fd
