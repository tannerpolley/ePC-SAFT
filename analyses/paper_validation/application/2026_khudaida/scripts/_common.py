from __future__ import annotations

import csv
import math
import os
import platform
import sys


from pathlib import Path
import sys as _bootstrap_sys
from pathlib import Path as _BootstrapPath

for _candidate in _BootstrapPath(__file__).resolve().parents:
    if (_candidate / "scripts" / "plot_outputs.py").is_file():
        if str(_candidate) not in _bootstrap_sys.path:
            _bootstrap_sys.path.insert(0, str(_candidate))
        break
else:
    raise ModuleNotFoundError("Could not locate repo root containing scripts/plot_outputs.py")
from scripts.plot_outputs import REPO_ROOT
from typing import Iterable

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parent
ANALYSIS_ROOT = ROOT.parent
SOURCE_INPUT_ROOT = ANALYSIS_ROOT / "data" / "input"
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts._env import require_epcsaft_install
from scripts.plot_outputs import analysis_data_path, analysis_plot_set_dir, paper_validation_output_path, save_plot_figure

require_epcsaft_install()


def _fast_machine() -> str:
    return os.environ.get("PROCESSOR_ARCHITECTURE", "AMD64")


platform.machine = _fast_machine

import scripts._epcsaft_oop as pcs
import epcsaft
from scripts.validation.equilibrium_core.confidence import explicit_to_formula
from scripts.validation.equilibrium_core.confidence import formula_to_explicit
from epcsaft.parameters import get_prop_dict

P_REF = 1.0e5
MODEL_SOLVE_MAX_ITERATIONS = 30
PARAMETER_DATASET = "2026_Khudaida"
SPECIES = ["H2O", "Ethanol", "Butanol", "Na+", "Cl-"]
FORMULA_SPECIES = ["H2O", "Ethanol", "Isobutanol", "NaCl"]
IDX5 = {name: i for i, name in enumerate(SPECIES)}
IDX4 = {name: i for i, name in enumerate(FORMULA_SPECIES)}
SQRT3_OVER_2 = math.sqrt(3.0) / 2.0

MW_ION = np.asarray([18.01528e-3, 46.068e-3, 74.1216e-3, 22.98e-3, 35.45e-3], dtype=float)
MW_FORMULA = np.asarray([18.01528e-3, 46.068e-3, 74.1216e-3, 58.43e-3], dtype=float)

BLACK = "#000000"
RED = "#d62728"
BLUE = "#1f4ed8"
LIGHT_GREEN = "#63c46b"
WATER_COLOR = "#2f6fb3"
ETHANOL_COLOR = "#c25a14"
ISOBUTANOL_COLOR = "#5f8f1f"

ORGANIC_SIDE_SCALED_RANGE_5WT = {
    "water_min": 0.30,
    "water_max": 0.50,
    "ethanol_min": 0.00,
    "ethanol_max": 0.20,
    "isobutanol_min": 0.40,
    "isobutanol_max": 0.60,
}
ORGANIC_SIDE_SCALED_RANGE_10WT = {
    "water_min": 0.20,
    "water_max": 0.40,
    "ethanol_min": 0.00,
    "ethanol_max": 0.20,
    "isobutanol_min": 0.55,
    "isobutanol_max": 0.75,
}
SCALED_FIGURE_OVERRIDES: dict[int, dict] = {
    2: dict(ORGANIC_SIDE_SCALED_RANGE_5WT),
    3: dict(ORGANIC_SIDE_SCALED_RANGE_5WT),
    4: dict(ORGANIC_SIDE_SCALED_RANGE_5WT),
    5: dict(ORGANIC_SIDE_SCALED_RANGE_10WT),
    6: dict(ORGANIC_SIDE_SCALED_RANGE_10WT),
    7: dict(ORGANIC_SIDE_SCALED_RANGE_10WT),
}


def configure_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Serif",
            "font.size": 10,
            "axes.linewidth": 1.0,
            "axes.grid": False,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": False,
            "ytick.right": False,
            "legend.frameon": False,
            "mathtext.default": "regular",
        }
    )


def save_figure(fig: plt.Figure, path: Path) -> None:
    path = paper_validation_output_path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    save_plot_figure(fig, path, dpi=300, svg_companion=True)


def add_figure_caption(
    fig: plt.Figure, caption: str, *, left: float = 0.09, y: float = 0.02, fontsize: float = 9.0
) -> None:
    fig.text(left, y, caption, ha="left", va="bottom", fontsize=fontsize, wrap=True)


def write_csv_rows(path: Path, fieldnames: Iterable[str], rows: Iterable[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(list(rows), columns=list(fieldnames)).to_csv(path, index=False)


def out_path(fig_dir: Path, filename: str) -> Path:
    plot_set_dir = analysis_plot_set_dir(fig_dir)
    return plot_set_dir / "data" / filename


def figure_root(figure_number: int) -> Path:
    return ANALYSIS_ROOT / "figures" / f"figure_{figure_number}"


NO_SALT_293_DIGITIZED = [
    {"x_ethanol_aq": 0.0046, "distribution": 4.8, "separation": 9.2, "x_isobutanol_aq": 0.0020},
    {"x_ethanol_aq": 0.0069, "distribution": 4.4, "separation": 8.3, "x_isobutanol_aq": 0.0022},
    {"x_ethanol_aq": 0.0096, "distribution": 4.2, "separation": 7.7, "x_isobutanol_aq": 0.0025},
    {"x_ethanol_aq": 0.0148, "distribution": 4.1, "separation": 7.0, "x_isobutanol_aq": 0.0030},
    {"x_ethanol_aq": 0.0205, "distribution": 3.7, "separation": 5.6, "x_isobutanol_aq": 0.0036},
    {"x_ethanol_aq": 0.0243, "distribution": 3.6, "separation": 4.8, "x_isobutanol_aq": 0.0040},
    {"x_ethanol_aq": 0.0265, "distribution": 3.4, "separation": 4.5, "x_isobutanol_aq": 0.0044},
    {"x_ethanol_aq": 0.0320, "distribution": 3.0, "separation": 4.2, "x_isobutanol_aq": 0.0050},
]

NOMINAL_ETHANOL_FEED_WT = {
    0.05: (0.12, 0.10, 0.08, 0.06, 0.05, 0.04, 0.03, 0.02),
    0.10: (0.06, 0.05, 0.04, 0.03, 0.02, 0.01),
}

EePCSAFT_AAD_REFERENCE = {
    0.05: {
        293.15: {
            "grand": 0.0161,
            "organic": (0.0287, 0.0031, 0.0313, 0.0009),
            "aqueous": (0.0482, 0.0064, 0.0021, 0.0083),
        },
        303.15: {
            "grand": 0.0168,
            "organic": (0.0305, 0.0027, 0.0336, 0.0005),
            "aqueous": (0.0495, 0.0061, 0.0026, 0.0088),
        },
        313.15: {
            "grand": 0.0203,
            "organic": (0.0245, 0.0033, 0.0276, 0.0006),
            "aqueous": (0.0487, 0.0052, 0.0025, 0.0081),
        },
    },
    0.10: {
        293.15: {
            "grand": 0.0342,
            "organic": (0.0729, 0.0096, 0.0739, 0.0005),
            "aqueous": (0.0975, 0.0020, 0.0007, 0.0161),
        },
        303.15: {
            "grand": 0.0341,
            "organic": (0.0760, 0.0021, 0.0774, 0.0005),
            "aqueous": (0.0973, 0.0023, 0.0008, 0.0165),
        },
        313.15: {
            "grand": 0.0357,
            "organic": (0.0801, 0.0039, 0.0783, 0.0002),
            "aqueous": (0.1021, 0.0025, 0.0003, 0.0179),
        },
    },
}

ENRTL_AAD_REFERENCE = {
    0.05: {
        293.15: {
            "grand": 0.0046,
            "organic": (0.0059, 0.0048, 0.0010, 0.0011),
            "aqueous": (0.0102, 0.0025, 0.0080, 0.0029),
        },
        303.15: {
            "grand": 0.0047,
            "organic": (0.0053, 0.0040, 0.0015, 0.0007),
            "aqueous": (0.0120, 0.0029, 0.0097, 0.0012),
        },
        313.15: {
            "grand": 0.0078,
            "organic": (0.0147, 0.0117, 0.0041, 0.0011),
            "aqueous": (0.0128, 0.0030, 0.0126, 0.0027),
        },
    },
    0.10: {
        293.15: {
            "grand": 0.0086,
            "organic": (0.0010, 0.0026, 0.0022, 0.0007),
            "aqueous": (0.0312, 0.0008, 0.0062, 0.0241),
        },
        303.15: {
            "grand": 0.0083,
            "organic": (0.0011, 0.0022, 0.0033, 0.0022),
            "aqueous": (0.0289, 0.0006, 0.0079, 0.0204),
        },
        313.15: {
            "grand": 0.0091,
            "organic": (0.0011, 0.0035, 0.0028, 0.0004),
            "aqueous": (0.0323, 0.0020, 0.0098, 0.0205),
        },
    },
}

SUPPORTING_FIGURE_PANELS = {
    2: [
        ("a", 2, 293.15, 0.05),
        ("b", 3, 303.15, 0.05),
        ("c", 4, 313.15, 0.05),
    ],
    3: [
        ("a", 5, 293.15, 0.10),
        ("b", 6, 303.15, 0.10),
        ("c", 7, 313.15, 0.10),
    ],
}


def source_input_path(filename: str) -> Path:
    return SOURCE_INPUT_ROOT / filename


def source_experimental_cases() -> dict[tuple[float, float], list[tuple[int, tuple[float, ...], tuple[float, ...]]]]:
    path = source_input_path("table_3_4_experimental_tielines.csv")
    cases: dict[tuple[float, float], list[tuple[int, tuple[float, ...], tuple[float, ...]]]] = {}
    with path.open("r", newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        for raw in reader:
            salt_wt = float(raw["salt_wtfrac"])
            temperature_k = float(raw["temperature_K"])
            key = (salt_wt, temperature_k)
            organic = tuple(
                float(raw[name]) for name in ("x_water_org", "x_ethanol_org", "x_isobutanol_org", "x_nacl_org")
            )
            aqueous = tuple(
                float(raw[name]) for name in ("x_water_aq", "x_ethanol_aq", "x_isobutanol_aq", "x_nacl_aq")
            )
            cases.setdefault(key, []).append((int(raw["tie_line"]), organic, aqueous))
    return cases


def _experimental_rows(salt_wt: float, temperature_k: float) -> list[dict]:
    rows = []
    for tie_line, organic, aqueous in source_experimental_cases()[(salt_wt, temperature_k)]:
        rows.append(
            {
                "tie_line": tie_line,
                "temperature_K": float(temperature_k),
                "salt_wtfrac": float(salt_wt),
                "organic_formula": np.asarray(organic, dtype=float),
                "aqueous_formula": np.asarray(aqueous, dtype=float),
            }
        )
    return rows


def _digitized_feed_rows_for_figure(figure_number: int, temperature_k: float, salt_wt: float) -> list[dict] | None:
    source_path = analysis_data_path(figure_root(figure_number), "feed_compositions_digitized.csv", kind="input")
    if not source_path.exists():
        return None
    rows = []
    with source_path.open("r", newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        for idx, raw in enumerate(reader, start=1):
            ethanol = float(raw["x_ethanol_salt_free"])
            isobutanol = float(raw["x_isobutanol_salt_free"])
            water = float(raw["x_water_salt_free"])
            feed_formula = build_feed_formula_from_salt_free_molefractions(
                np.asarray([water, ethanol, isobutanol], dtype=float), salt_wt
            )
            rows.append(
                {
                    "tie_line": idx,
                    "temperature_K": float(temperature_k),
                    "salt_wtfrac": float(salt_wt),
                    "feed_formula": feed_formula,
                    "ethanol_total_feed_wtfrac": "",
                    "source": raw.get("source", "digitized_user_supplied") or "digitized_user_supplied",
                }
            )
    return rows


def formula_to_ion_basis(x_formula: np.ndarray) -> np.ndarray:
    x_formula = np.asarray(x_formula, dtype=float)
    denom = 1.0 + float(x_formula[3])
    return np.asarray(
        [
            x_formula[0] / denom,
            x_formula[1] / denom,
            x_formula[2] / denom,
            x_formula[3] / denom,
            x_formula[3] / denom,
        ],
        dtype=float,
    )


def ion_to_formula_basis(x_ion: np.ndarray) -> np.ndarray:
    x_ion = np.asarray(x_ion, dtype=float)
    salt_formula = 0.5 * float(x_ion[IDX5["Na+"]] + x_ion[IDX5["Cl-"]])
    denom = float(x_ion[IDX5["H2O"]] + x_ion[IDX5["Ethanol"]] + x_ion[IDX5["Butanol"]] + salt_formula)
    return np.asarray(
        [
            x_ion[IDX5["H2O"]] / denom,
            x_ion[IDX5["Ethanol"]] / denom,
            x_ion[IDX5["Butanol"]] / denom,
            salt_formula / denom,
        ],
        dtype=float,
    )


def salt_free_from_formula(x_formula: np.ndarray) -> np.ndarray:
    x_formula = np.asarray(x_formula, dtype=float)
    total = float(np.sum(x_formula[:3]))
    return np.asarray(x_formula[:3], dtype=float) / total


def ternary_xy_from_formula(x_formula: np.ndarray) -> tuple[float, float]:
    salt_free = salt_free_from_formula(x_formula)
    x = float(salt_free[2] + 0.5 * salt_free[1])
    y = float(SQRT3_OVER_2 * salt_free[1])
    return x, y


def _scaled_xy_from_formula(
    x_formula: np.ndarray,
    *,
    ethanol_max: float,
    isobutanol_max: float,
) -> tuple[float, float]:
    salt_free = salt_free_from_formula(x_formula)
    ethanol = float(salt_free[1]) / ethanol_max
    isobutanol = float(salt_free[2]) / isobutanol_max
    x = isobutanol + 0.5 * ethanol
    y = SQRT3_OVER_2 * ethanol
    return x, y


def _display_scaled_xy_from_formula(
    x_formula: np.ndarray,
    *,
    water_min: float,
    water_max: float,
    ethanol_min: float,
    ethanol_max: float,
    isobutanol_min: float,
    isobutanol_max: float,
) -> tuple[float, float]:
    salt_free = salt_free_from_formula(x_formula)
    water = float(salt_free[0])
    ethanol = float(salt_free[1])
    isobutanol = float(salt_free[2])

    def _scale(value: float, lower: float, upper: float) -> float:
        if abs(upper - lower) < 1.0e-12:
            return 0.5
        return float(np.clip((value - lower) / (upper - lower), 0.0, 1.0))

    scaled = np.asarray(
        [
            _scale(water, water_min, water_max),
            _scale(ethanol, ethanol_min, ethanol_max),
            _scale(isobutanol, isobutanol_min, isobutanol_max),
        ],
        dtype=float,
    )
    scaled = scaled / np.sum(scaled)
    x = float(scaled[2] + 0.5 * scaled[1])
    y = float(SQRT3_OVER_2 * scaled[1])
    return x, y


def _scaled_triangle_vertices(axis_scale: float) -> np.ndarray:
    return np.asarray(
        [
            [0.0, 0.0],
            [axis_scale, 0.0],
            [0.5 * axis_scale, SQRT3_OVER_2 * axis_scale],
            [0.0, 0.0],
        ],
        dtype=float,
    )


def _draw_ternary_axes(ax: plt.Axes) -> None:
    triangle = np.asarray([[0.0, 0.0], [1.0, 0.0], [0.5, SQRT3_OVER_2], [0.0, 0.0]], dtype=float)
    ax.plot(triangle[:, 0], triangle[:, 1], color="black", linewidth=1.2, zorder=2)
    left_normal = np.asarray([-SQRT3_OVER_2, 0.5], dtype=float)
    right_normal = np.asarray([SQRT3_OVER_2, 0.5], dtype=float)
    for frac in np.linspace(0.05, 0.95, 19):
        is_major = abs(frac * 10.0 - round(frac * 10.0)) < 1.0e-9
        grid_linewidth = 0.6 if is_major else 0.35
        grid_alpha = 0.55 if is_major else 0.28
        # Constant ethanol: horizontal lines between left and right edges.
        ax.plot(
            [0.5 * frac, 1.0 - 0.5 * frac],
            [SQRT3_OVER_2 * frac, SQRT3_OVER_2 * frac],
            color=ETHANOL_COLOR,
            linewidth=grid_linewidth,
            alpha=grid_alpha,
            linestyle="--",
            zorder=0,
        )
        # Constant isobutanol: lines parallel to the left edge.
        ax.plot(
            [frac, 0.5 + 0.5 * frac],
            [0.0, SQRT3_OVER_2 * (1.0 - frac)],
            color=ISOBUTANOL_COLOR,
            linewidth=grid_linewidth,
            alpha=grid_alpha,
            linestyle="--",
            zorder=0,
        )
        # Constant water: lines parallel to the right edge.
        ax.plot(
            [1.0 - frac, 0.5 * (1.0 - frac)],
            [0.0, SQRT3_OVER_2 * (1.0 - frac)],
            color=WATER_COLOR,
            linewidth=grid_linewidth,
            alpha=grid_alpha,
            linestyle="--",
            zorder=0,
        )
    ax.text(0.50, -0.10, r"$x$ Isobutanol", color=ISOBUTANOL_COLOR, ha="center", va="top", fontsize=11)
    ax.text(0.15, 0.52, r"$x$ Water", color=WATER_COLOR, rotation=60, ha="center", va="center", fontsize=11)
    ax.text(0.85, 0.50, r"$x$ Ethanol", color=ETHANOL_COLOR, rotation=-60, ha="center", va="center", fontsize=11)
    left_tick_offset = 0.035
    right_tick_offset = 0.035
    for frac in np.arange(0.0, 1.01, 0.1):
        if frac < 1.0:
            ax.text(frac, -0.03, f"{frac:.1f}", color=ISOBUTANOL_COLOR, ha="center", va="top", fontsize=8)
        if 0.0 < frac < 1.0:
            # Left edge shows water decreasing from 1 at the lower-left corner to 0 at the apex.
            left_point = np.asarray([0.5 * frac, SQRT3_OVER_2 * frac], dtype=float) + left_tick_offset * left_normal
            ax.text(
                left_point[0],
                left_point[1],
                f"{1.0 - frac:.1f}",
                color=WATER_COLOR,
                ha="right",
                va="center",
                fontsize=8,
            )
            # Right edge shows ethanol increasing from 0 at the lower-right corner to 1 at the apex.
            right_point = (
                np.asarray([1.0 - 0.5 * frac, SQRT3_OVER_2 * frac], dtype=float) + right_tick_offset * right_normal
            )
            ax.text(
                right_point[0], right_point[1], f"{frac:.1f}", color=ETHANOL_COLOR, ha="left", va="center", fontsize=8
            )
    left_bottom = np.asarray([0.0, 0.0], dtype=float) + left_tick_offset * left_normal
    left_top = np.asarray([0.5, SQRT3_OVER_2], dtype=float) + left_tick_offset * left_normal
    right_bottom = np.asarray([1.0, 0.0], dtype=float) + right_tick_offset * right_normal
    right_top = np.asarray([0.5, SQRT3_OVER_2], dtype=float) + right_tick_offset * right_normal
    ax.text(left_bottom[0], left_bottom[1], "1.0", color=WATER_COLOR, ha="right", va="center", fontsize=8)
    ax.text(left_top[0], left_top[1], "0.0", color=WATER_COLOR, ha="right", va="center", fontsize=8)
    ax.text(right_bottom[0], right_bottom[1], "0.0", color=ETHANOL_COLOR, ha="left", va="center", fontsize=8)
    ax.text(right_top[0], right_top[1], "1.0", color=ETHANOL_COLOR, ha="left", va="center", fontsize=8)
    ax.set_xlim(-0.08, 1.08)
    ax.set_ylim(-0.08, SQRT3_OVER_2 + 0.06)
    ax.set_aspect("equal")
    ax.axis("off")


def _draw_scaled_ternary_axes(
    ax: plt.Axes,
    *,
    water_min: float,
    ethanol_max: float,
    isobutanol_max: float,
    tick_format: str = ".1f",
) -> None:
    triangle = _scaled_triangle_vertices(1.0)
    ax.plot(triangle[:, 0], triangle[:, 1], color="black", linewidth=1.2, zorder=2)
    left_normal = np.asarray([-SQRT3_OVER_2, 0.5], dtype=float)
    right_normal = np.asarray([SQRT3_OVER_2, 0.5], dtype=float)
    for frac in np.linspace(0.05, 0.95, 19):
        is_major = abs(frac * 10.0 - round(frac * 10.0)) < 1.0e-9
        grid_linewidth = 0.6 if is_major else 0.35
        grid_alpha = 0.55 if is_major else 0.28
        ethanol = frac
        ax.plot(
            [0.5 * ethanol, 1.0 - 0.5 * ethanol],
            [SQRT3_OVER_2 * ethanol, SQRT3_OVER_2 * ethanol],
            color=ETHANOL_COLOR,
            linewidth=grid_linewidth,
            alpha=grid_alpha,
            linestyle="--",
            zorder=0,
        )
        isobutanol = frac
        ax.plot(
            [isobutanol, 0.5 + 0.5 * isobutanol],
            [0.0, SQRT3_OVER_2 * (1.0 - isobutanol)],
            color=ISOBUTANOL_COLOR,
            linewidth=grid_linewidth,
            alpha=grid_alpha,
            linestyle="--",
            zorder=0,
        )
        water_drop = frac
        ax.plot(
            [1.0 - water_drop, 0.5 * (1.0 - water_drop)],
            [0.0, SQRT3_OVER_2 * (1.0 - water_drop)],
            color=WATER_COLOR,
            linewidth=grid_linewidth,
            alpha=grid_alpha,
            linestyle="--",
            zorder=0,
        )
    ax.text(0.50, -0.085, r"$x$ Isobutanol", color=ISOBUTANOL_COLOR, ha="center", va="top", fontsize=10)
    ax.text(-0.08, 0.52, r"$x$ Water", color=WATER_COLOR, rotation=60, ha="center", va="center", fontsize=10)
    ax.text(1.08, 0.52, r"$x$ Ethanol", color=ETHANOL_COLOR, rotation=-60, ha="center", va="center", fontsize=10)
    left_tick_offset = 0.018
    right_tick_offset = 0.02
    for frac in np.arange(0.0, 1.01, 0.1):
        isobutanol = isobutanol_max * frac
        ethanol = ethanol_max * frac
        water = 1.0 - (1.0 - water_min) * frac
        if frac < 1.0:
            ax.text(
                frac, -0.03, format(isobutanol, tick_format), color=ISOBUTANOL_COLOR, ha="center", va="top", fontsize=7
            )
        if 0.0 < frac < 1.0:
            left_point = np.asarray([0.5 * frac, SQRT3_OVER_2 * frac], dtype=float) + left_tick_offset * left_normal
            ax.text(
                left_point[0],
                left_point[1],
                format(water, tick_format),
                color=WATER_COLOR,
                ha="right",
                va="center",
                fontsize=7,
            )
            right_point = (
                np.asarray([1.0 - 0.5 * frac, SQRT3_OVER_2 * frac], dtype=float) + right_tick_offset * right_normal
            )
            ax.text(
                right_point[0],
                right_point[1],
                format(ethanol, tick_format),
                color=ETHANOL_COLOR,
                ha="left",
                va="center",
                fontsize=7,
            )
    left_bottom = np.asarray([0.0, 0.0], dtype=float) + left_tick_offset * left_normal
    left_top = np.asarray([0.5, SQRT3_OVER_2], dtype=float) + left_tick_offset * left_normal
    right_bottom = np.asarray([1.0, 0.0], dtype=float) + right_tick_offset * right_normal
    right_top = np.asarray([0.5, SQRT3_OVER_2], dtype=float) + right_tick_offset * right_normal
    ax.text(
        left_bottom[0], left_bottom[1], format(1.0, tick_format), color=WATER_COLOR, ha="right", va="center", fontsize=7
    )
    ax.text(
        left_top[0], left_top[1], format(water_min, tick_format), color=WATER_COLOR, ha="right", va="center", fontsize=7
    )
    ax.text(
        right_bottom[0],
        right_bottom[1],
        format(0.0, tick_format),
        color=ETHANOL_COLOR,
        ha="left",
        va="center",
        fontsize=7,
    )
    ax.text(
        right_top[0],
        right_top[1],
        format(ethanol_max, tick_format),
        color=ETHANOL_COLOR,
        ha="left",
        va="center",
        fontsize=7,
    )
    ax.set_xlim(-0.12, 1.14)
    ax.set_ylim(-0.10, SQRT3_OVER_2 + 0.08)
    ax.set_aspect("equal")
    ax.axis("off")


def _draw_display_scaled_axes(
    ax: plt.Axes,
    *,
    water_min: float,
    water_max: float,
    ethanol_min: float,
    ethanol_max: float,
    isobutanol_min: float,
    isobutanol_max: float,
    tick_format: str = ".3f",
) -> None:
    triangle = _scaled_triangle_vertices(1.0)
    ax.plot(triangle[:, 0], triangle[:, 1], color="black", linewidth=1.2, zorder=2)
    left_normal = np.asarray([-SQRT3_OVER_2, 0.5], dtype=float)
    right_normal = np.asarray([SQRT3_OVER_2, 0.5], dtype=float)
    for frac in np.linspace(0.05, 0.95, 19):
        is_major = abs(frac * 10.0 - round(frac * 10.0)) < 1.0e-9
        grid_linewidth = 0.6 if is_major else 0.35
        grid_alpha = 0.55 if is_major else 0.28
        ax.plot(
            [0.5 * frac, 1.0 - 0.5 * frac],
            [SQRT3_OVER_2 * frac, SQRT3_OVER_2 * frac],
            color=ETHANOL_COLOR,
            linewidth=grid_linewidth,
            alpha=grid_alpha,
            linestyle="--",
            zorder=0,
        )
        ax.plot(
            [frac, 0.5 + 0.5 * frac],
            [0.0, SQRT3_OVER_2 * (1.0 - frac)],
            color=ISOBUTANOL_COLOR,
            linewidth=grid_linewidth,
            alpha=grid_alpha,
            linestyle="--",
            zorder=0,
        )
        ax.plot(
            [1.0 - frac, 0.5 * (1.0 - frac)],
            [0.0, SQRT3_OVER_2 * (1.0 - frac)],
            color=WATER_COLOR,
            linewidth=grid_linewidth,
            alpha=grid_alpha,
            linestyle="--",
            zorder=0,
        )
    ax.text(0.50, -0.085, r"$x$ Isobutanol", color=ISOBUTANOL_COLOR, ha="center", va="top", fontsize=10)
    ax.text(-0.08, 0.52, r"$x$ Water", color=WATER_COLOR, rotation=60, ha="center", va="center", fontsize=10)
    ax.text(1.08, 0.52, r"$x$ Ethanol", color=ETHANOL_COLOR, rotation=-60, ha="center", va="center", fontsize=10)
    left_tick_offset = 0.018
    right_tick_offset = 0.02
    for frac in np.arange(0.0, 1.01, 0.2):
        water = water_max + frac * (water_min - water_max)
        ethanol = ethanol_min + frac * (ethanol_max - ethanol_min)
        isobutanol = isobutanol_min + frac * (isobutanol_max - isobutanol_min)
        if frac < 1.0:
            ax.text(
                frac, -0.03, format(isobutanol, tick_format), color=ISOBUTANOL_COLOR, ha="center", va="top", fontsize=7
            )
        if 0.0 < frac < 1.0:
            left_point = np.asarray([0.5 * frac, SQRT3_OVER_2 * frac], dtype=float) + left_tick_offset * left_normal
            ax.text(
                left_point[0],
                left_point[1],
                format(water, tick_format),
                color=WATER_COLOR,
                ha="right",
                va="center",
                fontsize=7,
            )
            right_point = (
                np.asarray([1.0 - 0.5 * frac, SQRT3_OVER_2 * frac], dtype=float) + right_tick_offset * right_normal
            )
            ax.text(
                right_point[0],
                right_point[1],
                format(ethanol, tick_format),
                color=ETHANOL_COLOR,
                ha="left",
                va="center",
                fontsize=7,
            )
    left_bottom = np.asarray([0.0, 0.0], dtype=float) + left_tick_offset * left_normal
    left_top = np.asarray([0.5, SQRT3_OVER_2], dtype=float) + left_tick_offset * left_normal
    right_bottom = np.asarray([1.0, 0.0], dtype=float) + right_tick_offset * right_normal
    right_top = np.asarray([0.5, SQRT3_OVER_2], dtype=float) + right_tick_offset * right_normal
    ax.text(
        left_bottom[0],
        left_bottom[1],
        format(water_max, tick_format),
        color=WATER_COLOR,
        ha="right",
        va="center",
        fontsize=7,
    )
    ax.text(
        left_top[0], left_top[1], format(water_min, tick_format), color=WATER_COLOR, ha="right", va="center", fontsize=7
    )
    ax.text(
        right_bottom[0],
        right_bottom[1],
        format(ethanol_min, tick_format),
        color=ETHANOL_COLOR,
        ha="left",
        va="center",
        fontsize=7,
    )
    ax.text(
        right_top[0],
        right_top[1],
        format(ethanol_max, tick_format),
        color=ETHANOL_COLOR,
        ha="left",
        va="center",
        fontsize=7,
    )
    ax.set_xlim(-0.12, 1.14)
    ax.set_ylim(-0.10, SQRT3_OVER_2 + 0.08)
    ax.set_aspect("equal")
    ax.axis("off")


def _scaled_axis_scale_from_rows(exp_rows: list[dict], model_rows: list[dict], feed_rows: list[dict]) -> float:
    ethanol_values = []
    for row in exp_rows:
        ethanol_values.extend(
            [
                float(salt_free_from_formula(row["organic_formula"])[1]),
                float(salt_free_from_formula(row["aqueous_formula"])[1]),
            ]
        )
    for row in model_rows:
        if np.all(np.isfinite(row["organic_formula"])) and np.all(np.isfinite(row["aqueous_formula"])):
            ethanol_values.extend(
                [
                    float(salt_free_from_formula(row["organic_formula"])[1]),
                    float(salt_free_from_formula(row["aqueous_formula"])[1]),
                ]
            )
    for row in feed_rows:
        ethanol_values.append(float(salt_free_from_formula(row["feed_formula"])[1]))
    max_ethanol = max(ethanol_values) if ethanol_values else 0.2
    return max(0.1, min(1.0, math.ceil(max_ethanol * 10.0) / 10.0))


def _scaled_water_min_from_rows(exp_rows: list[dict], model_rows: list[dict], feed_rows: list[dict]) -> float:
    water_values = []
    for row in exp_rows:
        water_values.extend(
            [
                float(salt_free_from_formula(row["organic_formula"])[0]),
                float(salt_free_from_formula(row["aqueous_formula"])[0]),
            ]
        )
    for row in model_rows:
        if np.all(np.isfinite(row["organic_formula"])) and np.all(np.isfinite(row["aqueous_formula"])):
            water_values.extend(
                [
                    float(salt_free_from_formula(row["organic_formula"])[0]),
                    float(salt_free_from_formula(row["aqueous_formula"])[0]),
                ]
            )
    for row in feed_rows:
        water_values.append(float(salt_free_from_formula(row["feed_formula"])[0]))
    min_water = min(water_values) if water_values else 0.0
    return max(0.0, min(0.9, math.floor(min_water * 10.0) / 10.0))


def _display_scaled_ranges_from_rows(
    exp_rows: list[dict], model_rows: list[dict], feed_rows: list[dict]
) -> dict[str, float]:
    water_values = []
    ethanol_values = []
    isobutanol_values = []

    def _append_from_formula(x_formula: np.ndarray) -> None:
        salt_free = salt_free_from_formula(x_formula)
        water_values.append(float(salt_free[0]))
        ethanol_values.append(float(salt_free[1]))
        isobutanol_values.append(float(salt_free[2]))

    for row in exp_rows:
        _append_from_formula(row["organic_formula"])
        _append_from_formula(row["aqueous_formula"])
    for row in model_rows:
        if np.all(np.isfinite(row["organic_formula"])) and np.all(np.isfinite(row["aqueous_formula"])):
            _append_from_formula(row["organic_formula"])
            _append_from_formula(row["aqueous_formula"])
    for row in feed_rows:
        _append_from_formula(row["feed_formula"])

    return {
        "water_min": min(water_values),
        "water_max": max(water_values),
        "ethanol_min": min(ethanol_values),
        "ethanol_max": max(ethanol_values),
        "isobutanol_min": min(isobutanol_values),
        "isobutanol_max": max(isobutanol_values),
    }


def _rounded_display_ranges(ranges: dict[str, float]) -> dict[str, float]:
    rounded = {
        "water_min": max(0.0, math.floor(ranges["water_min"] * 10.0) / 10.0),
        "water_max": min(1.0, math.ceil(ranges["water_max"] * 10.0) / 10.0),
        "ethanol_min": max(0.0, math.floor(ranges["ethanol_min"] * 10.0) / 10.0),
        "ethanol_max": min(1.0, math.ceil(ranges["ethanol_max"] * 10.0) / 10.0),
        "isobutanol_min": max(0.0, math.floor(ranges["isobutanol_min"] * 10.0) / 10.0),
        "isobutanol_max": min(1.0, math.ceil(ranges["isobutanol_max"] * 10.0) / 10.0),
    }
    for key_min, key_max in (
        ("water_min", "water_max"),
        ("ethanol_min", "ethanol_max"),
        ("isobutanol_min", "isobutanol_max"),
    ):
        if rounded[key_max] <= rounded[key_min]:
            rounded[key_max] = min(1.0, rounded[key_min] + 0.1)
    return rounded


def _plot_tie_lines(
    ax: plt.Axes,
    rows: list[dict],
    color: str,
    marker: str,
    label: str,
    *,
    linewidth: float = 1.0,
    markersize: float = 18.0,
    linestyle: str = "-",
    xy_transform=ternary_xy_from_formula,
) -> None:
    if not rows:
        return
    for idx, row in enumerate(rows):
        aq_xy = xy_transform(row["aqueous_formula"])
        org_xy = xy_transform(row["organic_formula"])
        ax.plot(
            [aq_xy[0], org_xy[0]],
            [aq_xy[1], org_xy[1]],
            color=color,
            linewidth=linewidth,
            linestyle=linestyle,
            zorder=3,
            label=label if idx == 0 else None,
        )
        ax.scatter(
            [aq_xy[0], org_xy[0]],
            [aq_xy[1], org_xy[1]],
            s=markersize,
            marker=marker,
            facecolors=color,
            edgecolors=color,
            linewidths=0.5,
            zorder=4,
        )


def _model_objective(exp_row: dict, pred_org: np.ndarray, pred_aq: np.ndarray) -> float:
    return float(
        np.sum(np.abs(pred_org - exp_row["organic_formula"])) + np.sum(np.abs(pred_aq - exp_row["aqueous_formula"]))
    )


def build_feed_formula_from_total_feed_weights(ethanol_wt: float, salt_wt: float) -> np.ndarray:
    water_wt = 1.0 - 0.50 - salt_wt - ethanol_wt
    weights = np.asarray([water_wt, ethanol_wt, 0.50, salt_wt], dtype=float)
    moles = weights / MW_FORMULA
    return moles / np.sum(moles)


def build_feed_formula_from_salt_free_molefractions(salt_free_xyz: np.ndarray, salt_wt: float) -> np.ndarray:
    salt_free_xyz = np.asarray(salt_free_xyz, dtype=float)
    salt_free_xyz = salt_free_xyz / np.sum(salt_free_xyz)
    neutral_average_mw = float(np.dot(salt_free_xyz, MW_FORMULA[:3]))
    salt_moles_per_neutral_mole = salt_wt * neutral_average_mw / ((1.0 - salt_wt) * MW_FORMULA[3])
    total_moles = 1.0 + salt_moles_per_neutral_mole
    return np.asarray(
        [
            salt_free_xyz[0] / total_moles,
            salt_free_xyz[1] / total_moles,
            salt_free_xyz[2] / total_moles,
            salt_moles_per_neutral_mole / total_moles,
        ],
        dtype=float,
    )


def _derived_feed_rows(salt_wt: float, temperature_k: float) -> list[dict]:
    nominal_ethanol = NOMINAL_ETHANOL_FEED_WT[salt_wt]
    return [
        {
            "tie_line": idx,
            "temperature_K": float(temperature_k),
            "salt_wtfrac": float(salt_wt),
            "feed_formula": build_feed_formula_from_total_feed_weights(ethanol_wt, salt_wt),
            "ethanol_total_feed_wtfrac": float(ethanol_wt),
            "source": "derived_from_methods_1to1_feed_rule",
        }
        for idx, ethanol_wt in enumerate(nominal_ethanol, start=1)
    ]


def _plot_feed_points(
    ax: plt.Axes,
    rows: list[dict],
    color: str = LIGHT_GREEN,
    marker: str = "^",
    markersize: float = 24.0,
    label: str | None = None,
    xy_transform=ternary_xy_from_formula,
    display_ranges: dict[str, float] | None = None,
) -> None:
    if not rows:
        return
    xs = []
    ys = []
    for row in rows:
        if display_ranges is not None and not _formula_in_display_scaled_ranges(row["feed_formula"], display_ranges):
            continue
        x_coord, y_coord = xy_transform(row["feed_formula"])
        xs.append(x_coord)
        ys.append(y_coord)
    if xs:
        ax.scatter(
            xs,
            ys,
            s=markersize,
            marker=marker,
            facecolors=color,
            edgecolors=color,
            linewidths=0.5,
            zorder=5,
            label=label,
        )


def _paper_epcsaft_digitized_rows(fig_dir: Path, temperature_k: float, salt_wt: float) -> list[dict]:
    path = analysis_data_path(fig_dir, "paper_epcsaft_digitized.csv", kind="input")
    if not path.exists():
        return []
    rows = []
    with path.open("r", newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        for raw in reader:
            salt_free = np.asarray(
                [
                    float(raw["x_water_salt_free"]),
                    float(raw["x_ethanol_salt_free"]),
                    float(raw["x_isobutanol_salt_free"]),
                ],
                dtype=float,
            )
            salt_free = salt_free / float(np.sum(salt_free))
            formula = np.asarray([salt_free[0], salt_free[1], salt_free[2], 0.0], dtype=float)
            rows.append(
                {
                    "tie_line": int(raw["point_id"]),
                    "temperature_K": float(temperature_k),
                    "salt_wtfrac": float(salt_wt),
                    "organic_formula": formula,
                    "source": raw.get("source", "digitized_user_supplied") or "digitized_user_supplied",
                }
            )
    return rows


def _plot_formula_series(
    ax: plt.Axes,
    rows: list[dict],
    key: str,
    color: str,
    marker: str,
    label: str,
    *,
    linewidth: float = 1.0,
    markersize: float = 18.0,
    linestyle: str = "-",
    xy_transform=ternary_xy_from_formula,
    display_ranges: dict[str, float] | None = None,
) -> None:
    xs = []
    ys = []
    for row in rows:
        x_formula = row[key]
        if not np.all(np.isfinite(x_formula)):
            continue
        if display_ranges is not None and not _formula_in_display_scaled_ranges(x_formula, display_ranges):
            continue
        x_coord, y_coord = xy_transform(x_formula)
        xs.append(x_coord)
        ys.append(y_coord)
    if not xs:
        return
    if len(xs) > 1 and linestyle:
        ax.plot(xs, ys, color=color, linewidth=linewidth, linestyle=linestyle, zorder=3, label=label)
        scatter_label = None
    else:
        scatter_label = label
    ax.scatter(
        xs,
        ys,
        s=markersize,
        marker=marker,
        facecolors=color,
        edgecolors=color,
        linewidths=0.5,
        zorder=4,
        label=scatter_label,
    )


def _scaled_ranges_for_figure(
    figure_number: int,
    exp_rows: list[dict],
    model_rows: list[dict],
    feed_rows: list[dict],
) -> dict[str, float]:
    override = SCALED_FIGURE_OVERRIDES.get(figure_number)
    if override is not None:
        return dict(override)
    return _rounded_display_ranges(_display_scaled_ranges_from_rows(exp_rows, model_rows, feed_rows))


def _formula_in_display_scaled_ranges(
    x_formula: np.ndarray, ranges: dict[str, float], *, tolerance: float = 0.005
) -> bool:
    salt_free = salt_free_from_formula(x_formula)
    checks = (
        ("water", float(salt_free[0])),
        ("ethanol", float(salt_free[1])),
        ("isobutanol", float(salt_free[2])),
    )
    return all(
        ranges[f"{name}_min"] - tolerance <= value <= ranges[f"{name}_max"] + tolerance for name, value in checks
    )


def _plot_phase_points(
    ax: plt.Axes,
    rows: list[dict],
    phase: str,
    color: str,
    marker: str,
    *,
    markersize: float = 22.0,
    label: str | None = None,
    xy_transform=ternary_xy_from_formula,
    display_ranges: dict[str, float] | None = None,
) -> None:
    xs = []
    ys = []
    key = f"{phase}_formula"
    for row in rows:
        x_formula = row[key]
        if not np.all(np.isfinite(x_formula)):
            continue
        if display_ranges is not None and not _formula_in_display_scaled_ranges(x_formula, display_ranges):
            continue
        x_coord, y_coord = xy_transform(x_formula)
        xs.append(x_coord)
        ys.append(y_coord)
    if xs:
        ax.scatter(
            xs,
            ys,
            s=markersize,
            marker=marker,
            facecolors=color,
            edgecolors=color,
            linewidths=0.5,
            zorder=4,
            label=label,
        )


def _candidate_formula_feeds(exp_row: dict, target_feed_formula: np.ndarray | None = None) -> list[np.ndarray]:
    feeds = []
    if target_feed_formula is not None and np.all(np.isfinite(target_feed_formula)):
        feeds.append(np.asarray(target_feed_formula, dtype=float) / np.sum(target_feed_formula))
    midpoint = 0.5 * (exp_row["organic_formula"] + exp_row["aqueous_formula"])
    feeds.append(midpoint / np.sum(midpoint))
    unique = []
    for feed in feeds:
        if not any(np.allclose(feed, prior, atol=1e-10) for prior in unique):
            unique.append(feed)
    return unique


def _solve_formula_feed(temperature_k: float, feed_formula: np.ndarray) -> dict | None:
    z_feed = formula_to_ion_basis(feed_formula)
    if not np.all(np.isfinite(z_feed)):
        return None

    mixture = epcsaft.ePCSAFTMixture.from_dataset(PARAMETER_DATASET, SPECIES, z_feed, temperature_k)
    try:
        result = mixture.equilibrium(
            kind="electrolyte_lle",
            T=float(temperature_k),
            P=P_REF,
            z=z_feed,
            options=epcsaft.EquilibriumOptions(max_iterations=MODEL_SOLVE_MAX_ITERATIONS, tolerance=1.0e-8),
        )
    except epcsaft.SolutionError as exc:
        diagnostics = getattr(exc, "diagnostics", None) or (
            exc.args[1] if len(exc.args) > 1 and isinstance(exc.args[1], dict) else {}
        )
        return {
            "converged": False,
            "status": None,
            "message": str(exc.args[0] if exc.args else exc),
            "residual_norm": float(diagnostics.get("solver_residual_norm", np.nan)),
            "route_status": diagnostics.get("route_status"),
            "solver_status": diagnostics.get("solver_status"),
            "application_status": diagnostics.get("application_status"),
            "feed_formula": ion_to_formula_basis(z_feed),
            "organic_formula": np.full(4, np.nan),
            "aqueous_formula": np.full(4, np.nan),
            "beta_organic": np.nan,
            "beta_aqueous": np.nan,
            "split_norm": float(diagnostics.get("phase_distance", np.nan)),
            "objective": np.nan,
            "source": "epcsaft_native_v5",
        }

    phases = {phase.label: phase for phase in result.phases}
    if "org" not in phases or "aq" not in phases:
        return None
    org_formula = explicit_to_formula(phases["org"].composition)
    aq_formula = explicit_to_formula(phases["aq"].composition)
    residual_norm = float(result.diagnostics.get("solver_residual_norm", np.nan))
    return {
        "converged": True,
        "status": "accepted",
        "message": None,
        "residual_norm": residual_norm,
        "route_status": result.diagnostics.get("route_status"),
        "solver_status": result.diagnostics.get("solver_status"),
        "application_status": result.diagnostics.get("application_status"),
        "feed_formula": ion_to_formula_basis(z_feed),
        "organic_formula": org_formula,
        "aqueous_formula": aq_formula,
        "beta_organic": float(phases["org"].phase_fraction),
        "beta_aqueous": float(phases["aq"].phase_fraction),
        "split_norm": float(result.diagnostics.get("phase_distance", np.max(np.abs(org_formula - aq_formula)))),
        "objective": np.nan,
        "source": "epcsaft_native_v5",
    }


def solve_model_rows(exp_rows: list[dict], feed_rows: list[dict] | None = None) -> list[dict]:
    feed_map = {int(row["tie_line"]): np.asarray(row["feed_formula"], dtype=float) for row in (feed_rows or [])}
    solved = []
    for exp_row in exp_rows:
        best = None
        target_feed_formula = feed_map.get(int(exp_row["tie_line"]))
        for feed_idx, feed_formula in enumerate(_candidate_formula_feeds(exp_row, target_feed_formula)):
            candidate = _solve_formula_feed(exp_row["temperature_K"], feed_formula)
            if candidate is None:
                continue
            score = _model_objective(exp_row, candidate["organic_formula"], candidate["aqueous_formula"])
            if best is None or score < best["objective"]:
                best = {**candidate, "objective": score}
            if best is not None and target_feed_formula is not None and feed_idx == 0 and best["converged"]:
                break
        if best is None:
            best = {
                "converged": False,
                "status": None,
                "message": "No converged model candidate.",
                "residual_norm": np.nan,
                "route_status": None,
                "solver_status": None,
                "application_status": None,
                "feed_formula": np.full(4, np.nan),
                "organic_formula": np.full(4, np.nan),
                "aqueous_formula": np.full(4, np.nan),
                "beta_organic": np.nan,
                "beta_aqueous": np.nan,
                "split_norm": np.nan,
                "objective": np.nan,
            }
        solved.append(
            {
                "tie_line": exp_row["tie_line"],
                "temperature_K": exp_row["temperature_K"],
                "salt_wtfrac": exp_row["salt_wtfrac"],
                **best,
            }
        )
    return solved


def _parse_float(value: str | float | int | None) -> float:
    if value in (None, ""):
        return np.nan
    return float(value)


def _load_model_rows_from_csv(path: Path) -> list[dict]:
    grouped: dict[int, dict] = {}
    with path.open("r", newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        for raw in reader:
            tie_line = int(raw["tie_line"])
            entry = grouped.setdefault(
                tie_line,
                {
                    "tie_line": tie_line,
                    "temperature_K": _parse_float(raw.get("temperature_K")),
                    "salt_wtfrac": _parse_float(raw.get("salt_wtfrac")),
                    "organic_formula": np.full(4, np.nan),
                    "aqueous_formula": np.full(4, np.nan),
                    "feed_formula": np.full(4, np.nan),
                    "beta_organic": np.nan,
                    "beta_aqueous": np.nan,
                    "residual_norm": np.nan,
                    "objective": np.nan,
                    "converged": False,
                    "status": None,
                    "message": None,
                    "split_norm": np.nan,
                },
            )
            phase_vec = np.asarray(
                [
                    _parse_float(raw.get("x_water")),
                    _parse_float(raw.get("x_ethanol")),
                    _parse_float(raw.get("x_isobutanol")),
                    _parse_float(raw.get("x_nacl")),
                ],
                dtype=float,
            )
            phase_name = str(raw.get("phase", "")).strip().lower()
            if phase_name == "organic":
                entry["organic_formula"] = phase_vec
                entry["beta_organic"] = _parse_float(raw.get("beta"))
            elif phase_name == "aqueous":
                entry["aqueous_formula"] = phase_vec
                entry["beta_aqueous"] = _parse_float(raw.get("beta"))
            entry["residual_norm"] = _parse_float(raw.get("residual_norm"))
            entry["objective"] = _parse_float(raw.get("objective"))
            converged_raw = str(raw.get("converged", "")).strip().lower()
            entry["converged"] = converged_raw in {"true", "1", "yes"}
    return [grouped[key] for key in sorted(grouped)]


def _model_cache_path(fig_dir: Path) -> Path:
    return out_path(fig_dir, "model_tielines.csv")


def get_or_build_model_rows(
    fig_dir: Path, exp_rows: list[dict], feed_rows: list[dict] | None = None, force_recompute: bool = False
) -> list[dict]:
    cache_path = _model_cache_path(fig_dir)
    if not force_recompute and cache_path.exists():
        cached = _load_model_rows_from_csv(cache_path)
        cached_ties = {row["tie_line"] for row in cached}
        expected_ties = {row["tie_line"] for row in exp_rows}
        if cached_ties == expected_ties:
            return cached
    model_rows = solve_model_rows(exp_rows, feed_rows=feed_rows)
    write_case_data(fig_dir, exp_rows, model_rows)
    return model_rows


def _phase_rows_for_csv(rows: list[dict], source: str) -> list[dict]:
    out = []
    for row in rows:
        for phase_name, x in (("organic", row["organic_formula"]), ("aqueous", row["aqueous_formula"])):
            out.append(
                {
                    "tie_line": row["tie_line"],
                    "phase": phase_name,
                    "temperature_K": row["temperature_K"],
                    "salt_wtfrac": row["salt_wtfrac"],
                    "x_water": float(x[0]),
                    "x_ethanol": float(x[1]),
                    "x_isobutanol": float(x[2]),
                    "x_nacl": float(x[3]),
                    "source": source,
                }
            )
    return out


def _model_rows_for_csv(rows: list[dict]) -> list[dict]:
    out = []
    for row in rows:
        for phase_name, x, beta in (
            ("organic", row["organic_formula"], row["beta_organic"]),
            ("aqueous", row["aqueous_formula"], row["beta_aqueous"]),
        ):
            out.append(
                {
                    "tie_line": row["tie_line"],
                    "phase": phase_name,
                    "temperature_K": row["temperature_K"],
                    "salt_wtfrac": row["salt_wtfrac"],
                    "x_water": float(x[0]),
                    "x_ethanol": float(x[1]),
                    "x_isobutanol": float(x[2]),
                    "x_nacl": float(x[3]),
                    "beta": float(beta),
                    "residual_norm": float(row["residual_norm"]) if np.isfinite(row["residual_norm"]) else "",
                    "objective": float(row["objective"]) if np.isfinite(row["objective"]) else "",
                    "converged": bool(row["converged"]),
                    "source": str(row.get("source", "epcsaft_native_v5")),
                }
            )
    return out


def _feed_rows_for_csv(rows: list[dict]) -> list[dict]:
    out = []
    ordered_rows = sorted(rows, key=lambda row: float(salt_free_from_formula(row["feed_formula"])[1]))
    for row in ordered_rows:
        x = row["feed_formula"]
        x_salt_free = salt_free_from_formula(x)
        out.append(
            {
                "tie_line": row["tie_line"],
                "temperature_K": row["temperature_K"],
                "salt_wtfrac": row["salt_wtfrac"],
                "x_water_total": float(x[0]),
                "x_ethanol_total": float(x[1]),
                "x_isobutanol_total": float(x[2]),
                "x_nacl_total": float(x[3]),
                "x_water_salt_free": float(x_salt_free[0]),
                "x_ethanol_salt_free": float(x_salt_free[1]),
                "x_isobutanol_salt_free": float(x_salt_free[2]),
                "ethanol_total_feed_wtfrac": (
                    float(row["ethanol_total_feed_wtfrac"])
                    if row["ethanol_total_feed_wtfrac"] not in ("", None)
                    else ""
                ),
                "source": row["source"],
            }
        )
    return out


def write_case_data(fig_dir: Path, exp_rows: list[dict], model_rows: list[dict]) -> None:
    fieldnames = [
        "tie_line",
        "phase",
        "temperature_K",
        "salt_wtfrac",
        "x_water",
        "x_ethanol",
        "x_isobutanol",
        "x_nacl",
        "source",
    ]
    write_csv_rows(
        out_path(fig_dir, "experimental_tielines.csv"), fieldnames, _phase_rows_for_csv(exp_rows, "paper_table")
    )
    model_fieldnames = fieldnames[:-1] + ["beta", "residual_norm", "objective", "converged", "source"]
    write_csv_rows(out_path(fig_dir, "model_tielines.csv"), model_fieldnames, _model_rows_for_csv(model_rows))


def plot_lle_figure(fig_dir: Path, figure_number: int, temperature_k: float, salt_wt: float) -> None:
    configure_style()
    exp_rows = _experimental_rows(salt_wt, temperature_k)
    feed_rows = _digitized_feed_rows_for_figure(figure_number, temperature_k, salt_wt) or _derived_feed_rows(
        salt_wt, temperature_k
    )
    force_recompute = os.environ.get("KHUDAIDA_FORCE_RECOMPUTE", "").strip().lower() in {"1", "true", "yes", "on"}
    model_rows = get_or_build_model_rows(fig_dir, exp_rows, feed_rows=feed_rows, force_recompute=force_recompute)
    paper_rows = _paper_epcsaft_digitized_rows(fig_dir, temperature_k, salt_wt)
    write_csv_rows(
        out_path(fig_dir, "feed_compositions.csv"),
        [
            "tie_line",
            "temperature_K",
            "salt_wtfrac",
            "x_water_total",
            "x_ethanol_total",
            "x_isobutanol_total",
            "x_nacl_total",
            "x_water_salt_free",
            "x_ethanol_salt_free",
            "x_isobutanol_salt_free",
            "ethanol_total_feed_wtfrac",
            "source",
        ],
        _feed_rows_for_csv(feed_rows),
    )

    fig, ax = plt.subplots(figsize=(6.1, 5.9))
    fig.subplots_adjust(left=0.08, right=0.98, top=0.98, bottom=0.18)
    _draw_ternary_axes(ax)
    _plot_tie_lines(ax, exp_rows, BLACK, "o", "Exp.", linestyle="-")
    valid_model_rows = [
        row
        for row in model_rows
        if np.all(np.isfinite(row["organic_formula"])) and np.all(np.isfinite(row["aqueous_formula"]))
    ]
    _plot_tie_lines(ax, valid_model_rows, RED, "o", "model ePC-SAFT", linestyle="--")
    _plot_formula_series(ax, paper_rows, "organic_formula", BLUE, "s", "paper ePC-SAFT", linestyle="-")
    _plot_feed_points(ax, feed_rows, label="Feed")
    ax.legend(loc="upper right", fontsize=8)
    if valid_model_rows:
        model_caption = f"red dashed (accepted native ePC-SAFT, {len(valid_model_rows)}/{len(model_rows)} tie-lines)"
    else:
        model_caption = "native ePC-SAFT rejected all model tie-lines; no red model tie-lines are drawn"
    paper_caption = f"blue squares (digitized paper ePC-SAFT organic branch, {len(paper_rows)} points)"
    add_figure_caption(
        fig,
        f"Figure {figure_number}. LLE for the system water + ethanol + isobutanol + {int(round(salt_wt * 100))} wt % NaCl at {temperature_k:.2f} K and atmospheric pressure expressed as salt-free composition: black (exp), {model_caption}, {paper_caption}, and green (feed compositions).",
    )
    save_figure(fig, fig_dir / f"figure_{figure_number}.png")
    plt.close(fig)

    fig_scaled, ax_scaled = plt.subplots(figsize=(6.1, 5.9))
    fig_scaled.subplots_adjust(left=0.08, right=0.98, top=0.98, bottom=0.18)
    scaled_ranges = _scaled_ranges_for_figure(figure_number, exp_rows, valid_model_rows, feed_rows)
    _draw_display_scaled_axes(
        ax_scaled,
        water_min=scaled_ranges["water_min"],
        water_max=scaled_ranges["water_max"],
        ethanol_min=scaled_ranges["ethanol_min"],
        ethanol_max=scaled_ranges["ethanol_max"],
        isobutanol_min=scaled_ranges["isobutanol_min"],
        isobutanol_max=scaled_ranges["isobutanol_max"],
        tick_format=".2f",
    )
    scaled_xy_transform = lambda x_formula, ranges=scaled_ranges: _display_scaled_xy_from_formula(
        x_formula,
        water_min=ranges["water_min"],
        water_max=ranges["water_max"],
        ethanol_min=ranges["ethanol_min"],
        ethanol_max=ranges["ethanol_max"],
        isobutanol_min=ranges["isobutanol_min"],
        isobutanol_max=ranges["isobutanol_max"],
    )
    _plot_phase_points(
        ax_scaled,
        exp_rows,
        "organic",
        BLACK,
        "o",
        label="Exp. organic",
        xy_transform=scaled_xy_transform,
        display_ranges=scaled_ranges,
    )
    _plot_phase_points(
        ax_scaled,
        valid_model_rows,
        "organic",
        RED,
        "o",
        label="model ePC-SAFT organic",
        xy_transform=scaled_xy_transform,
        display_ranges=scaled_ranges,
    )
    _plot_phase_points(
        ax_scaled,
        paper_rows,
        "organic",
        BLUE,
        "s",
        label="paper ePC-SAFT organic",
        xy_transform=scaled_xy_transform,
        display_ranges=scaled_ranges,
    )
    _plot_feed_points(
        ax_scaled,
        feed_rows,
        label="Feed",
        xy_transform=scaled_xy_transform,
        display_ranges=scaled_ranges,
    )
    ax_scaled.legend(loc="upper right", fontsize=8)
    add_figure_caption(
        fig_scaled,
        f"Figure {figure_number} (scaled). Organic-phase LLE compositions for the system water + ethanol + isobutanol + {int(round(salt_wt * 100))} wt % NaCl at {temperature_k:.2f} K and atmospheric pressure expressed as salt-free composition, zoomed to water {scaled_ranges['water_min']:.2f}-{scaled_ranges['water_max']:.2f}, ethanol {scaled_ranges['ethanol_min']:.2f}-{scaled_ranges['ethanol_max']:.2f}, and isobutanol {scaled_ranges['isobutanol_min']:.2f}-{scaled_ranges['isobutanol_max']:.2f}: black (exp organic), red ({len(valid_model_rows)}/{len(model_rows)} accepted native ePC-SAFT organic compositions), blue (paper ePC-SAFT organic), and green (feed compositions).",
    )
    save_figure(fig_scaled, fig_dir / f"figure_{figure_number}_scaled.png")
    plt.close(fig_scaled)


def _no_salt_digitized_rows() -> list[dict]:
    rows = []
    for idx, item in enumerate(NO_SALT_293_DIGITIZED, start=1):
        x2_aq = float(item["x_ethanol_aq"])
        x3_aq = float(item["x_isobutanol_aq"])
        x1_aq = 1.0 - x2_aq - x3_aq
        D = float(item["distribution"])
        S = float(item["separation"])
        x2_org = D * x2_aq
        x1_org = x1_aq * D / S
        x3_org = max(1.0 - x1_org - x2_org, 1.0e-6)
        org = np.asarray([x1_org, x2_org, x3_org, 0.0], dtype=float)
        org = org / np.sum(org)
        aq = np.asarray([x1_aq, x2_aq, x3_aq, 0.0], dtype=float)
        aq = aq / np.sum(aq)
        rows.append(
            {
                "tie_line": idx,
                "temperature_K": 293.15,
                "salt_wtfrac": 0.0,
                "organic_formula": org,
                "aqueous_formula": aq,
            }
        )
    return rows


def plot_figure_1(fig_dir: Path) -> None:
    configure_style()
    no_salt_rows = _no_salt_digitized_rows()
    five_rows = _experimental_rows(0.05, 293.15)
    ten_rows = _experimental_rows(0.10, 293.15)
    fieldnames = [
        "tie_line",
        "phase",
        "temperature_K",
        "salt_wtfrac",
        "x_water",
        "x_ethanol",
        "x_isobutanol",
        "x_nacl",
        "source",
    ]
    write_csv_rows(
        out_path(fig_dir, "without_nacl_29315_digitized.csv"),
        fieldnames,
        _phase_rows_for_csv(no_salt_rows, "digitized_local_paper"),
    )
    write_csv_rows(
        out_path(fig_dir, "with_5wt_nacl_29315.csv"), fieldnames, _phase_rows_for_csv(five_rows, "paper_table")
    )
    write_csv_rows(
        out_path(fig_dir, "with_10wt_nacl_29315.csv"), fieldnames, _phase_rows_for_csv(ten_rows, "paper_table")
    )

    fig, ax = plt.subplots(figsize=(6.1, 5.9))
    fig.subplots_adjust(left=0.08, right=0.98, top=0.98, bottom=0.18)
    _draw_ternary_axes(ax)
    _plot_tie_lines(ax, no_salt_rows, BLACK, "o", "without NaCl", linestyle="-")
    _plot_tie_lines(ax, five_rows, RED, "o", "5 wt% NaCl", linestyle="-")
    _plot_tie_lines(ax, ten_rows, BLUE, "o", "10 wt% NaCl", linestyle="-")
    add_figure_caption(
        fig,
        "Figure 1. LLE for the system water + ethanol + isobutanol without and with NaCl at 293.15 K and atmospheric pressure expressed as salt-free composition: black (without NaCl), red (5 wt % NaCl), and blue (10 wt % NaCl).",
    )
    save_figure(fig, fig_dir / "figure_1.png")
    plt.close(fig)


def _metric_from_row(row: dict) -> dict:
    org = row["organic_formula"]
    aq = row["aqueous_formula"]
    distribution = float(org[1] / aq[1]) if aq[1] > 0.0 else np.nan
    separation = float((org[1] / aq[1]) / (org[0] / aq[0])) if aq[1] > 0.0 and org[0] > 0.0 else np.nan
    return {
        "temperature_K": row["temperature_K"],
        "salt_wtfrac": row["salt_wtfrac"],
        "x_ethanol_aq": float(aq[1]),
        "distribution": distribution,
        "separation": separation,
    }


def _metric_dataset() -> list[dict]:
    rows = []
    for salt_wt, temperature_k in sorted(source_experimental_cases()):
        for exp_row in _experimental_rows(salt_wt, temperature_k):
            rows.append({**_metric_from_row(exp_row), "source": "paper_table"})
    for item in NO_SALT_293_DIGITIZED:
        rows.append(
            {
                "temperature_K": 293.15,
                "salt_wtfrac": 0.0,
                "x_ethanol_aq": float(item["x_ethanol_aq"]),
                "distribution": float(item["distribution"]),
                "separation": float(item["separation"]),
                "source": "digitized_local_paper",
            }
        )
    rows.sort(key=lambda item: (item["temperature_K"], item["salt_wtfrac"], item["x_ethanol_aq"]))
    return rows


def _metric_color(temperature_k: float) -> str:
    if abs(temperature_k - 293.15) < 1e-6:
        return BLACK
    if abs(temperature_k - 303.15) < 1e-6:
        return RED
    return BLUE


def _metric_marker(salt_wt: float) -> str:
    if salt_wt <= 1e-12:
        return "s"
    if abs(salt_wt - 0.05) < 1e-9:
        return "o"
    return "^"


def _plot_metric_figure(fig_dir: Path, figure_number: int, y_key: str, y_label: str, y_max: float) -> None:
    configure_style()
    rows = _metric_dataset()
    write_csv_rows(
        out_path(fig_dir, f"figure_{figure_number}_metrics.csv"),
        ["temperature_K", "salt_wtfrac", "x_ethanol_aq", "distribution", "separation", "source"],
        rows,
    )

    fig, ax = plt.subplots(figsize=(6.0, 4.4))
    fig.subplots_adjust(left=0.13, right=0.98, top=0.97, bottom=0.24)
    for row in rows:
        ax.scatter(
            row["x_ethanol_aq"],
            row[y_key],
            color=_metric_color(row["temperature_K"]),
            marker=_metric_marker(row["salt_wtfrac"]),
            s=22,
            linewidths=0.6,
        )
    ax.set_xlim(0.0, 0.042)
    ax.set_ylim(0.0, y_max)
    ax.set_xlabel(r"$x$ ethanol in aqueous phase")
    ax.set_ylabel(y_label)
    if figure_number == 8:
        caption = "Figure 8. Separation factor of water over ethanol for the mixture water + ethanol + isobutanol + NaCl; black (293.15 K), red (303.15 K), blue (313.15 K), squares (without NaCl), circles (5 wt % NaCl), and triangles (10 wt % NaCl)."
    else:
        caption = "Figure 9. Distribution coefficient of ethanol over water of the mixture water + ethanol + isobutanol + NaCl; black (293.15 K), red (303.15 K), blue (313.15 K), squares (without NaCl), circles (5 wt % NaCl), and triangles (10 wt % NaCl)."
    add_figure_caption(fig, caption, left=0.13, y=0.02)
    save_figure(fig, fig_dir / f"figure_{figure_number}.png")
    plt.close(fig)


def plot_figure_8(fig_dir: Path) -> None:
    _plot_metric_figure(fig_dir, 8, "separation", "Separation Factor", 260.0)


def plot_figure_9(fig_dir: Path) -> None:
    _plot_metric_figure(fig_dir, 9, "distribution", "Distribution Coefficient of Ethanol", 70.0)


def _supporting_panel_rows(figure_number: int) -> list[dict]:
    rows = []
    for panel, main_figure, temperature_k, salt_wt in SUPPORTING_FIGURE_PANELS[figure_number]:
        exp_rows = _experimental_rows(salt_wt, temperature_k)
        feed_rows = _digitized_feed_rows_for_figure(main_figure, temperature_k, salt_wt) or _derived_feed_rows(
            salt_wt, temperature_k
        )
        main_fig_dir = figure_root(main_figure)
        model_rows = get_or_build_model_rows(main_fig_dir, exp_rows, feed_rows)
        paper_rows = _paper_epcsaft_digitized_rows(
            main_fig_dir, temperature_k, salt_wt
        )
        rows.append(
            {
                "panel": panel,
                "main_figure": main_figure,
                "temperature_K": temperature_k,
                "salt_wtfrac": salt_wt,
                "experimental": exp_rows,
                "feed": feed_rows,
                "model": model_rows,
                "paper": paper_rows,
            }
        )
    return rows


def write_supporting_figure_data(fig_dir: Path, figure_number: int) -> None:
    if figure_number not in SUPPORTING_FIGURE_PANELS:
        raise ValueError(f"Unsupported Khudaida supporting figure: S{figure_number}")
    phase_rows: list[dict] = []
    feed_rows_out: list[dict] = []
    for panel in _supporting_panel_rows(figure_number):
        for source, rows, phase_key in (
            ("experimental_table_3_4", panel["experimental"], None),
            ("package_public_native_electrolyte_lle", panel["model"], None),
            ("paper_epcsaft_digitized", panel["paper"], "organic_formula"),
        ):
            for row in rows:
                if phase_key is not None:
                    phases = (("organic", row[phase_key]),)
                else:
                    phases = (("organic", row["organic_formula"]), ("aqueous", row["aqueous_formula"]))
                for phase, x in phases:
                    phase_rows.append(
                        {
                            "panel": panel["panel"],
                            "main_figure": panel["main_figure"],
                            "tie_line": row["tie_line"],
                            "phase": phase,
                            "temperature_K": panel["temperature_K"],
                            "salt_wtfrac": panel["salt_wtfrac"],
                            "x_water": float(x[0]),
                            "x_ethanol": float(x[1]),
                            "x_isobutanol": float(x[2]),
                            "x_nacl": float(x[3]),
                            "source": source,
                        }
                    )
        for row in panel["feed"]:
            x = row["feed_formula"]
            salt_free = salt_free_from_formula(x)
            feed_rows_out.append(
                {
                    "panel": panel["panel"],
                    "main_figure": panel["main_figure"],
                    "tie_line": row["tie_line"],
                    "temperature_K": panel["temperature_K"],
                    "salt_wtfrac": panel["salt_wtfrac"],
                    "x_water_total": float(x[0]),
                    "x_ethanol_total": float(x[1]),
                    "x_isobutanol_total": float(x[2]),
                    "x_nacl_total": float(x[3]),
                    "x_water_salt_free": float(salt_free[0]),
                    "x_ethanol_salt_free": float(salt_free[1]),
                    "x_isobutanol_salt_free": float(salt_free[2]),
                    "source": row["source"],
                }
            )
    write_csv_rows(
        out_path(fig_dir, f"figure_s{figure_number}_phase_points.csv"),
        [
            "panel",
            "main_figure",
            "tie_line",
            "phase",
            "temperature_K",
            "salt_wtfrac",
            "x_water",
            "x_ethanol",
            "x_isobutanol",
            "x_nacl",
            "source",
        ],
        phase_rows,
    )
    write_csv_rows(
        out_path(fig_dir, f"figure_s{figure_number}_feed_points.csv"),
        [
            "panel",
            "main_figure",
            "tie_line",
            "temperature_K",
            "salt_wtfrac",
            "x_water_total",
            "x_ethanol_total",
            "x_isobutanol_total",
            "x_nacl_total",
            "x_water_salt_free",
            "x_ethanol_salt_free",
            "x_isobutanol_salt_free",
            "source",
        ],
        feed_rows_out,
    )


def plot_supporting_figure_grid(fig_dir: Path, figure_number: int) -> None:
    configure_style()
    panel_rows = _supporting_panel_rows(figure_number)
    write_supporting_figure_data(fig_dir, figure_number)
    fig, axes = plt.subplots(1, 3, figsize=(14.0, 4.9))
    fig.subplots_adjust(left=0.04, right=0.99, top=0.88, bottom=0.18, wspace=0.18)
    for ax, panel in zip(axes, panel_rows, strict=True):
        _draw_ternary_axes(ax)
        model_rows = [
            row
            for row in panel["model"]
            if np.all(np.isfinite(row["organic_formula"])) and np.all(np.isfinite(row["aqueous_formula"]))
        ]
        _plot_tie_lines(ax, panel["experimental"], BLACK, "o", "Exp.", linestyle="-")
        _plot_tie_lines(ax, model_rows, RED, "o", "package ePC-SAFT", linestyle="--")
        _plot_formula_series(ax, panel["paper"], "organic_formula", BLUE, "s", "paper ePC-SAFT", linestyle="-")
        _plot_feed_points(ax, panel["feed"], label="Feed", marker="^")
        ax.set_title(f"({panel['panel']}) {panel['temperature_K']:.2f} K", fontsize=10)
    axes[0].legend(loc="upper right", fontsize=7)
    salt_label = "5 wt % NaCl" if figure_number == 2 else "10 wt % NaCl"
    add_figure_caption(
        fig,
        f"Figure S{figure_number}. LLE for water + ethanol + isobutanol + {salt_label} at 293.15, 303.15, and 313.15 K: black (experimental), red dashed (package public native electrolyte_lle solve), blue (digitized paper ePC-SAFT), and green triangles (feed compositions).",
        left=0.04,
        y=0.02,
    )
    save_figure(fig, fig_dir / f"figure_s{figure_number}.png")
    plt.close(fig)


def _aad_summary(exp_rows: list[dict], model_rows: list[dict]) -> dict:
    valid = [
        row
        for row in model_rows
        if np.all(np.isfinite(row["organic_formula"])) and np.all(np.isfinite(row["aqueous_formula"]))
    ]
    if not valid:
        nan4 = (np.nan, np.nan, np.nan, np.nan)
        return {"organic": nan4, "aqueous": nan4, "grand": np.nan}
    exp_map = {row["tie_line"]: row for row in exp_rows}
    organic_delta = []
    aqueous_delta = []
    for row in valid:
        exp_row = exp_map[row["tie_line"]]
        organic_delta.append(np.abs(row["organic_formula"] - exp_row["organic_formula"]))
        aqueous_delta.append(np.abs(row["aqueous_formula"] - exp_row["aqueous_formula"]))
    organic_arr = np.vstack(organic_delta)
    aqueous_arr = np.vstack(aqueous_delta)
    grand = float((organic_arr.sum() + aqueous_arr.sum()) / (8.0 * organic_arr.shape[0]))
    return {
        "organic": tuple(np.mean(organic_arr, axis=0).tolist()),
        "aqueous": tuple(np.mean(aqueous_arr, axis=0).tolist()),
        "grand": grand,
    }


def _table_rows_for_png(salt_wt: float) -> list[list[str]]:
    rows = []
    temps = sorted({key[1] for key in source_experimental_cases() if abs(key[0] - salt_wt) < 1e-9})
    figure_number_map = {
        (0.05, 293.15): 2,
        (0.05, 303.15): 3,
        (0.05, 313.15): 4,
        (0.10, 293.15): 5,
        (0.10, 303.15): 6,
        (0.10, 313.15): 7,
    }
    force_recompute = os.environ.get("KHUDAIDA_FORCE_RECOMPUTE", "").strip().lower() in {"1", "true", "yes", "on"}
    for temperature_k in temps:
        exp_rows = _experimental_rows(salt_wt, temperature_k)
        feed_rows = _digitized_feed_rows_for_figure(
            figure_number_map[(salt_wt, temperature_k)], temperature_k, salt_wt
        ) or _derived_feed_rows(salt_wt, temperature_k)
        fig_dir = figure_root(figure_number_map[(salt_wt, temperature_k)])
        model_rows = get_or_build_model_rows(fig_dir, exp_rows, feed_rows=feed_rows, force_recompute=force_recompute)
        ours = _aad_summary(exp_rows, model_rows)
        paper_epc = EePCSAFT_AAD_REFERENCE[salt_wt][temperature_k]
        paper_enrtl = ENRTL_AAD_REFERENCE[salt_wt][temperature_k]
        for model_name, summary in (
            ("ePC-SAFT (package)", ours),
            ("ePC-SAFT (paper)", paper_epc),
            ("eNRTL (paper)", paper_enrtl),
        ):
            rows.append(
                [
                    f"{temperature_k:.2f}",
                    model_name,
                    f"{summary['organic'][0]:.4f}",
                    f"{summary['organic'][1]:.4f}",
                    f"{summary['organic'][2]:.4f}",
                    f"{summary['organic'][3]:.4f}",
                    f"{summary['aqueous'][0]:.4f}",
                    f"{summary['aqueous'][1]:.4f}",
                    f"{summary['aqueous'][2]:.4f}",
                    f"{summary['aqueous'][3]:.4f}",
                    f"{summary['grand']:.4f}",
                ]
            )
    return rows


def plot_tables_9_10(fig_dir: Path) -> None:
    configure_style()
    columns = [
        "T / K",
        "Model",
        "Org x1",
        "Org x2",
        "Org x3",
        "Org x4",
        "Aq x1",
        "Aq x2",
        "Aq x3",
        "Aq x4",
        "Grand AAD",
    ]
    for table_number, salt_wt in ((9, 0.05), (10, 0.10)):
        rows = _table_rows_for_png(salt_wt)
        csv_rows = [dict(zip(columns, row)) for row in rows]
        write_csv_rows(out_path(fig_dir, f"table_{table_number}.csv"), columns, csv_rows)
        fig, ax = plt.subplots(figsize=(13.0, 3.0 + 0.28 * len(rows)))
        ax.axis("off")
        table = ax.table(cellText=rows, colLabels=columns, loc="center", cellLoc="center")
        table.auto_set_font_size(False)
        table.set_fontsize(8)
        table.scale(1.0, 1.2)
        ax.set_title(
            f"Table {table_number}: AAD comparison for {int(round(salt_wt * 100))} wt% NaCl", fontsize=11, pad=12
        )
        save_figure(fig, fig_dir / f"table_{table_number}.png")
        plt.close(fig)


def write_provenance_notes() -> None:
    notes = """# 2026 Khudaida analysis provenance

- Tables 3 and 4 in the local Khudaida markdown/PDF are treated as the canonical experimental tie-line source for Figures 2-7 and for the salted points in Figures 8-9.
- Analysis-owned source CSVs in `data/input/` retain paper Tables 3-7 and Supporting Information Tables S1-S4. The runtime dataset `2026_Khudaida` uses the paper's Table 5 pure-component parameters, Table 6 dielectric constants, Table 7 binary interaction parameters, and the Figiel 2025 SSM+DS Born option family.
- Figures S2 and S3 are represented as three-panel figure-owned workflows that reuse the public non-reactive `mix.equilibrium(kind="electrolyte_lle", ...)` path and snapshot the exact plotted phase/feed points.
- Figure 1 salt-free data and the no-salt points in Figures 8-9 were reconstructed from the local paper figures because the Zotero baseline source remained inaccessible in this session.
- The no-salt baseline is therefore marked as `digitized_local_paper` in the emitted CSV files.
- Tables 9 and 10 include package-generated ePC-SAFT AAD values and paper-copied eNRTL/ePC-SAFT reference values for comparison.
- The package model evidence is non-reactive electrolyte LLE through the public native Ipopt route with the `2026_Khudaida` Born SSM+DS dataset options.
"""
    (ANALYSIS_ROOT / "provenance_notes.md").write_text(notes, encoding="utf-8")
