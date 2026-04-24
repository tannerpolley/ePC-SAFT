from __future__ import annotations

import csv
import math
import os
import platform
import sys
from pathlib import Path
from typing import Iterable

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parent
REPO_ROOT = ROOT.parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts._env import require_epcsaft_install
from scripts.plot_outputs import paper_validation_output_path

require_epcsaft_install()


def _fast_machine() -> str:
    return os.environ.get("PROCESSOR_ARCHITECTURE", "AMD64")


platform.machine = _fast_machine

import scripts._epcsaft_oop as pcs
from epcsaft.parameters import get_prop_dict


P_REF = 1.0e5
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

SCALED_FIGURE_OVERRIDES: dict[int, dict] = {}


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
    fig.savefig(path, dpi=300, bbox_inches="tight")


def add_figure_caption(fig: plt.Figure, caption: str, *, left: float = 0.09, y: float = 0.02, fontsize: float = 9.0) -> None:
    fig.text(left, y, caption, ha="left", va="bottom", fontsize=fontsize, wrap=True)


def write_csv_rows(path: Path, fieldnames: Iterable[str], rows: Iterable[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fieldnames))
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


EXPERIMENTAL_CASES: dict[tuple[float, float], list[tuple[int, tuple[float, ...], tuple[float, ...]]]] = {
    (0.05, 293.15): [
        (1, (0.3705, 0.1624, 0.4650, 0.0021), (0.9331, 0.0159, 0.0028, 0.0482)),
        (2, (0.3703, 0.1400, 0.4880, 0.0017), (0.9404, 0.0121, 0.0027, 0.0448)),
        (3, (0.3691, 0.1144, 0.5153, 0.0012), (0.9446, 0.0110, 0.0022, 0.0422)),
        (4, (0.3695, 0.0927, 0.5366, 0.0012), (0.9486, 0.0092, 0.0037, 0.0385)),
        (5, (0.3692, 0.0759, 0.5537, 0.0012), (0.9509, 0.0073, 0.0034, 0.0384)),
        (6, (0.3685, 0.0646, 0.5658, 0.0011), (0.9539, 0.0056, 0.0034, 0.0371)),
        (7, (0.3676, 0.0478, 0.5836, 0.0010), (0.9567, 0.0044, 0.0036, 0.0353)),
        (8, (0.3641, 0.0336, 0.6015, 0.0008), (0.9580, 0.0030, 0.0037, 0.0353)),
    ],
    (0.05, 303.15): [
        (1, (0.3906, 0.1597, 0.4480, 0.0017), (0.9345, 0.0149, 0.0025, 0.0481)),
        (2, (0.3874, 0.1371, 0.4740, 0.0015), (0.9397, 0.0124, 0.0029, 0.0450)),
        (3, (0.3887, 0.1105, 0.4995, 0.0013), (0.9437, 0.0099, 0.0036, 0.0428)),
        (4, (0.3849, 0.0885, 0.5258, 0.0008), (0.9474, 0.0081, 0.0044, 0.0401)),
        (5, (0.3821, 0.0759, 0.5410, 0.0010), (0.9489, 0.0069, 0.0049, 0.0393)),
        (6, (0.3768, 0.0623, 0.5601, 0.0008), (0.9514, 0.0055, 0.0050, 0.0381)),
        (7, (0.3753, 0.0476, 0.5767, 0.0004), (0.9550, 0.0043, 0.0041, 0.0366)),
        (8, (0.3740, 0.0316, 0.5940, 0.0004), (0.9572, 0.0026, 0.0046, 0.0356)),
    ],
    (0.05, 313.15): [
        (1, (0.4154, 0.1313, 0.4513, 0.0020), (0.9406, 0.0122, 0.0037, 0.0435)),
        (2, (0.4130, 0.1094, 0.4759, 0.0017), (0.9426, 0.0101, 0.0041, 0.0432)),
        (3, (0.4077, 0.0866, 0.5045, 0.0012), (0.9457, 0.0079, 0.0043, 0.0421)),
        (4, (0.4091, 0.0701, 0.5192, 0.0016), (0.9497, 0.0064, 0.0042, 0.0397)),
        (5, (0.4079, 0.0577, 0.5333, 0.0011), (0.9520, 0.0051, 0.0044, 0.0385)),
        (6, (0.4013, 0.0487, 0.5489, 0.0011), (0.9536, 0.0043, 0.0045, 0.0376)),
        (7, (0.3954, 0.0342, 0.5694, 0.0010), (0.9560, 0.0029, 0.0046, 0.0365)),
    ],
    (0.10, 293.15): [
        (1, (0.2456, 0.0919, 0.6612, 0.0013), (0.9223, 0.0044, 0.0013, 0.0720)),
        (2, (0.2480, 0.0767, 0.6739, 0.0014), (0.9253, 0.0033, 0.0010, 0.0704)),
        (3, (0.2476, 0.0615, 0.6896, 0.0013), (0.9277, 0.0025, 0.0013, 0.0685)),
        (4, (0.2543, 0.0383, 0.7062, 0.0012), (0.9311, 0.0014, 0.0014, 0.0661)),
        (5, (0.2560, 0.0182, 0.7248, 0.0010), (0.9353, 0.0003, 0.0001, 0.0643)),
    ],
    (0.10, 303.15): [
        (1, (0.2587, 0.0885, 0.6521, 0.0007), (0.9235, 0.0049, 0.0015, 0.0701)),
        (2, (0.2566, 0.0735, 0.6593, 0.0006), (0.9249, 0.0038, 0.0016, 0.0697)),
        (3, (0.2632, 0.0554, 0.6811, 0.0003), (0.9280, 0.0025, 0.0016, 0.0679)),
        (4, (0.2681, 0.0359, 0.6957, 0.0003), (0.9309, 0.0014, 0.0012, 0.0665)),
        (5, (0.2736, 0.0179, 0.7078, 0.0007), (0.9363, 0.0005, 0.0003, 0.0629)),
    ],
    (0.10, 313.15): [
        (1, (0.2684, 0.1031, 0.6268, 0.0017), (0.9210, 0.0057, 0.0013, 0.0720)),
        (2, (0.2669, 0.0954, 0.6358, 0.0019), (0.9242, 0.0042, 0.0008, 0.0708)),
        (3, (0.2693, 0.0759, 0.6535, 0.0013), (0.9271, 0.0028, 0.0008, 0.0694)),
        (4, (0.2695, 0.0615, 0.6677, 0.0013), (0.9283, 0.0022, 0.0013, 0.0682)),
        (5, (0.2733, 0.0408, 0.6845, 0.0014), (0.9295, 0.0015, 0.0014, 0.0676)),
        (6, (0.2831, 0.0176, 0.6982, 0.0011), (0.9333, 0.0005, 0.0010, 0.0652)),
    ],
}

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
        293.15: {"grand": 0.0161, "organic": (0.0287, 0.0031, 0.0313, 0.0009), "aqueous": (0.0482, 0.0064, 0.0021, 0.0083)},
        303.15: {"grand": 0.0168, "organic": (0.0305, 0.0027, 0.0336, 0.0005), "aqueous": (0.0495, 0.0061, 0.0026, 0.0088)},
        313.15: {"grand": 0.0203, "organic": (0.0245, 0.0033, 0.0276, 0.0006), "aqueous": (0.0487, 0.0052, 0.0025, 0.0081)},
    },
    0.10: {
        293.15: {"grand": 0.0342, "organic": (0.0729, 0.0096, 0.0739, 0.0005), "aqueous": (0.0975, 0.0020, 0.0007, 0.0161)},
        303.15: {"grand": 0.0341, "organic": (0.0760, 0.0021, 0.0774, 0.0005), "aqueous": (0.0973, 0.0023, 0.0008, 0.0165)},
        313.15: {"grand": 0.0357, "organic": (0.0801, 0.0039, 0.0783, 0.0002), "aqueous": (0.1021, 0.0025, 0.0003, 0.0179)},
    },
}

ENRTL_AAD_REFERENCE = {
    0.05: {
        293.15: {"grand": 0.0046, "organic": (0.0059, 0.0048, 0.0010, 0.0011), "aqueous": (0.0102, 0.0025, 0.0080, 0.0029)},
        303.15: {"grand": 0.0047, "organic": (0.0053, 0.0040, 0.0015, 0.0007), "aqueous": (0.0120, 0.0029, 0.0097, 0.0012)},
        313.15: {"grand": 0.0078, "organic": (0.0147, 0.0117, 0.0041, 0.0011), "aqueous": (0.0128, 0.0030, 0.0126, 0.0027)},
    },
    0.10: {
        293.15: {"grand": 0.0086, "organic": (0.0010, 0.0026, 0.0022, 0.0007), "aqueous": (0.0312, 0.0008, 0.0062, 0.0241)},
        303.15: {"grand": 0.0083, "organic": (0.0011, 0.0022, 0.0033, 0.0022), "aqueous": (0.0289, 0.0006, 0.0079, 0.0204)},
        313.15: {"grand": 0.0091, "organic": (0.0011, 0.0035, 0.0028, 0.0004), "aqueous": (0.0323, 0.0020, 0.0098, 0.0205)},
    },
}


def _experimental_rows(salt_wt: float, temperature_k: float) -> list[dict]:
    rows = []
    for tie_line, organic, aqueous in EXPERIMENTAL_CASES[(salt_wt, temperature_k)]:
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
    source_path = ROOT / f"figure_{figure_number}" / "data" / "feed_compositions_digitized.csv"
    if not source_path.exists():
        return None
    rows = []
    with source_path.open("r", newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        for idx, raw in enumerate(reader, start=1):
            ethanol = float(raw["x_ethanol_salt_free"])
            isobutanol = float(raw["x_isobutanol_salt_free"])
            water = float(raw["x_water_salt_free"])
            feed_formula = build_feed_formula_from_salt_free_molefractions(np.asarray([water, ethanol, isobutanol], dtype=float), salt_wt)
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
            ax.text(left_point[0], left_point[1], f"{1.0 - frac:.1f}", color=WATER_COLOR, ha="right", va="center", fontsize=8)
            # Right edge shows ethanol increasing from 0 at the lower-right corner to 1 at the apex.
            right_point = np.asarray([1.0 - 0.5 * frac, SQRT3_OVER_2 * frac], dtype=float) + right_tick_offset * right_normal
            ax.text(right_point[0], right_point[1], f"{frac:.1f}", color=ETHANOL_COLOR, ha="left", va="center", fontsize=8)
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
            ax.text(frac, -0.03, format(isobutanol, tick_format), color=ISOBUTANOL_COLOR, ha="center", va="top", fontsize=7)
        if 0.0 < frac < 1.0:
            left_point = np.asarray([0.5 * frac, SQRT3_OVER_2 * frac], dtype=float) + left_tick_offset * left_normal
            ax.text(left_point[0], left_point[1], format(water, tick_format), color=WATER_COLOR, ha="right", va="center", fontsize=7)
            right_point = np.asarray([1.0 - 0.5 * frac, SQRT3_OVER_2 * frac], dtype=float) + right_tick_offset * right_normal
            ax.text(right_point[0], right_point[1], format(ethanol, tick_format), color=ETHANOL_COLOR, ha="left", va="center", fontsize=7)
    left_bottom = np.asarray([0.0, 0.0], dtype=float) + left_tick_offset * left_normal
    left_top = np.asarray([0.5, SQRT3_OVER_2], dtype=float) + left_tick_offset * left_normal
    right_bottom = np.asarray([1.0, 0.0], dtype=float) + right_tick_offset * right_normal
    right_top = np.asarray([0.5, SQRT3_OVER_2], dtype=float) + right_tick_offset * right_normal
    ax.text(left_bottom[0], left_bottom[1], format(1.0, tick_format), color=WATER_COLOR, ha="right", va="center", fontsize=7)
    ax.text(left_top[0], left_top[1], format(water_min, tick_format), color=WATER_COLOR, ha="right", va="center", fontsize=7)
    ax.text(right_bottom[0], right_bottom[1], format(0.0, tick_format), color=ETHANOL_COLOR, ha="left", va="center", fontsize=7)
    ax.text(right_top[0], right_top[1], format(ethanol_max, tick_format), color=ETHANOL_COLOR, ha="left", va="center", fontsize=7)
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
            ax.text(frac, -0.03, format(isobutanol, tick_format), color=ISOBUTANOL_COLOR, ha="center", va="top", fontsize=7)
        if 0.0 < frac < 1.0:
            left_point = np.asarray([0.5 * frac, SQRT3_OVER_2 * frac], dtype=float) + left_tick_offset * left_normal
            ax.text(left_point[0], left_point[1], format(water, tick_format), color=WATER_COLOR, ha="right", va="center", fontsize=7)
            right_point = np.asarray([1.0 - 0.5 * frac, SQRT3_OVER_2 * frac], dtype=float) + right_tick_offset * right_normal
            ax.text(right_point[0], right_point[1], format(ethanol, tick_format), color=ETHANOL_COLOR, ha="left", va="center", fontsize=7)
    left_bottom = np.asarray([0.0, 0.0], dtype=float) + left_tick_offset * left_normal
    left_top = np.asarray([0.5, SQRT3_OVER_2], dtype=float) + left_tick_offset * left_normal
    right_bottom = np.asarray([1.0, 0.0], dtype=float) + right_tick_offset * right_normal
    right_top = np.asarray([0.5, SQRT3_OVER_2], dtype=float) + right_tick_offset * right_normal
    ax.text(left_bottom[0], left_bottom[1], format(water_max, tick_format), color=WATER_COLOR, ha="right", va="center", fontsize=7)
    ax.text(left_top[0], left_top[1], format(water_min, tick_format), color=WATER_COLOR, ha="right", va="center", fontsize=7)
    ax.text(right_bottom[0], right_bottom[1], format(ethanol_min, tick_format), color=ETHANOL_COLOR, ha="left", va="center", fontsize=7)
    ax.text(right_top[0], right_top[1], format(ethanol_max, tick_format), color=ETHANOL_COLOR, ha="left", va="center", fontsize=7)
    ax.set_xlim(-0.12, 1.14)
    ax.set_ylim(-0.10, SQRT3_OVER_2 + 0.08)
    ax.set_aspect("equal")
    ax.axis("off")


def _scaled_axis_scale_from_rows(exp_rows: list[dict], model_rows: list[dict], feed_rows: list[dict]) -> float:
    ethanol_values = []
    for row in exp_rows:
        ethanol_values.extend([float(salt_free_from_formula(row["organic_formula"])[1]), float(salt_free_from_formula(row["aqueous_formula"])[1])])
    for row in model_rows:
        if np.all(np.isfinite(row["organic_formula"])) and np.all(np.isfinite(row["aqueous_formula"])):
            ethanol_values.extend([float(salt_free_from_formula(row["organic_formula"])[1]), float(salt_free_from_formula(row["aqueous_formula"])[1])])
    for row in feed_rows:
        ethanol_values.append(float(salt_free_from_formula(row["feed_formula"])[1]))
    max_ethanol = max(ethanol_values) if ethanol_values else 0.2
    return max(0.1, min(1.0, math.ceil(max_ethanol * 10.0) / 10.0))


def _scaled_water_min_from_rows(exp_rows: list[dict], model_rows: list[dict], feed_rows: list[dict]) -> float:
    water_values = []
    for row in exp_rows:
        water_values.extend([float(salt_free_from_formula(row["organic_formula"])[0]), float(salt_free_from_formula(row["aqueous_formula"])[0])])
    for row in model_rows:
        if np.all(np.isfinite(row["organic_formula"])) and np.all(np.isfinite(row["aqueous_formula"])):
            water_values.extend([float(salt_free_from_formula(row["organic_formula"])[0]), float(salt_free_from_formula(row["aqueous_formula"])[0])])
    for row in feed_rows:
        water_values.append(float(salt_free_from_formula(row["feed_formula"])[0]))
    min_water = min(water_values) if water_values else 0.0
    return max(0.0, min(0.9, math.floor(min_water * 10.0) / 10.0))


def _display_scaled_ranges_from_rows(exp_rows: list[dict], model_rows: list[dict], feed_rows: list[dict]) -> dict[str, float]:
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
    for key_min, key_max in (("water_min", "water_max"), ("ethanol_min", "ethanol_max"), ("isobutanol_min", "isobutanol_max")):
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
        ax.scatter([aq_xy[0], org_xy[0]], [aq_xy[1], org_xy[1]], s=markersize, marker=marker, facecolors=color, edgecolors=color, linewidths=0.5, zorder=4)


def _model_objective(exp_row: dict, pred_org: np.ndarray, pred_aq: np.ndarray) -> float:
    return float(np.sum(np.abs(pred_org - exp_row["organic_formula"])) + np.sum(np.abs(pred_aq - exp_row["aqueous_formula"])))


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
    xy_transform=ternary_xy_from_formula,
) -> None:
    if not rows:
        return
    xs = []
    ys = []
    for row in rows:
        x_coord, y_coord = xy_transform(row["feed_formula"])
        xs.append(x_coord)
        ys.append(y_coord)
    ax.scatter(xs, ys, s=markersize, marker=marker, facecolors=color, edgecolors=color, linewidths=0.5, zorder=5)


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


def _solve_formula_feed(temperature_k: float, feed_formula: np.ndarray, seed_formula_candidates: list[np.ndarray]) -> dict | None:
    z_feed = formula_to_ion_basis(feed_formula)
    params = get_prop_dict(PARAMETER_DATASET, SPECIES, z_feed, temperature_k)
    unique_seed_candidates: list[np.ndarray] = []
    for seed_formula in seed_formula_candidates:
        seed_formula = np.asarray(seed_formula, dtype=float)
        if not np.all(np.isfinite(seed_formula)):
            continue
        seed_formula = seed_formula / np.sum(seed_formula)
        if not any(np.allclose(seed_formula, prior, atol=1e-10) for prior in unique_seed_candidates):
            unique_seed_candidates.append(seed_formula)
    attempts = []
    for seed_formula in unique_seed_candidates[:2]:
        attempts.append(
            {
                "seed_x": formula_to_ion_basis(seed_formula),
                "force_seed_solve": True,
                "tpdf_global_trials": 300,
                "tpdf_local_trials": 120,
                "solver_tol": 1.0e-9,
                "max_nfev": 180,
                "charge_weight": 2500.0,
                "solver_accept_norm": 0.5,
                "split_tol": 1.0e-4,
                "debug": False,
            }
        )
    best = None
    for options in attempts:
        raise NotImplementedError("The legacy multiphase LLE workflow has been removed and will be rewritten later.")
    if best is None:
        fallback_options = {
            "tpdf_global_trials": 300,
            "tpdf_local_trials": 120,
            "solver_tol": 1.0e-9,
            "max_nfev": 180,
            "charge_weight": 2500.0,
            "solver_accept_norm": 0.5,
            "split_tol": 1.0e-4,
            "debug": False,
        }
        raise NotImplementedError("The legacy multiphase LLE workflow has been removed and will be rewritten later.")
    return best


def solve_model_rows(exp_rows: list[dict], feed_rows: list[dict] | None = None) -> list[dict]:
    feed_map = {int(row["tie_line"]): np.asarray(row["feed_formula"], dtype=float) for row in (feed_rows or [])}
    solved = []
    for exp_row in exp_rows:
        best = None
        target_feed_formula = feed_map.get(int(exp_row["tie_line"]))
        midpoint = 0.5 * (exp_row["organic_formula"] + exp_row["aqueous_formula"])
        seed_candidates = [exp_row["organic_formula"], exp_row["aqueous_formula"], midpoint]
        if target_feed_formula is not None:
            seed_candidates.append(target_feed_formula)
        for feed_idx, feed_formula in enumerate(_candidate_formula_feeds(exp_row, target_feed_formula)):
            candidate = _solve_formula_feed(exp_row["temperature_K"], feed_formula, seed_candidates)
            if candidate is None:
                continue
            score = _model_objective(exp_row, candidate["organic_formula"], candidate["aqueous_formula"])
            if best is None or score < best["objective"]:
                best = {"objective": score, **candidate}
            if best is not None and target_feed_formula is not None and feed_idx == 0 and best["converged"]:
                break
        if best is None:
            best = {
                "converged": False,
                "status": None,
                "message": "No converged model candidate.",
                "residual_norm": np.nan,
                "feed_formula": np.full(4, np.nan),
                "organic_formula": np.full(4, np.nan),
                "aqueous_formula": np.full(4, np.nan),
                "beta_organic": np.nan,
                "beta_aqueous": np.nan,
                "split_norm": np.nan,
                "objective": np.nan,
            }
        solved.append({"tie_line": exp_row["tie_line"], "temperature_K": exp_row["temperature_K"], "salt_wtfrac": exp_row["salt_wtfrac"], **best})
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
    return fig_dir / "data" / "model_tielines.csv"


def get_or_build_model_rows(fig_dir: Path, exp_rows: list[dict], feed_rows: list[dict] | None = None, force_recompute: bool = False) -> list[dict]:
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
                    "source": "epcsaft_pressureackage",
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
    fieldnames = ["tie_line", "phase", "temperature_K", "salt_wtfrac", "x_water", "x_ethanol", "x_isobutanol", "x_nacl", "source"]
    write_csv_rows(fig_dir / "data" / "experimental_tielines.csv", fieldnames, _phase_rows_for_csv(exp_rows, "paper_table"))
    model_fieldnames = fieldnames[:-1] + ["beta", "residual_norm", "objective", "converged", "source"]
    write_csv_rows(fig_dir / "data" / "model_tielines.csv", model_fieldnames, _model_rows_for_csv(model_rows))


def plot_lle_figure(fig_dir: Path, figure_number: int, temperature_k: float, salt_wt: float) -> None:
    configure_style()
    exp_rows = _experimental_rows(salt_wt, temperature_k)
    feed_rows = _digitized_feed_rows_for_figure(figure_number, temperature_k, salt_wt) or _derived_feed_rows(salt_wt, temperature_k)
    force_recompute = os.environ.get("KHUDAIDA_FORCE_RECOMPUTE", "").strip().lower() in {"1", "true", "yes", "on"}
    model_rows = get_or_build_model_rows(fig_dir, exp_rows, feed_rows=feed_rows, force_recompute=force_recompute)
    write_csv_rows(
        fig_dir / "data" / "feed_compositions.csv",
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
    valid_model_rows = [row for row in model_rows if np.all(np.isfinite(row["organic_formula"])) and np.all(np.isfinite(row["aqueous_formula"]))]
    _plot_tie_lines(ax, valid_model_rows, RED, "o", "ePC-SAFT", linestyle="--")
    _plot_feed_points(ax, feed_rows)
    add_figure_caption(
        fig,
        f"Figure {figure_number}. LLE for the system water + ethanol + isobutanol + {int(round(salt_wt * 100))} wt % NaCl at {temperature_k:.2f} K and atmospheric pressure expressed as salt-free composition: black (exp), red (ePC-SAFT), and green (feed compositions).",
    )
    save_figure(fig, fig_dir / f"figure_{figure_number}.png")
    plt.close(fig)

    fig_scaled, ax_scaled = plt.subplots(figsize=(6.1, 5.9))
    fig_scaled.subplots_adjust(left=0.08, right=0.98, top=0.98, bottom=0.18)
    water_min = _scaled_water_min_from_rows(exp_rows, valid_model_rows, feed_rows)
    scaled_axis_scale = 1.0 - water_min
    _draw_scaled_ternary_axes(
        ax_scaled,
        water_min=water_min,
        ethanol_max=scaled_axis_scale,
        isobutanol_max=scaled_axis_scale,
        tick_format=".3f",
    )
    scaled_xy_transform = lambda x_formula, scale=scaled_axis_scale: _scaled_xy_from_formula(
        x_formula,
        ethanol_max=scale,
        isobutanol_max=scale,
    )
    _plot_tie_lines(ax_scaled, exp_rows, BLACK, "o", "Exp.", linestyle="-", xy_transform=scaled_xy_transform)
    _plot_tie_lines(ax_scaled, valid_model_rows, RED, "o", "ePC-SAFT", linestyle="--", xy_transform=scaled_xy_transform)
    _plot_feed_points(ax_scaled, feed_rows, xy_transform=scaled_xy_transform)
    add_figure_caption(
        fig_scaled,
        f"Figure {figure_number} (scaled). LLE for the system water + ethanol + isobutanol + {int(round(salt_wt * 100))} wt % NaCl at {temperature_k:.2f} K and atmospheric pressure expressed as salt-free composition, zoomed to the water-rich corner: black (exp), red (ePC-SAFT), and green (feed compositions).",
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
        rows.append({"tie_line": idx, "temperature_K": 293.15, "salt_wtfrac": 0.0, "organic_formula": org, "aqueous_formula": aq})
    return rows


def plot_figure_1(fig_dir: Path) -> None:
    configure_style()
    no_salt_rows = _no_salt_digitized_rows()
    five_rows = _experimental_rows(0.05, 293.15)
    ten_rows = _experimental_rows(0.10, 293.15)
    fieldnames = ["tie_line", "phase", "temperature_K", "salt_wtfrac", "x_water", "x_ethanol", "x_isobutanol", "x_nacl", "source"]
    write_csv_rows(fig_dir / "data" / "without_nacl_29315_digitized.csv", fieldnames, _phase_rows_for_csv(no_salt_rows, "digitized_local_paper"))
    write_csv_rows(fig_dir / "data" / "with_5wt_nacl_29315.csv", fieldnames, _phase_rows_for_csv(five_rows, "paper_table"))
    write_csv_rows(fig_dir / "data" / "with_10wt_nacl_29315.csv", fieldnames, _phase_rows_for_csv(ten_rows, "paper_table"))

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
    for salt_wt, temperature_k in sorted(EXPERIMENTAL_CASES):
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
    write_csv_rows(fig_dir / "data" / f"figure_{figure_number}_metrics.csv", ["temperature_K", "salt_wtfrac", "x_ethanol_aq", "distribution", "separation", "source"], rows)

    fig, ax = plt.subplots(figsize=(6.0, 4.4))
    fig.subplots_adjust(left=0.13, right=0.98, top=0.97, bottom=0.24)
    for row in rows:
        ax.scatter(row["x_ethanol_aq"], row[y_key], color=_metric_color(row["temperature_K"]), marker=_metric_marker(row["salt_wtfrac"]), s=22, linewidths=0.6)
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


def _aad_summary(exp_rows: list[dict], model_rows: list[dict]) -> dict:
    valid = [row for row in model_rows if np.all(np.isfinite(row["organic_formula"])) and np.all(np.isfinite(row["aqueous_formula"]))]
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
    temps = sorted({key[1] for key in EXPERIMENTAL_CASES if abs(key[0] - salt_wt) < 1e-9})
    figure_number_map = {(0.05, 293.15): 2, (0.05, 303.15): 3, (0.05, 313.15): 4, (0.10, 293.15): 5, (0.10, 303.15): 6, (0.10, 313.15): 7}
    force_recompute = os.environ.get("KHUDAIDA_FORCE_RECOMPUTE", "").strip().lower() in {"1", "true", "yes", "on"}
    for temperature_k in temps:
        exp_rows = _experimental_rows(salt_wt, temperature_k)
        feed_rows = _digitized_feed_rows_for_figure(figure_number_map[(salt_wt, temperature_k)], temperature_k, salt_wt) or _derived_feed_rows(salt_wt, temperature_k)
        fig_dir = ROOT / f"figure_{figure_number_map[(salt_wt, temperature_k)]}"
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
    columns = ["T / K", "Model", "Org x1", "Org x2", "Org x3", "Org x4", "Aq x1", "Aq x2", "Aq x3", "Aq x4", "Grand AAD"]
    for table_number, salt_wt in ((9, 0.05), (10, 0.10)):
        rows = _table_rows_for_png(salt_wt)
        csv_rows = [dict(zip(columns, row)) for row in rows]
        write_csv_rows(fig_dir / "data" / f"table_{table_number}.csv", columns, csv_rows)
        fig, ax = plt.subplots(figsize=(13.0, 3.0 + 0.28 * len(rows)))
        ax.axis("off")
        table = ax.table(cellText=rows, colLabels=columns, loc="center", cellLoc="center")
        table.auto_set_font_size(False)
        table.set_fontsize(8)
        table.scale(1.0, 1.2)
        ax.set_title(f"Table {table_number}: AAD comparison for {int(round(salt_wt * 100))} wt% NaCl", fontsize=11, pad=12)
        save_figure(fig, fig_dir / f"table_{table_number}.png")
        plt.close(fig)


def write_provenance_notes() -> None:
    notes = """# 2026 Khudaida analysis provenance

- Tables 3 and 4 in the local Khudaida markdown/PDF are treated as the canonical experimental tie-line source for Figures 2-7 and for the salted points in Figures 8-9.
- `2026_Khudaida` uses the paper's Table 5 pure-component parameters, Table 7 binary interaction parameters, and an exact copy of `2025_Figiel/user_options.json` as requested.
- Figure 1 salt-free data and the no-salt points in Figures 8-9 were reconstructed from the local paper figures because the Zotero baseline source remained inaccessible in this session.
- The no-salt baseline is therefore marked as `digitized_local_paper` in the emitted CSV files.
- Tables 9 and 10 include package-generated ePC-SAFT AAD values and paper-copied eNRTL/ePC-SAFT reference values for comparison.
- The legacy package solver note is retained here only as historical context; the multiphase LLE workflow is removed from the active Python package and will be rewritten later in native code.
"""
    (ROOT / "provenance_notes.md").write_text(notes, encoding="utf-8")


