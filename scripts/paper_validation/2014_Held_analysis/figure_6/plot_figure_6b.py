from __future__ import annotations

import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

import _shared

import _common as common

OUTPUT = Path(__file__).with_name("figure_6b.png")
SQRT3_OVER_2 = math.sqrt(3.0) / 2.0


def _to_xy(phase_mass_fraction: np.ndarray) -> tuple[float, float]:
    w_h2o = float(phase_mass_fraction[_shared.IDX["H2O"]])
    w_nh4cl = float(phase_mass_fraction[_shared.IDX["NH4+"]]) + float(phase_mass_fraction[_shared.IDX["Cl-"]])
    w_buoh = float(phase_mass_fraction[_shared.IDX["Butanol"]])
    x = w_buoh + 0.5 * w_nh4cl
    y = SQRT3_OVER_2 * w_nh4cl
    return x, y


def _draw_grid(ax: plt.Axes) -> None:
    for frac in np.linspace(0.2, 0.8, 4):
        p1 = _to_xy(np.asarray([1.0 - frac, 0.0, frac * (_shared.MW[2] / _shared.MW_NH4CL), frac * (_shared.MW[3] / _shared.MW_NH4CL)], dtype=float))
        p2 = _to_xy(np.asarray([0.0, 1.0 - frac, frac * (_shared.MW[2] / _shared.MW_NH4CL), frac * (_shared.MW[3] / _shared.MW_NH4CL)], dtype=float))
        ax.plot([p1[0], p2[0]], [p1[1], p2[1]], color="0.82", linewidth=1.0, zorder=0)

        p1 = _to_xy(np.asarray([frac, 1.0 - frac, 0.0, 0.0], dtype=float))
        p2 = _to_xy(np.asarray([frac, 0.0, (1.0 - frac) * (_shared.MW[2] / _shared.MW_NH4CL), (1.0 - frac) * (_shared.MW[3] / _shared.MW_NH4CL)], dtype=float))
        ax.plot([p1[0], p2[0]], [p1[1], p2[1]], color="0.82", linewidth=1.0, zorder=0)

        p1 = _to_xy(np.asarray([1.0 - frac, frac, 0.0, 0.0], dtype=float))
        p2 = _to_xy(np.asarray([0.0, frac, (1.0 - frac) * (_shared.MW[2] / _shared.MW_NH4CL), (1.0 - frac) * (_shared.MW[3] / _shared.MW_NH4CL)], dtype=float))
        ax.plot([p1[0], p2[0]], [p1[1], p2[1]], color="0.82", linewidth=1.0, zorder=0)


def _plot_tie_lines(ax: plt.Axes, rows: tuple[dict, ...], line_color: str, marker_face: str, marker_edge: str, label_prefix: str, zbase: int) -> None:
    aqueous_pts = []
    organic_pts = []
    for idx, row in enumerate(rows):
        aq_xy = _to_xy(np.asarray(row["aqueous_mass_fraction"], dtype=float))
        org_xy = _to_xy(np.asarray(row["organic_mass_fraction"], dtype=float))
        aqueous_pts.append(aq_xy)
        organic_pts.append(org_xy)
        ax.plot(
            [aq_xy[0], org_xy[0]],
            [aq_xy[1], org_xy[1]],
            color=line_color,
            linewidth=1.0 if zbase < 5 else 1.2,
            zorder=zbase,
            label=f"{label_prefix} tie lines" if idx == 0 else None,
        )

    aqueous = np.asarray(aqueous_pts, dtype=float)
    organic = np.asarray(organic_pts, dtype=float)
    envelope = np.vstack([aqueous, organic[::-1], aqueous[:1]])
    ax.plot(
        envelope[:, 0],
        envelope[:, 1],
        color=line_color,
        linewidth=1.2 if zbase < 5 else 1.6,
        zorder=zbase + 1,
        label=f"{label_prefix} binodal",
    )
    ax.scatter(
        aqueous[:, 0],
        aqueous[:, 1],
        marker="o",
        s=42,
        facecolors=marker_face,
        edgecolors=marker_edge,
        linewidths=0.8,
        zorder=zbase + 2,
        label=None,
    )
    ax.scatter(
        organic[:, 0],
        organic[:, 1],
        marker="o",
        s=42,
        facecolors=marker_face,
        edgecolors=marker_edge,
        linewidths=0.8,
        zorder=zbase + 2,
        label=None,
    )


def main() -> None:
    common.configure_style()
    exp_rows = _shared.load_experimental_rows()
    model_rows = _shared.solve_model_rows()

    fig, ax = plt.subplots(figsize=(7.0, 6.2))

    tri = np.asarray(
        [
            _to_xy(np.asarray([1.0, 0.0, 0.0, 0.0], dtype=float)),
            _to_xy(np.asarray([0.0, 0.0, _shared.MW[2] / _shared.MW_NH4CL, _shared.MW[3] / _shared.MW_NH4CL], dtype=float)),
            _to_xy(np.asarray([0.0, 1.0, 0.0, 0.0], dtype=float)),
            _to_xy(np.asarray([1.0, 0.0, 0.0, 0.0], dtype=float)),
        ],
        dtype=float,
    )
    ax.plot(tri[:, 0], tri[:, 1], color="black", linewidth=1.4, zorder=2)
    _draw_grid(ax)

    _plot_tie_lines(ax, exp_rows, line_color="0.60", marker_face="white", marker_edge="0.25", label_prefix="Data", zbase=3)
    if model_rows:
        _plot_tie_lines(ax, model_rows, line_color="black", marker_face="0.75", marker_edge="0.10", label_prefix="Model", zbase=6)

    ax.text(-0.06, -0.03, r"$H_2O$", fontsize=12)
    ax.text(1.00, -0.03, "BuOH", fontsize=12, ha="left")
    ax.text(0.50, SQRT3_OVER_2 + 0.03, r"$NH_4Cl$", fontsize=12, ha="center")

    for frac in np.linspace(0.2, 0.8, 4):
        left = _to_xy(np.asarray([1.0 - frac, 0.0, frac * (_shared.MW[2] / _shared.MW_NH4CL), frac * (_shared.MW[3] / _shared.MW_NH4CL)], dtype=float))
        right = _to_xy(np.asarray([0.0, 1.0 - frac, frac * (_shared.MW[2] / _shared.MW_NH4CL), frac * (_shared.MW[3] / _shared.MW_NH4CL)], dtype=float))
        bottom = _to_xy(np.asarray([1.0 - frac, frac, 0.0, 0.0], dtype=float))
        ax.text(left[0] - 0.03, left[1], f"{1.0 - frac:.1f}", fontsize=8, ha="right", va="center")
        ax.text(right[0] + 0.03, right[1], f"{frac:.1f}", fontsize=8, ha="left", va="center")
        ax.text(bottom[0], -0.045, f"{frac:.1f}", fontsize=8, ha="center")

    ax.set_aspect("equal")
    ax.set_xlim(-0.08, 1.08)
    ax.set_ylim(-0.08, SQRT3_OVER_2 + 0.08)
    ax.axis("off")
    ax.set_title("Held 2014 Fig. 6b ternary reproduction")
    ax.legend(loc="lower center", fontsize=8, ncol=2)

    common.save_figure(fig, OUTPUT)
    plt.close(fig)


if __name__ == "__main__":
    main()
