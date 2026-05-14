from __future__ import annotations

import csv
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np

ANALYSIS_SCRIPTS = Path(__file__).resolve().parents[3] / "scripts"
if str(ANALYSIS_SCRIPTS) not in sys.path:
    sys.path.insert(0, str(ANALYSIS_SCRIPTS))

import _common as common

OUTPUT = Path(__file__).with_name("figure_4b.png")
DIGITIZED = common.analysis_data_path(__file__, "figure_4b_digitized.csv", kind="input")


def _load() -> dict[tuple[str, str], tuple[np.ndarray, np.ndarray]]:
    grouped: dict[tuple[str, str], list[tuple[float, float]]] = {}
    with DIGITIZED.open("r", newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            kind = str(row.get("kind", "")).strip().lower()
            solute = str(row.get("solute", "")).strip().lower()
            if kind not in {"exp", "strategy1", "strategy2"}:
                continue
            if solute not in {"glycine", "alanine"}:
                continue
            try:
                x = float(row["m_kcl"])
                y = float(row["solubility_molkg"])
            except (KeyError, TypeError, ValueError):
                continue
            grouped.setdefault((kind, solute), []).append((x, y))

    out: dict[tuple[str, str], tuple[np.ndarray, np.ndarray]] = {}
    for key, pts in grouped.items():
        pts_sorted = sorted(pts, key=lambda p: p[0])
        out[key] = (
            np.asarray([p[0] for p in pts_sorted], dtype=float),
            np.asarray([p[1] for p in pts_sorted], dtype=float),
        )
    return out


def _dense_xy(x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    xd = np.linspace(float(np.min(x)), float(np.max(x)), 250)
    yd = np.interp(xd, x, y)
    return xd, yd


def main() -> None:
    common.configure_style()
    data = _load()

    fig, ax = plt.subplots(figsize=(7.0, 5.0))

    x, y = data[("exp", "glycine")]
    ax.scatter(x, y, marker="s", s=40, color="black", label="Glycine data", zorder=6)
    x, y = data[("exp", "alanine")]
    ax.scatter(
        x,
        y,
        marker="o",
        s=54,
        facecolors="0.85",
        edgecolors="black",
        linewidths=0.9,
        label="Alanine data",
        zorder=6,
    )

    for kind, color, width, style in (
        ("strategy2", "black", 2.2, "-"),
        ("strategy1", "0.45", 1.2, "--"),
    ):
        x, y = data[(kind, "glycine")]
        xd, yd = _dense_xy(x, y)
        ax.plot(xd, yd, color=color, linewidth=width, linestyle=style, label=f"Glycine model ({kind})")

        x, y = data[(kind, "alanine")]
        xd, yd = _dense_xy(x, y)
        ax.plot(
            xd,
            yd,
            color="0.55" if kind == "strategy1" else "0.45",
            linewidth=width,
            linestyle=style,
            label=f"Alanine model ({kind})",
        )

    ax.set_xlim(0.0, 3.2)
    ax.set_ylim(1.5, 3.8)
    ax.set_xlabel(r"KCl molality, $m_{KCl}$ / mol kg$^{-1}$")
    ax.set_ylabel(r"amino-acid solubility / mol kg$^{-1}$")
    ax.set_title("Held 2014 Fig. 4b digitized reproduction")
    ax.grid(True, alpha=0.24)
    ax.legend(fontsize=8, ncol=2)

    common.save_figure(fig, OUTPUT)
    plt.close(fig)


if __name__ == "__main__":
    main()


