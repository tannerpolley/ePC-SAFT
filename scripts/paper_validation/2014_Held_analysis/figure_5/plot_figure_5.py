from __future__ import annotations

import csv
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import _common as common

OUTPUT = Path(__file__).with_name("figure_5.png")
DIGITIZED = Path(__file__).with_name("figure_5_digitized.csv")


def _load() -> dict[tuple[str, str], tuple[np.ndarray, np.ndarray]]:
    grouped: dict[tuple[str, str], list[tuple[float, float]]] = {}
    with DIGITIZED.open("r", newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            kind = str(row.get("kind", "")).strip().lower()
            salt = str(row.get("salt", "")).strip()
            if kind not in {"exp", "model"} or salt not in {"NaCl", "NaI"}:
                continue
            try:
                m = float(row["m_salt"])
                y = float(row["m_benzene_mmolkg"])
            except (TypeError, ValueError, KeyError):
                continue
            grouped.setdefault((kind, salt), []).append((m, y))

    out: dict[tuple[str, str], tuple[np.ndarray, np.ndarray]] = {}
    for key, pts in grouped.items():
        pts_sorted = sorted(pts, key=lambda p: p[0])
        out[key] = (
            np.asarray([p[0] for p in pts_sorted], dtype=float),
            np.asarray([p[1] for p in pts_sorted], dtype=float),
        )
    return out


def _interp_dense(x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    xd = np.linspace(float(np.min(x)), float(np.max(x)), 250)
    yd = np.interp(xd, x, y)
    return xd, yd


def main() -> None:
    common.configure_style()
    data = _load()

    fig, ax = plt.subplots(figsize=(7.2, 4.9))

    # Experimental points
    x, y = data[("exp", "NaCl")]
    ax.scatter(x, y, marker="o", s=44, color="black", label="NaCl data", zorder=7)
    x, y = data[("exp", "NaI")]
    ax.scatter(x, y, marker="s", s=42, facecolors="0.65", edgecolors="0.25", linewidths=0.8, label="NaI data", zorder=7)

    # Digitized model curves from the same panel.
    x, y = data[("model", "NaCl")]
    xd, yd = _interp_dense(x, y)
    ax.plot(xd, yd, color="black", linewidth=2.0, label="NaCl model (digitized)")

    x, y = data[("model", "NaI")]
    xd, yd = _interp_dense(x, y)
    ax.plot(xd, yd, color="0.45", linewidth=1.9, label="NaI model (digitized)")

    ax.set_xlim(0.0, 3.4)
    ax.set_ylim(0.0, 24.0)
    ax.set_xlabel(r"salt molality, $m_{salt}$ / mol kg$^{-1}$")
    ax.set_ylabel(r"benzene molality, $m_{benzene}$ / mmol kg$^{-1}$")
    ax.set_title("Held 2014 Fig. 5 digitized reproduction")
    ax.grid(True, alpha=0.24)
    ax.legend(fontsize=8, ncol=2)

    common.save_figure(fig, OUTPUT)
    plt.close(fig)


if __name__ == "__main__":
    main()

