from __future__ import annotations

import csv
from pathlib import Path
import sys

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import _common as common

OUTPUT = Path(__file__).with_name("figure_2.png")
DIGITIZED_SOURCES = (
    common.REPO_ROOT / "data" / "osmotic" / "water" / "2014_Held_figure_2_digitized.csv",
    common.analysis_data_path(__file__, "figure_2_digitized.csv", kind="input"),
)


def _load_optional_digitized() -> dict[tuple[str, float], tuple[np.ndarray, np.ndarray]]:
    out: dict[tuple[str, float], tuple[np.ndarray, np.ndarray]] = {}
    data_path = next((path for path in DIGITIZED_SOURCES if path.exists()), None)
    if data_path is None:
        return out

    grouped: dict[tuple[str, float], list[tuple[float, float]]] = {}
    with data_path.open("r", newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            salt = str(row.get("salt", "")).strip()
            if salt not in {"NaCl", "KCl"}:
                continue
            try:
                temp_k = float(row.get("temperature_K", ""))
                m = float(row.get("molality", ""))
                phi = float(row.get("osmotic", ""))
            except (TypeError, ValueError):
                continue
            key = (salt, temp_k)
            grouped.setdefault(key, []).append((m, phi))

    for key, pts in grouped.items():
        pts_sorted = sorted(pts, key=lambda p: p[0])
        m = np.asarray([p[0] for p in pts_sorted], dtype=float)
        phi = np.asarray([p[1] for p in pts_sorted], dtype=float)
        out[key] = (m, phi)
    return out


def main() -> None:
    common.configure_style()

    m_plot = np.linspace(0.01, 4.0, 401)
    digitized = _load_optional_digitized()

    fig, ax = plt.subplots(figsize=(7.9, 4.9))

    temp_styles = {
        298.15: {"line_color": "0.55", "marker_face": "none", "label_suffix": "298.15 K"},
        273.15: {"line_color": "black", "marker_face": "black", "label_suffix": "273.15 K"},
    }
    marker_map = {"NaCl": "o", "KCl": "s"}
    ls_map = {"NaCl": "-", "KCl": "--"}

    for temp_k in (298.15, 273.15):
        style = temp_styles[temp_k]
        for salt in ("NaCl", "KCl"):
            phi_model = common.calc_osmotic_curve(salt, m_plot, "2014", T=temp_k)
            ax.plot(
                m_plot,
                phi_model,
                color=style["line_color"],
                linewidth=1.9,
                linestyle=ls_map[salt],
                zorder=3,
            )

            data_key = (salt, temp_k)
            if data_key in digitized:
                m_exp, phi_exp = digitized[data_key]
                ax.scatter(
                    m_exp,
                    phi_exp,
                    marker=marker_map[salt],
                    s=34,
                    facecolors=style["marker_face"],
                    edgecolors="black",
                    linewidths=0.95,
                    zorder=6,
                )

    # Exact NaCl 298.15 K points are available in machine-readable CSV.
    m_nacl, phi_nacl = common.load_osmotic_data("NaCl")
    ax.scatter(
        m_nacl,
        phi_nacl,
        marker="o",
        s=34,
        facecolors="none",
        edgecolors="black",
        linewidths=1.0,
        zorder=7,
    )

    ax.set_xlim(0.0, 4.0)
    ax.set_ylim(0.84, 1.16)
    ax.set_xlabel(r"molality, $m$ / mol kg$^{-1}$")
    ax.set_ylabel(r"molal osmotic coefficient, $\phi_m$")
    ax.set_title("Held 2014 Fig. 2-style reproduction (NaCl and KCl)")
    ax.grid(True, alpha=0.25)

    handles_298 = [
        Line2D([0], [0], color=temp_styles[298.15]["line_color"], linewidth=1.9, linestyle="-", label="NaCl fit"),
        Line2D(
            [0],
            [0],
            linestyle="None",
            marker="o",
            markersize=5.8,
            markerfacecolor="none",
            markeredgecolor="black",
            markeredgewidth=1.0,
            label="NaCl data",
        ),
        Line2D([0], [0], color=temp_styles[298.15]["line_color"], linewidth=1.9, linestyle="--", label="KCl fit"),
        Line2D(
            [0],
            [0],
            linestyle="None",
            marker="s",
            markersize=5.8,
            markerfacecolor="none",
            markeredgecolor="black",
            markeredgewidth=1.0,
            label="KCl data",
        ),
    ]
    handles_273 = [
        Line2D([0], [0], color=temp_styles[273.15]["line_color"], linewidth=1.9, linestyle="-", label="NaCl fit"),
        Line2D(
            [0],
            [0],
            linestyle="None",
            marker="o",
            markersize=5.8,
            markerfacecolor="black",
            markeredgecolor="black",
            markeredgewidth=1.0,
            label="NaCl data",
        ),
        Line2D([0], [0], color=temp_styles[273.15]["line_color"], linewidth=1.9, linestyle="--", label="KCl fit"),
        Line2D(
            [0],
            [0],
            linestyle="None",
            marker="s",
            markersize=5.8,
            markerfacecolor="black",
            markeredgecolor="black",
            markeredgewidth=1.0,
            label="KCl data",
        ),
    ]

    legend_298 = ax.legend(
        handles=handles_298,
        title="298.15 K",
        fontsize=8,
        title_fontsize=8,
        loc="upper left",
        bbox_to_anchor=(0.02, 0.98),
        labelspacing=0.35,
        handlelength=2.2,
        borderaxespad=0.2,
    )
    ax.add_artist(legend_298)
    ax.legend(
        handles=handles_273,
        title="273.15 K",
        fontsize=8,
        title_fontsize=8,
        loc="upper left",
        bbox_to_anchor=(0.22, 0.98),
        labelspacing=0.35,
        handlelength=2.2,
        borderaxespad=0.2,
    )

    common.save_figure(fig, OUTPUT)
    plt.close(fig)


if __name__ == "__main__":
    main()
