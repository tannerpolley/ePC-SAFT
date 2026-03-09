from __future__ import annotations

from pathlib import Path
import sys

import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import _common as common

OUTPUT = Path(__file__).with_name("figure_7.png")
DATA = common.REPO_ROOT / "data" / "MIAC" / "methanol" / "methanol-NaBr.csv"


def main() -> None:
    common.configure_style()
    data = common.read_miac_dataset(DATA, "methanol")
    x_data = [row["molality"] for row in data]
    y_data = [row["miac_m"] for row in data]
    m_max = max(float(row["molality"]) for row in data)
    m_grid, y_default = common.mean_ionic_activity_curve("figiel_2025", "NaBr", "methanol", {"methanol": 1.0}, m_max, points=500)
    _, y_linear = common.mean_ionic_activity_curve(
        "figiel_2025",
        "NaBr",
        "methanol",
        {"methanol": 1.0},
        m_max,
        points=500,
        user_options={"elec_model": {"rel_perm": {"rule": 1}}},
    )

    fig, ax = plt.subplots(figsize=(3.4, 2.8))
    ax.scatter(x_data, y_data, s=24, facecolor=common.LIGHT_GRAY, edgecolor=common.GRAY_COLOR, linewidth=0.8, label="Literature")
    ax.plot(m_grid, y_default, color="black", linewidth=1.5, label="Figiel 2025")
    ax.plot(m_grid, y_linear, color="black", linewidth=1.3, linestyle="--", label="Rule 1")
    ax.set_xlim(0.0, m_max)
    ax.set_ylim(0.0, 1.0)
    ax.set_xlabel(r"$\bar{m}_{NaBr}$ / mol kg$^{-1}$")
    ax.set_ylabel(r"$\gamma_{\pm}^{m,*}$ / -")
    ax.legend(loc="upper right", fontsize=8)
    common.save_figure(fig, OUTPUT)
    plt.close(fig)


if __name__ == "__main__":
    main()
