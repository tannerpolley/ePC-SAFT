from __future__ import annotations

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
from scripts.plot_outputs import analysis_root
import sys



import matplotlib.pyplot as plt

ANALYSIS_ROOT = analysis_root(__file__)
if str(ANALYSIS_ROOT) not in sys.path:
    sys.path.insert(0, str(ANALYSIS_ROOT))

from shared import _common as common
from shared import figure_data

OUTPUT = Path(__file__).with_name("figure_7.png")


def main() -> None:
    common.configure_style()
    rows = figure_data.read_payload("figure_7")
    x_data, y_data = figure_data.xy(figure_data.select_rows(rows, series_id="data_NaBr_methanol"))
    m_grid, y_default = figure_data.xy(figure_data.select_rows(rows, series_id="model_default"))
    _, y_linear = figure_data.xy(figure_data.select_rows(rows, series_id="model_rule1"))
    m_max = max(x_data)

    fig, ax = plt.subplots(figsize=(3.4, 2.8))
    ax.scatter(
        x_data,
        y_data,
        s=24,
        facecolor=common.LIGHT_GRAY,
        edgecolor=common.GRAY_COLOR,
        linewidth=0.8,
        label="Literature",
    )
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
