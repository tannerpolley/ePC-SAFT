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
from scripts.plot_outputs import analysis_scripts_dir
import sys



import matplotlib.pyplot as plt

ANALYSIS_SCRIPTS = analysis_scripts_dir(__file__)
if str(ANALYSIS_SCRIPTS) not in sys.path:
    sys.path.insert(0, str(ANALYSIS_SCRIPTS))

import _common as common

OUTPUT = Path(__file__).with_name("figure_3_regressed.png")


def main() -> None:
    common.configure_style()
    output_dir = common.output_path(OUTPUT).parent
    parameter_rows = common.regressed_parameter_rows(common.regressed_parameters_path())
    parameter_sets = common.parameter_sets_from_rows(parameter_rows)
    diagnostics_path = output_dir / "data" / "figure_3_regressed_diagnostics.csv"
    dippr, model_rows, diagnostics = common.saturation_rows_for_plot(diagnostics_path, parameter_sets=parameter_sets)

    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    ax.plot(
        [float(row["T_K"]) for row in dippr],
        [float(row["rho_mol_m3"]) for row in dippr],
        color="black",
        linewidth=2.0,
        label="DIPPR",
    )
    for assoc_scheme, style in common.SCHEME_STYLES.items():
        rows = [row for row in model_rows if row["series"] == f"MEA {assoc_scheme}"]
        rows = sorted(rows, key=lambda row: float(row["T_K"]))
        if not rows:
            continue
        ax.plot(
            [float(row["T_K"]) for row in rows],
            [float(row["rho_mol_m3"]) for row in rows],
            color=style["color"],
            linestyle=style["linestyle"],
            linewidth=1.8,
            label=f"MEA {assoc_scheme}",
        )

    ax.set_xlabel(r"temperature, $T$ / K")
    ax.set_ylabel(r"saturated liquid density, $\rho^L_{sat}$ / mol m$^{-3}$")
    ax.set_title("Baygi 2015 Figure 3: ePC-SAFT-regressed pure MEA parameters")
    ax.grid(True, alpha=0.25)
    ax.legend(loc="best", fontsize=8)

    output_path = common.save_figure(fig, OUTPUT)
    plt.close(fig)

    common.write_csv_rows(
        diagnostics_path,
        ("series", "T_K", "P_Pa", "rho_mol_m3", "status", "contribution_terms"),
        [*dippr, *diagnostics],
    )
    metrics_path = output_path.parent / "data" / "figure_3_regressed_metrics.csv"
    common.write_csv_rows(
        metrics_path,
        ("series", "field", "n_solved", "aard_percent", "reported_table2_aad_percent", "within_reported_aad"),
        common.metric_rows([*dippr, *diagnostics], "rho_mol_m3"),
    )


if __name__ == "__main__":
    main()

