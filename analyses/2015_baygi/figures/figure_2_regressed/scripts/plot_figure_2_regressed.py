from __future__ import annotations

from pathlib import Path
import sys

import matplotlib.pyplot as plt

ANALYSIS_SCRIPTS = Path(__file__).resolve().parents[3] / "scripts"
if str(ANALYSIS_SCRIPTS) not in sys.path:
    sys.path.insert(0, str(ANALYSIS_SCRIPTS))

import _common as common

OUTPUT = Path(__file__).with_name("figure_2_regressed.png")


def _rho_kg_m3(row: dict) -> float:
    return float(row["rho_mol_m3"]) * common.MEA_MW


def main() -> None:
    common.configure_style()
    output_dir = common.output_path(OUTPUT).parent
    parameter_rows = common.regressed_parameter_rows(common.regressed_parameters_path())
    parameter_sets = common.parameter_sets_from_rows(parameter_rows)
    diagnostics_path = output_dir / "data" / "figure_2_regressed_diagnostics.csv"
    dippr, model_rows, diagnostics = common.saturation_rows_for_plot(diagnostics_path, parameter_sets=parameter_sets)

    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    ax.plot(
        [_rho_kg_m3(row) for row in dippr],
        [float(row["P_Pa"]) / 1000.0 for row in dippr],
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
            [_rho_kg_m3(row) for row in rows],
            [float(row["P_Pa"]) / 1000.0 for row in rows],
            color=style["color"],
            linestyle=style["linestyle"],
            linewidth=1.8,
            label=f"MEA {assoc_scheme}",
        )

    ax.set_xlim(870.0, 1030.0)
    ax.set_ylim(0.0, 120.0)
    ax.set_xticks([880, 920, 960, 1000])
    ax.set_yticks([0, 20, 40, 60, 80, 100, 120])
    ax.set_xlabel(r"saturated liquid density, $\rho^L_{sat}$ / kg m$^{-3}$")
    ax.set_ylabel(r"saturated vapor pressure, $p_{sat}$ / kPa")
    ax.set_title("Baygi 2015 Figure 2: ePC-SAFT-regressed pure MEA parameters")
    ax.grid(True, alpha=0.25)
    ax.legend(loc="best", fontsize=8)

    output_path = common.save_figure(fig, OUTPUT)
    plt.close(fig)

    common.write_csv_rows(
        diagnostics_path,
        ("series", "T_K", "P_Pa", "rho_mol_m3", "status", "contribution_terms"),
        [*dippr, *diagnostics],
    )
    metrics_path = output_path.parent / "data" / "figure_2_regressed_metrics.csv"
    common.write_csv_rows(
        metrics_path,
        ("series", "field", "n_solved", "aard_percent", "reported_table2_aad_percent", "within_reported_aad"),
        common.metric_rows([*dippr, *diagnostics], "P_Pa"),
    )


if __name__ == "__main__":
    main()

