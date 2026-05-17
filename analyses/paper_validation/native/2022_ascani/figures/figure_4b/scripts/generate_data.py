from __future__ import annotations

import sys
from pathlib import Path

for _candidate in Path(__file__).resolve().parents:
    if (_candidate / "scripts" / "plot_outputs.py").is_file():
        if str(_candidate) not in sys.path:
            sys.path.insert(0, str(_candidate))
        break
else:
    raise ModuleNotFoundError("Could not locate repo root containing scripts/plot_outputs.py")

from scripts.plot_outputs import analysis_scripts_dir, figure_output_path

ANALYSIS_SCRIPTS = analysis_scripts_dir(__file__)
if str(ANALYSIS_SCRIPTS) not in sys.path:
    sys.path.insert(0, str(ANALYSIS_SCRIPTS))

import _common as common


def main() -> None:
    rows = common.phase_diagram_rows()
    common.write_rows(
        figure_output_path(__file__, "data/figure_4b_phase_diagram.csv"),
        ["series", "phase", "w_water", "w_butanol", "w_total_salt", "x_plot", "y_plot"],
        rows,
    )

    _mix, result = common.current_solution()
    common.write_rows(
        figure_output_path(__file__, "data/current_phase_compositions.csv"),
        [
            "source",
            "phase",
            "component",
            "formula_mole_fraction",
            "explicit_mole_fraction",
            "phase_fraction",
            "density_mol_m3",
        ],
        common.current_phase_formula_rows(result),
    )
    common.write_rows(
        figure_output_path(__file__, "data/paper_case2_phase_compositions.csv"),
        ["source", "phase", "component", "formula_mole_fraction"],
        common.paper_phase_formula_rows(),
    )
    common.write_rows(
        figure_output_path(__file__, "data/current_feed_formula_composition.csv"),
        ["source", "phase", "component", "formula_mole_fraction"],
        common.current_feed_formula_rows(),
    )
    print(f"[write] {figure_output_path(__file__, 'data/figure_4b_phase_diagram.csv')}")


if __name__ == "__main__":
    main()
