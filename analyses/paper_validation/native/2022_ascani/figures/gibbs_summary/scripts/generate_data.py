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
    comparison, phases = common.gibbs_rows()
    common.write_rows(
        figure_output_path(__file__, "data/gibbs_comparison.csv"),
        ["quantity", "paper", "current_native_ln_fugacity_basis_j_per_mol", "current_minus_paper"],
        comparison,
    )
    common.write_rows(
        figure_output_path(__file__, "data/current_phase_gibbs.csv"),
        ["phase", "phase_fraction", "g_hat_phase_j_per_mol", "density_mol_m3"],
        phases,
    )
    print(f"[write] {figure_output_path(__file__, 'data/gibbs_comparison.csv')}")


if __name__ == "__main__":
    main()
