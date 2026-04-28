from __future__ import annotations

import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
ANALYSIS_ROOT = SCRIPT_DIR.parent
REPO_ROOT = ANALYSIS_ROOT.parents[3]
DIAGNOSTICS_DIR = SCRIPT_DIR / "diagnostics"

if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
if str(DIAGNOSTICS_DIR) not in sys.path:
    sys.path.insert(0, str(DIAGNOSTICS_DIR))

from scripts.plot_outputs import paper_validation_path
from figure6b_libr_ethanol_contributions import run_analysis


def main() -> None:
    run_analysis(
        data_path=REPO_ROOT / "data" / "MIAC" / "ethanol" / "ethanol-LiBr.csv",
        output_path=paper_validation_path(__file__, "figure_6b.png"),
        x_min=0.0,
        x_max=0.2,
        y_min=-3.0,
        y_max=4.0,
        grid_points=1201,
        max_molality=None,
        plot_title="LiBr in ethanol at 298.15 K and 1 bar",
        method="mu",
    )


if __name__ == "__main__":
    main()

