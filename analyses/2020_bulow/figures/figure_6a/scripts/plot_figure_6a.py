from __future__ import annotations

import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
ANALYSIS_ROOT = SCRIPT_DIR.parents[2]
REPO_ROOT = ANALYSIS_ROOT.parents[1]
DIAGNOSTICS_DIR = ANALYSIS_ROOT / "figure_6b" / "diagnostics"

if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
if str(DIAGNOSTICS_DIR) not in sys.path:
    sys.path.insert(0, str(DIAGNOSTICS_DIR))

from scripts.plot_outputs import paper_validation_path

try:
    from figure6a_libr_ethanol_analysis import run_analysis
except ModuleNotFoundError:
    run_analysis = None


def main() -> None:
    output_path = paper_validation_path(__file__, "figure_6a.png")
    if run_analysis is None:
        if output_path.exists():
            print("[skip] figure_6a helper source is unavailable; " f"keeping existing output at {output_path}")
            return
        raise ModuleNotFoundError(
            "figure6a_libr_ethanol_analysis is unavailable and no existing figure_6a.png was found."
        )

    run_analysis(
        data_path=REPO_ROOT / "data" / "MIAC" / "ethanol" / "ethanol-LiBr.csv",
        output_path=output_path,
        x_min=0.0,
        x_max=0.2,
        y_min=0.0,
        y_max=4.0,
        grid_points=1201,
        max_molality=None,
    )


if __name__ == "__main__":
    main()

