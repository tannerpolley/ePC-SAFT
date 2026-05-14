from __future__ import annotations

import sys


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
from scripts.plot_outputs import REPO_ROOT

SCRIPT_DIR = Path(__file__).resolve().parent
ANALYSIS_ROOT = SCRIPT_DIR.parents[2]
DIAGNOSTICS_DIR = SCRIPT_DIR / "diagnostics"

if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
if str(DIAGNOSTICS_DIR) not in sys.path:
    sys.path.insert(0, str(DIAGNOSTICS_DIR))

from scripts.plot_outputs import paper_validation_path
from figure6b_libr_ethanol_contributions import run_analysis


def main() -> None:
    run_analysis(
        data_path=ANALYSIS_ROOT / "figures" / "figure_6b" / "input" / "ethanol-LiBr.csv",
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

