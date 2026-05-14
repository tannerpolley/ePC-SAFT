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
from scripts.plot_outputs import analysis_root

ANALYSIS_ROOT = analysis_root(__file__)
if str(ANALYSIS_ROOT) not in sys.path:
    sys.path.insert(0, str(ANALYSIS_ROOT))

from shared import figure_data


def main() -> None:
    path = figure_data.write_figure("figure_6")
    print(path.relative_to(figure_data.REPO_ROOT))


if __name__ == "__main__":
    main()
