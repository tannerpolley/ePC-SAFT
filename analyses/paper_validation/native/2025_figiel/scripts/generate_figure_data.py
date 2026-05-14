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
import subprocess

ROOT = Path(__file__).resolve().parent
ANALYSIS_ROOT = ROOT.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(ANALYSIS_ROOT) not in sys.path:
    sys.path.insert(0, str(ANALYSIS_ROOT))

from shared import figure_data

GENERATE_SCRIPTS = [ANALYSIS_ROOT / "figures" / figure_id / "scripts" / "generate_data.py" for figure_id in figure_data.GENERATORS]


def main() -> None:
    for script in GENERATE_SCRIPTS:
        subprocess.run([sys.executable, str(script.relative_to(REPO_ROOT))], cwd=REPO_ROOT, check=True)


if __name__ == "__main__":
    main()
