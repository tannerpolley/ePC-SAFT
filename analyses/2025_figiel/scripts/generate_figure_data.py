from __future__ import annotations

import sys
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parent
ANALYSIS_ROOT = ROOT.parent
REPO_ROOT = ROOT.parents[2]
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
