from __future__ import annotations

import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def _run(path: Path) -> None:
    print(f"[run] {path}")
    subprocess.run([sys.executable, str(path)], check=True)


def main() -> None:
    figure_ids = (
        "baygi_outputs",
        "api_parity",
        "contribution_outputs",
        "equilibrium_confidence",
        "equilibrium_outputs",
        "native_outputs",
        "property_outputs",
        "regression_outputs",
    )
    for figure_id in figure_ids:
        figure_scripts = ROOT.parent / "figures" / figure_id / "scripts"
        _run(figure_scripts / "generate_data.py")
        _run(figure_scripts / f"plot_{figure_id}.py")


if __name__ == "__main__":
    main()
