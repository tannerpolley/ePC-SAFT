from __future__ import annotations

import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def run_script(path: Path) -> None:
    if not path.exists():
        print(f"[skip] missing script: {path}")
        return
    cmd = [sys.executable, str(path)]
    print(f"[run] {' '.join(cmd)}")
    subprocess.run(cmd, check=True)


def main() -> None:
    figure_ids = ("figure_2d", "figure_3", "figure_5", "figure_6", "figure_7")
    for figure_id in figure_ids:
        figure_scripts = ROOT.parent / "figures" / figure_id / "scripts"
        run_script(figure_scripts / "generate_data.py")
        run_script(figure_scripts / f"plot_{figure_id}.py")
    print("[done] 2012 figure scripts completed.")
    print(f"[note] Data-gap report: {ROOT / 'data_gap_report.md'}")


if __name__ == "__main__":
    main()
