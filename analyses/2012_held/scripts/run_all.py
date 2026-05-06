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
    run_script(ROOT / "figure_2d" / "plot_figure_2d.py")
    run_script(ROOT / "figure_3" / "plot_figure_3.py")
    run_script(ROOT / "figure_5" / "plot_figure_5.py")
    run_script(ROOT / "figure_6" / "plot_figure_6.py")
    run_script(ROOT / "figure_7" / "plot_figure_7.py")
    print("[done] 2012 figure scripts completed.")
    print(f"[note] Data-gap report: {ROOT / 'data_gap_report.md'}")


if __name__ == "__main__":
    main()
