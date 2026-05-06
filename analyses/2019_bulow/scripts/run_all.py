from __future__ import annotations

import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def run_script(path: Path) -> None:
    cmd = [sys.executable, str(path)]
    print(f"[run] {' '.join(cmd)}")
    subprocess.run(cmd, check=True)


def main() -> None:
    run_script(ROOT / "figure_1" / "plot_figure_1.py")
    run_script(ROOT / "figure_2" / "plot_figure_2.py")
    run_script(ROOT / "figure_3" / "plot_figure_3.py")
    run_script(ROOT / "figure_4" / "plot_figure_4.py")
    run_script(ROOT / "figure_5" / "plot_figure_5.py")
    run_script(ROOT / "figure_6" / "plot_figure_6a.py")
    run_script(ROOT / "figure_6" / "plot_figure_6b.py")
    run_script(ROOT / "figure_7" / "plot_figure_7.py")
    print("[done] 2019 figure scripts completed.")
    print(f"[note] Data-gap report: {ROOT / 'data_gap_report.md'}")


if __name__ == "__main__":
    main()
