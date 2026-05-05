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
    run_script(ROOT / "figure_2" / "plot_figure_2.py")
    run_script(ROOT / "figure_3" / "plot_figure_3.py")
    run_script(ROOT / "figure_2_regressed" / "plot_figure_2_regressed.py")
    run_script(ROOT / "figure_3_regressed" / "plot_figure_3_regressed.py")
    print("[done] 2015 Baygi figure scripts completed.")


if __name__ == "__main__":
    main()
