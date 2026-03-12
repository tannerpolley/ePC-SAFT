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
    run_script(ROOT / "figure_4" / "plot_figure_4.py")
    run_script(ROOT / "figure_5" / "plot_figure_5.py")
    run_script(ROOT / "figure_6" / "figure_6a" / "plot_figure_6a.py")
    run_script(ROOT / "figure_6" / "figure_6b" / "plot_figure_6b.py")
    print("[done] 2020 figure scripts completed.")
    print(f"[note] Diagnostics: {ROOT / 'diagnostics'}")


if __name__ == "__main__":
    main()
