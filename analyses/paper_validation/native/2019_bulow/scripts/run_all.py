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
    workflows = (
        ("figure_1", "plot_figure_1.py"),
        ("figure_2", "plot_figure_2.py"),
        ("figure_3", "plot_figure_3.py"),
        ("figure_4", "plot_figure_4.py"),
        ("figure_5", "plot_figure_5.py"),
        ("figure_6a", "plot_figure_6a.py"),
        ("figure_6b", "plot_figure_6b.py"),
        ("figure_7", "plot_figure_7.py"),
    )
    for figure_id, plot_script in workflows:
        figure_scripts = ROOT.parent / "figures" / figure_id / "scripts"
        run_script(figure_scripts / "generate_data.py")
        run_script(figure_scripts / plot_script)
    print("[done] 2019 figure scripts completed.")
    print(f"[note] Data-gap report: {ROOT / 'data_gap_report.md'}")


if __name__ == "__main__":
    main()
