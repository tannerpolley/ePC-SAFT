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
        ("figure_2", "plot_figure_2.py"),
        ("figure_3", "plot_figure_3.py"),
        ("figure_2_regressed", "plot_figure_2_regressed.py"),
        ("figure_3_regressed", "plot_figure_3_regressed.py"),
    )
    for figure_id, plot_script in workflows:
        figure_scripts = ROOT.parent / "figures" / figure_id / "scripts"
        run_script(figure_scripts / "generate_data.py")
        run_script(figure_scripts / plot_script)
    print("[done] 2015 Baygi figure scripts completed.")


if __name__ == "__main__":
    main()
