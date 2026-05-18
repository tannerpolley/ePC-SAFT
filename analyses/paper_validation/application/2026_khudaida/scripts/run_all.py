from __future__ import annotations

import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def _run(path: Path) -> None:
    print(f"[run] {path}")
    subprocess.run([sys.executable, str(path)], check=True)


def main() -> None:
    for idx in range(1, 10):
        figure_scripts = ROOT.parent / "figures" / f"figure_{idx}" / "scripts"
        _run(figure_scripts / "generate_data.py")
        _run(figure_scripts / f"plot_figure_{idx}.py")
    table_scripts = ROOT.parent / "figures" / "tables_9_10" / "scripts"
    _run(table_scripts / "generate_data.py")
    _run(table_scripts / "plot_tables_9_10.py")
    for idx in (2, 3):
        supporting_scripts = ROOT.parent / "figures" / f"figure_s{idx}" / "scripts"
        _run(supporting_scripts / "generate_data.py")
        _run(supporting_scripts / f"plot_figure_s{idx}.py")
    print("[done] 2026 Khudaida analysis complete.")


if __name__ == "__main__":
    main()
