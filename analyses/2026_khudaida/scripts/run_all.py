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
        _run(ROOT / f"figure_{idx}" / f"plot_figure_{idx}.py")
    _run(ROOT / "tables_9_10" / "plot_tables_9_10.py")
    print("[done] 2026 Khudaida analysis complete.")


if __name__ == "__main__":
    main()
