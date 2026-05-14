from __future__ import annotations

from pathlib import Path
import sys


import sys as _bootstrap_sys
from pathlib import Path as _BootstrapPath

for _candidate in _BootstrapPath(__file__).resolve().parents:
    if (_candidate / "scripts" / "plot_outputs.py").is_file():
        if str(_candidate) not in _bootstrap_sys.path:
            _bootstrap_sys.path.insert(0, str(_candidate))
        break
else:
    raise ModuleNotFoundError("Could not locate repo root containing scripts/plot_outputs.py")
from scripts.plot_outputs import figure_root_dir


def main() -> None:
    figure_root = figure_root_dir(__file__)
    output_dir = figure_root / "output"
    output_dir.mkdir(parents=True, exist_ok=True)
    print(f"[skip] retained parameter exports live under {output_dir}.")


if __name__ == "__main__":
    main()
