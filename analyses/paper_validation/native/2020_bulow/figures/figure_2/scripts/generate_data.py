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
import shutil


REFERENCE_INPUTS = (("water_comparisons.csv", Path("data/reference/G_solv/water_comparisons.csv")),)


def main() -> None:
    figure_root = figure_root_dir(__file__)
    repo_root = figure_root.parents[3]
    input_dir = figure_root / "input"
    output_dir = figure_root / "output"
    input_dir.mkdir(parents=True, exist_ok=True)
    output_dir.mkdir(parents=True, exist_ok=True)
    for filename, source_relpath in REFERENCE_INPUTS:
        source = repo_root / source_relpath
        dest = input_dir / filename
        shutil.copyfile(source, dest)
        print(f"[copy] {source_relpath} -> {dest.relative_to(repo_root)}")


if __name__ == "__main__":
    main()

