from __future__ import annotations

from pathlib import Path
import shutil


REFERENCE_INPUTS = (("LiAc-NaAc-KAc.csv", Path("data/reference/osmotic/water/LiAc-NaAc-KAc.csv")),)


def main() -> None:
    figure_root = Path(__file__).resolve().parents[1]
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

