from __future__ import annotations

from pathlib import Path


def main() -> None:
    figure_root = Path(__file__).resolve().parents[1]
    output_dir = figure_root / "output"
    output_dir.mkdir(parents=True, exist_ok=True)
    print(f"[skip] retained parameter exports live under {output_dir}.")


if __name__ == "__main__":
    main()
