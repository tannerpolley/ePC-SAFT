from __future__ import annotations

from pathlib import Path


def main() -> None:
    figure_root = Path(__file__).resolve().parents[1]
    (figure_root / "input").mkdir(parents=True, exist_ok=True)
    (figure_root / "output").mkdir(parents=True, exist_ok=True)
    print(f"[skip] no standalone data-generation step for {figure_root.parent.name}/{figure_root.name}.")


if __name__ == "__main__":
    main()
