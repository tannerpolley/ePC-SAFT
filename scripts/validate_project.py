from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]


CHECK_COMMANDS: dict[str, tuple[tuple[str, ...], ...]] = {
    "quick": (
        ("scripts/doctor.py",),
        ("run_pytest.py", "-q"),
    ),
    "confidence": (
        ("scripts/doctor.py",),
        ("run_pytest.py", "--confidence", "-q"),
    ),
    "docs": (("-m", "sphinx", "-b", "html", "docs", "build/docs-html"),),
    "full": (
        ("scripts/doctor.py",),
        ("run_pytest.py", "--all", "-q"),
    ),
    "ceres-cppad": (
        ("scripts/build_epcsaft.py", "--enable-ceres", "--enable-cppad"),
        (
            "run_pytest.py",
            "tests/native/test_ceres_pure_regression.py",
            "tests/native/test_ceres_binary_regression.py",
            "tests/regression/test_literature_binary_kij_regression.py",
            "-q",
        ),
    ),
}


def _python_command(args: tuple[str, ...]) -> list[str]:
    return [sys.executable, *args]


def run_mode(mode: str) -> int:
    for args in CHECK_COMMANDS[mode]:
        cmd = _python_command(args)
        print("Running:", " ".join(cmd), flush=True)
        completed = subprocess.run(cmd, cwd=REPO_ROOT, check=False)
        if completed.returncode != 0:
            return int(completed.returncode)
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Run the repo's standard validation modes.",
    )
    parser.add_argument(
        "mode",
        choices=tuple(CHECK_COMMANDS),
        help="Validation bundle to run.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    return run_mode(args.mode)


if __name__ == "__main__":
    raise SystemExit(main())
