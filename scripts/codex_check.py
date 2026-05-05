from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]


CHECK_COMMANDS: dict[str, tuple[tuple[str, ...], ...]] = {
    "quick": (
        ("scripts/codex_doctor.py",),
        ("run_pytest.py", "--generic", "-q"),
    ),
    "confidence": (
        ("scripts/codex_doctor.py",),
        ("run_pytest.py", "--confidence", "-q"),
        ("scripts/build_plot_manifest.py", "--check"),
    ),
    "docs": (
        ("scripts/build_plot_manifest.py", "--check"),
        ("-m", "sphinx", "-b", "html", "docs", "build/docs-html"),
    ),
    "plots": (
        ("run_pytest.py", "--plots", "-q"),
        ("scripts/build_plot_manifest.py", "--refresh"),
        ("scripts/paper_validation/tools/build_analysis_galleries.py",),
    ),
    "full": (
        ("scripts/codex_doctor.py",),
        ("run_pytest.py", "--confidence", "-q"),
        ("run_pytest.py", "--equilibrium-confidence", "-q"),
        ("scripts/build_plot_manifest.py", "--check"),
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
        description="Run the repo's Codex-friendly validation modes.",
    )
    parser.add_argument(
        "mode",
        choices=tuple(CHECK_COMMANDS),
        help="Validation bundle to run: quick, confidence, docs, plots, or full.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    return run_mode(args.mode)


if __name__ == "__main__":
    raise SystemExit(main())
