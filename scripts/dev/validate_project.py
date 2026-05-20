from __future__ import annotations

import argparse
import os
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = REPO_ROOT / "src"
for import_root in (REPO_ROOT, SRC_ROOT):
    import_path = str(import_root)
    if import_path not in sys.path:
        sys.path.insert(0, import_path)

from scripts.dev.native_runtime_env import apply_native_runtime_env

apply_native_runtime_env(os.environ)

from epcsaft.capability_evidence import VALIDATION_LANES, validation_lane_commands

CHECK_COMMANDS: dict[str, tuple[tuple[str, ...], ...]] = {
    name: validation_lane_commands(name) for name in VALIDATION_LANES
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
