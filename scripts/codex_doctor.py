from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

import build_pcsaft


REPO_ROOT = Path(__file__).resolve().parents[1]


def _git_output(*args: str) -> str | None:
    try:
        result = subprocess.run(
            ["git", "-C", str(REPO_ROOT), *args],
            check=True,
            capture_output=True,
            text=True,
        )
    except (OSError, subprocess.CalledProcessError):
        return None
    return result.stdout.strip()


def main() -> int:
    branch = _git_output("branch", "--show-current") or "<unknown>"
    head = _git_output("rev-parse", "--short", "HEAD") or "<unknown>"
    install_state, install_reason = build_pcsaft._editable_install_state()
    rebuild_plan = build_pcsaft._rebuild_plan()
    imported = build_pcsaft._installed_module_path()

    print(f"repo_root: {REPO_ROOT}")
    print(f"python: {sys.executable}")
    print(f"conda_env: {os.environ.get('CONDA_DEFAULT_ENV', '<unset>')}")
    print(f"git_branch: {branch}")
    print(f"git_head: {head}")
    print(f"pcsaft_import: {imported if imported is not None else '<missing>'}")
    print(f"install_state: {install_state}")
    print(f"install_reason: {install_reason}")

    if rebuild_plan is None:
        print("rebuild_plan: none")
        return 0

    action, reason = rebuild_plan
    print(f"rebuild_plan: {action}")
    print(f"rebuild_reason: {reason}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
