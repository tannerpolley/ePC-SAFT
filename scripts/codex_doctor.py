from __future__ import annotations

import os
import importlib.util
import subprocess
import sys
from pathlib import Path

import build_epcsaft


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


def _module_path(module_name: str) -> Path | None:
    spec = importlib.util.find_spec(module_name)
    if spec is None or spec.origin is None:
        return None
    return Path(spec.origin).resolve()


def _next_command(rebuild_plan: tuple[str, str] | None) -> str:
    if rebuild_plan is None:
        return "none"
    action, _ = rebuild_plan
    if action == "install-dev":
        return "python scripts/install_dev.py"
    if action == "build-ext-inplace":
        return "python scripts/build_epcsaft.py"
    return "inspect scripts/codex_doctor.py output"


def _python_env_name() -> str:
    parts = Path(sys.executable).resolve().parts
    for idx, part in enumerate(parts[:-1]):
        if part.lower() == "envs":
            return parts[idx + 1]
    return Path(sys.prefix).resolve().name or "<unknown>"


def main() -> int:
    branch = _git_output("branch", "--show-current") or "<unknown>"
    head = _git_output("rev-parse", "--short", "HEAD") or "<unknown>"
    install_state, install_reason = build_epcsaft._editable_install_state()
    rebuild_plan = build_epcsaft._rebuild_plan()
    imported = build_epcsaft._installed_module_path()
    extension = _module_path("epcsaft.epcsaft")

    print(f"repo_root: {REPO_ROOT}")
    print(f"python: {sys.executable}")
    print(f"python_env: {_python_env_name()}")
    print(f"conda_env_var: {os.environ.get('CONDA_DEFAULT_ENV', '<unset>')}")
    print(f"git_branch: {branch}")
    print(f"git_head: {head}")
    print(f"epcsaft_import: {imported if imported is not None else '<missing>'}")
    print(f"epcsaft_extension: {extension if extension is not None else '<missing>'}")
    print(f"install_state: {install_state}")
    print(f"install_reason: {install_reason}")

    if rebuild_plan is None:
        print("rebuild_plan: none")
        print("next_command: none")
        return 0

    action, reason = rebuild_plan
    print(f"rebuild_plan: {action}")
    print(f"rebuild_reason: {reason}")
    print(f"next_command: {_next_command(rebuild_plan)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

