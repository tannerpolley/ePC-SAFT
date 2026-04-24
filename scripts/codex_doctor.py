from __future__ import annotations

import importlib.util
import os
import shutil
import subprocess
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = REPO_ROOT / "src"


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


def _module_path(module_name: str) -> tuple[Path | None, str | None]:
    try:
        spec = importlib.util.find_spec(module_name)
    except Exception as exc:
        return None, str(exc)
    if spec is None or spec.origin is None:
        return None, None
    return Path(spec.origin).resolve(), None


def _tool_path(name: str) -> str:
    return shutil.which(name) or "<missing>"


def main() -> int:
    if str(SRC_ROOT) not in sys.path:
        sys.path.insert(0, str(SRC_ROOT))

    branch = _git_output("branch", "--show-current") or "<unknown>"
    head = _git_output("rev-parse", "--short", "HEAD") or "<unknown>"
    package_path, package_error = _module_path("epcsaft")
    core_path, core_error = _module_path("epcsaft._core")

    print(f"repo_root: {REPO_ROOT}")
    print(f"python: {sys.executable}")
    print(f"python_prefix: {sys.prefix}")
    print(f"conda_env_var: {os.environ.get('CONDA_DEFAULT_ENV', '<unset>')}")
    print(f"git_branch: {branch}")
    print(f"git_head: {head}")
    print(f"uv: {_tool_path('uv')}")
    print(f"cmake: {_tool_path('cmake')}")
    print(f"ninja: {_tool_path('ninja')}")
    print(f"epcsaft_import: {package_path if package_path else '<missing>'}")
    print(f"epcsaft_import_error: {package_error or '<none>'}")
    print(f"epcsaft_core: {core_path if core_path else '<missing>'}")
    print(f"epcsaft_core_error: {core_error or '<none>'}")

    if package_path is None:
        print("install_state: missing-package")
        print("next_command: uv sync --no-install-project")
        return 1
    if core_path is None:
        print("install_state: missing-core")
        print("next_command: uv run python scripts\\build_epcsaft.py --clean")
        return 1

    print("install_state: current")
    print("next_command: none")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
