from __future__ import annotations

import os
import subprocess
import sys
import tempfile
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]


def _run_command(cmd: list[str], env: dict[str, str]) -> None:
    print("Running:", " ".join(cmd), flush=True)
    subprocess.run(cmd, cwd=str(REPO_ROOT), check=True, env=env)


def _pip_temp_env() -> tuple[dict[str, str], tempfile.TemporaryDirectory[str], tempfile.TemporaryDirectory[str]]:
    pip_temp = tempfile.TemporaryDirectory(prefix="epcsaft-pip-")
    pip_cache = tempfile.TemporaryDirectory(prefix="epcsaft-pip-cache-")
    env = os.environ.copy()
    env["TMP"] = pip_temp.name
    env["TEMP"] = pip_temp.name
    env["TMPDIR"] = pip_temp.name
    env["PIP_CACHE_DIR"] = pip_cache.name
    env.setdefault("PIP_DISABLE_PIP_VERSION_CHECK", "1")
    return env, pip_temp, pip_cache


def main() -> int:
    pyproject = REPO_ROOT / "pyproject.toml"
    if not pyproject.exists():
        print(f"error: expected {pyproject} to exist")
        return 1

    env, pip_temp, pip_cache = _pip_temp_env()
    with pip_temp, pip_cache:
        _run_command(
            [
                sys.executable,
                "-m",
                "pip",
                "install",
                "-e",
                ".",
                "--no-build-isolation",
                "--config-settings",
                "editable_mode=compat",
            ],
            env,
        )
        _run_command([sys.executable, "scripts/build_epcsaft.py", "--force"], env)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
