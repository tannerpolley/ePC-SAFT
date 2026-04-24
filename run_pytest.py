import argparse
import os
import shutil
import sys
import uuid
from pathlib import Path


def _repo_root() -> Path:
    return Path(__file__).resolve().parent


def _pytest_temp(repo_root: Path) -> Path:
    path = repo_root / "build" / "pytest-temp" / f"run-{os.getpid()}-{uuid.uuid4().hex[:8]}"
    path.mkdir(parents=True, exist_ok=True)
    return path


def _pytest_env(pytest_temp: Path) -> dict[str, str]:
    env = os.environ.copy()
    env["TMP"] = str(pytest_temp.resolve())
    env["TEMP"] = str(pytest_temp.resolve())
    env["TMPDIR"] = str(pytest_temp.resolve())
    return env


def _pytest_args(pytest_args: list[str], pytest_temp: Path) -> list[str]:
    cmd: list[str] = []
    has_positional_target = any(not arg.startswith("-") for arg in pytest_args)
    if has_positional_target:
        cmd.extend(pytest_args)
    else:
        cmd.append("tests")
        cmd.extend(pytest_args)
    if not any(arg == "--basetemp" or arg.startswith("--basetemp=") for arg in pytest_args):
        cmd.extend(["--basetemp", str(pytest_temp)])
    return cmd


def _patch_windows_pytest_temp_acl() -> None:
    if os.name != "nt":
        return

    original_mkdir = Path.mkdir

    def sandbox_safe_mkdir(self, mode=0o777, parents=False, exist_ok=False):
        return original_mkdir(self, mode=0o777, parents=parents, exist_ok=exist_ok)

    Path.mkdir = sandbox_safe_mkdir


def _patch_pytest_cleanup() -> None:
    # Pytest's Windows cleanup hook can trip over Codex sandbox ACLs after
    # tmp_path tests pass. Keep the per-run basetemp, but skip that hook.
    try:
        import _pytest.pathlib as pytest_pathlib
        import _pytest.tmpdir as pytest_tmpdir

        pytest_pathlib.cleanup_dead_symlinks = lambda root: None
        pytest_tmpdir.cleanup_dead_symlinks = lambda root: None
    except Exception:
        pass


def main() -> int:
    parser = argparse.ArgumentParser()
    args, pytest_args = parser.parse_known_args()

    repo_root = _repo_root()
    pytest_temp = _pytest_temp(repo_root)
    env = _pytest_env(pytest_temp)
    src_root = repo_root / "src"
    sys.path.insert(0, str(src_root))
    env["PYTHONPATH"] = str(src_root)

    cmd = _pytest_args(pytest_args, pytest_temp)
    print("Running:", f"{sys.executable} -m pytest", " ".join(cmd), flush=True)
    os.environ.update(env)

    _patch_windows_pytest_temp_acl()
    _patch_pytest_cleanup()

    import pytest

    exit_code = int(pytest.main(cmd))
    if exit_code == 0:
        shutil.rmtree(pytest_temp, ignore_errors=True)
    else:
        print(f"Keeping pytest temp directory for failed run: {pytest_temp}", flush=True)
    return exit_code


if __name__ == "__main__":
    raise SystemExit(main())

