import argparse
import os
import subprocess
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


def _run_install_dev(repo_root: Path, env: dict[str, str]) -> None:
    install_script = repo_root / "scripts" / "install_dev.py"
    print("Running:", f"{sys.executable} {install_script}", flush=True)
    subprocess.run([sys.executable, str(install_script)], cwd=str(repo_root), check=True, env=env)


def _run_build(repo_root: Path, env: dict[str, str], *, force: bool) -> None:
    build_script = repo_root / "scripts" / "build_epcsaft.py"
    if not build_script.exists():
        print(f"warning: build script not found at {build_script}")
        return

    build_cmd = [sys.executable, str(build_script)]
    if force:
        build_cmd.append("--force")
    print("Running:", " ".join(build_cmd), flush=True)
    subprocess.run(build_cmd, cwd=str(repo_root), check=True, env=env)


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
    parser.add_argument(
        "--skip-build",
        action="store_true",
        help="Run pytest without checking whether the editable install needs a rebuild.",
    )
    parser.add_argument(
        "--force-build",
        action="store_true",
        help="Force a native extension rebuild before running pytest.",
    )
    parser.add_argument(
        "--reinstall-editable",
        action="store_true",
        help="Repair/recreate the editable install before running pytest.",
    )
    args, pytest_args = parser.parse_known_args()

    repo_root = _repo_root()
    pytest_temp = _pytest_temp(repo_root)
    env = _pytest_env(pytest_temp)

    if args.reinstall_editable:
        _run_install_dev(repo_root, env)
    elif not args.skip_build:
        _run_build(repo_root, env, force=args.force_build)

    cmd = _pytest_args(pytest_args, pytest_temp)
    print("Running:", f"{sys.executable} -m pytest", " ".join(cmd), flush=True)
    os.environ.update(env)

    _patch_windows_pytest_temp_acl()
    _patch_pytest_cleanup()

    import pytest

    return int(pytest.main(cmd))


if __name__ == "__main__":
    raise SystemExit(main())

