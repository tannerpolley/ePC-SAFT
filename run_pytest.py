import argparse
import os
import subprocess
import sys
import uuid
from pathlib import Path


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
        help="Force an editable reinstall before running pytest.",
    )
    args, pytest_args = parser.parse_known_args()

    repo_root = Path(__file__).resolve().parent
    pytest_temp = repo_root / "build" / "pytest-temp" / f"run-{os.getpid()}-{uuid.uuid4().hex[:8]}"
    pytest_temp.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    env["TMP"] = str(pytest_temp.resolve())
    env["TEMP"] = str(pytest_temp.resolve())
    env["TMPDIR"] = str(pytest_temp.resolve())

    build_script = repo_root / "scripts" / "build_epcsaft.py"
    if build_script.exists() and not args.skip_build:
        build_cmd = [sys.executable, str(build_script)]
        if args.force_build:
            build_cmd.append("--force")
        print("Running:", " ".join(build_cmd), flush=True)
        subprocess.run(build_cmd, cwd=str(repo_root), check=True, env=env)
    elif not build_script.exists():
        print(f"warning: build script not found at {build_script}")

    has_positional_target = any(not arg.startswith("-") for arg in pytest_args)
    cmd = []
    if has_positional_target:
        cmd.extend(pytest_args)
    else:
        cmd.append("tests")
        cmd.extend(pytest_args)
    if not any(arg == "--basetemp" or arg.startswith("--basetemp=") for arg in pytest_args):
        cmd.extend(["--basetemp", str(pytest_temp)])
    print("Running:", f"{sys.executable} -m pytest", " ".join(cmd), flush=True)
    os.environ.update(env)

    if os.name == "nt":
        original_mkdir = Path.mkdir

        def sandbox_safe_mkdir(self, mode=0o777, parents=False, exist_ok=False):
            return original_mkdir(self, mode=0o777, parents=parents, exist_ok=exist_ok)

        Path.mkdir = sandbox_safe_mkdir

    # Pytest's Windows cleanup hook can trip over Codex sandbox ACLs after
    # tmp_path tests pass. Keep the per-run basetemp, but skip that hook.
    try:
        import _pytest.pathlib as pytest_pathlib
        import _pytest.tmpdir as pytest_tmpdir

        pytest_pathlib.cleanup_dead_symlinks = lambda root: None
        pytest_tmpdir.cleanup_dead_symlinks = lambda root: None
    except Exception:
        pass

    import pytest

    return int(pytest.main(cmd))


if __name__ == "__main__":
    raise SystemExit(main())

