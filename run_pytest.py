import argparse
import os
import shutil
import sys
import uuid
from pathlib import Path


GENERIC_TEST_TARGETS = (
    "tests/test_runtime.py",
    "tests/test_parameter_templates.py",
    "tests/test_equation_registry.py",
    "tests/test_regression.py",
    "tests/test_regression_api.py",
)
CONFIDENCE_TEST_TARGETS = GENERIC_TEST_TARGETS + ("tests/test_native_runtime_contracts.py",)
RUNTIME_TEST_TARGETS = ("tests/test_runtime.py", "tests/test_native_runtime_contracts.py")
API_TEST_TARGETS = ("tests/test_runtime.py", "tests/test_parameter_templates.py", "tests/test_regression_api.py")
NATIVE_TEST_TARGETS = ("tests/test_native_runtime_contracts.py",)
PROFILE_TEST_TARGETS = ("tests/test_runtime_profile.py",)


def _repo_root() -> Path:
    return Path(__file__).resolve().parent


def _pytest_temp(repo_root: Path) -> Path:
    configured_root = os.environ.get("EPCSAFT_PYTEST_TEMP_ROOT")
    if configured_root is not None:
        root = Path(configured_root).expanduser()
        if not root.is_absolute():
            root = (repo_root / root).resolve()
        root = root / "pytest-temp"
    else:
        root = repo_root / "build" / "pytest-temp"

    path = root / f"run-{os.getpid()}-{uuid.uuid4().hex[:8]}"
    path.mkdir(parents=True, exist_ok=True)
    return path


def _pytest_env(pytest_temp: Path, profile: bool = False) -> dict[str, str]:
    env = os.environ.copy()
    env["TMP"] = str(pytest_temp.resolve())
    env["TEMP"] = str(pytest_temp.resolve())
    env["TMPDIR"] = str(pytest_temp.resolve())
    if profile:
        env["EPCSAFT_RUN_PERF"] = "1"
        env["ePCSAFT_RUN_PERF"] = "1"
    return env


def _pytest_args(
    pytest_args: list[str],
    pytest_temp: Path,
    generic: bool = False,
    confidence: bool = False,
    runtime: bool = False,
    api: bool = False,
    native: bool = False,
    profile: bool = False,
) -> list[str]:
    cmd: list[str] = []
    has_predefined_targets = generic or confidence or runtime or api or native or profile
    if confidence:
        cmd.extend(CONFIDENCE_TEST_TARGETS)
    elif generic:
        cmd.extend(GENERIC_TEST_TARGETS)
    elif runtime:
        cmd.extend(RUNTIME_TEST_TARGETS)
    elif api:
        cmd.extend(API_TEST_TARGETS)
    elif native:
        cmd.extend(NATIVE_TEST_TARGETS)
    elif profile:
        cmd.extend(PROFILE_TEST_TARGETS)

    if has_predefined_targets:
        cmd.extend(pytest_args)
    else:
        has_positional_target = any(not arg.startswith("-") for arg in pytest_args)
        if has_positional_target:
            cmd.extend(pytest_args)
        else:
            cmd.append("tests")
            cmd.extend(pytest_args)

    if not any(arg == "--basetemp" or arg.startswith("--basetemp=") for arg in cmd):
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


def _failure_message(pytest_temp: Path) -> str:
    cleanup_path = str(pytest_temp.resolve())
    return (
        "Pytest failed; keeping temp directory for triage: "
        f"{cleanup_path}\n"
        "Cleanup with: "
        f"Remove-Item -Recurse -Force '{cleanup_path}' (PowerShell)\n"
        f"or rm -rf {cleanup_path} (POSIX shells)"
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    predefined = parser.add_mutually_exclusive_group()
    predefined.add_argument("--generic", action="store_true", help="Run the core generic test slice")
    predefined.add_argument(
        "--confidence",
        action="store_true",
        help="Run generic tests plus native runtime contract tests",
    )
    predefined.add_argument("--runtime", action="store_true", help="Run runtime API and native contract tests")
    predefined.add_argument("--api", action="store_true", help="Run public API and regression API tests")
    predefined.add_argument("--native", action="store_true", help="Run native runtime contract tests")
    predefined.add_argument("--profile", action="store_true", help="Run the opt-in runtime profile test")
    args, pytest_args = parser.parse_known_args()

    repo_root = _repo_root()
    pytest_temp = _pytest_temp(repo_root)
    env = _pytest_env(pytest_temp, profile=args.profile)
    src_root = repo_root / "src"
    sys.path.insert(0, str(src_root))
    env["PYTHONPATH"] = str(src_root)

    cmd = _pytest_args(
        pytest_args,
        pytest_temp,
        args.generic,
        confidence=args.confidence,
        runtime=args.runtime,
        api=args.api,
        native=args.native,
        profile=args.profile,
    )
    print("Running:", f"{sys.executable} -m pytest", " ".join(cmd), flush=True)
    os.environ.update(env)

    _patch_windows_pytest_temp_acl()
    _patch_pytest_cleanup()

    import pytest

    exit_code = int(pytest.main(cmd))
    if exit_code == 0:
        shutil.rmtree(pytest_temp, ignore_errors=True)
    else:
        print(_failure_message(pytest_temp), flush=True)
    return exit_code


if __name__ == "__main__":
    raise SystemExit(main())

