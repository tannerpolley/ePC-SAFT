import argparse
import os
import shutil
import sys
import uuid
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from scripts.dev.native_runtime_env import apply_native_runtime_env

apply_native_runtime_env(os.environ)

from epcsaft.capability_evidence import TEST_SLICES, registry_targets

GENERIC_TEST_TARGETS = registry_targets("generic")
FAST_TEST_TARGETS = GENERIC_TEST_TARGETS
CONFIDENCE_TEST_TARGETS = registry_targets("confidence")
EQUILIBRIUM_CONFIDENCE_TEST_TARGETS = registry_targets("equilibrium-confidence")
EQUILIBRIUM_API_TEST_TARGETS = registry_targets("equilibrium-api")
ALL_TEST_TARGETS = registry_targets("all")
RUNTIME_TEST_TARGETS = registry_targets("runtime")
API_TEST_TARGETS = registry_targets("api")
NATIVE_TEST_TARGETS = registry_targets("native")
NATIVE_CONTRACT_TEST_TARGETS = registry_targets("native-contracts")
PROFILE_TEST_TARGETS = registry_targets("profile")
FULL_PROFILE_TEST_TARGETS = registry_targets("profile-full")
SLICE_TARGETS = {name: registry_targets(name) for name in TEST_SLICES}
LONG_NATIVE_TARGETS = {
    "tests/native/equilibrium",
    "tests/native/equilibrium/test_route_builders.py",
}
LONG_NATIVE_TARGETS_NOTE = (
    "Broad native equilibrium route-builder targets are intentionally guarded because they can take a long time. "
    "Use `uv run python run_pytest.py --native-contracts -q` for metadata/result-contract checks, "
    "or target a single test node. Pass `--allow-long-native-tests` only when you intentionally need the broad route suite."
)
FULL_PROFILE_MIN_TIMEOUT_SECONDS = 120
FULL_PROFILE_RUNTIME_NOTE = (
    "--profile-full runs runtime, MIAC, and regression profiling; "
    f"it can take about a minute locally, so allow at least {FULL_PROFILE_MIN_TIMEOUT_SECONDS} seconds."
)
SLICE_SELECTION_NOTE = (
    "Slice flags are mutually exclusive. Developers should normally start with "
    "`uv run python scripts/dev/validate_project.py quick` or `uv run python run_pytest.py -q`. "
    "Use `--all` only when you explicitly need the exhaustive historical suite. "
    "Extra positional pytest targets after a slice are appended and will run in addition to that slice."
)


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
    try:
        from scripts.dev.native_runtime_env import apply_native_runtime_env
    except ModuleNotFoundError:
        apply_native_runtime_env = None
    if apply_native_runtime_env is not None:
        apply_native_runtime_env(env)
    if profile:
        env["EPCSAFT_RUN_PERF"] = "1"
        env["ePCSAFT_RUN_PERF"] = "1"
    return env


def _pytest_args(
    pytest_args: list[str],
    pytest_temp: Path,
    generic: bool = False,
    confidence: bool = False,
    equilibrium_confidence: bool = False,
    equilibrium_api: bool = False,
    runtime: bool = False,
    api: bool = False,
    native: bool = False,
    native_contracts: bool = False,
    profile: bool = False,
    profile_full: bool = False,
    all_tests: bool = False,
    allow_long_native_tests: bool = False,
) -> list[str]:
    _reject_unbounded_native_targets(
        pytest_args,
        all_tests=all_tests,
        allow_long_native_tests=allow_long_native_tests,
    )
    cmd: list[str] = []
    has_predefined_targets = (
        generic
        or confidence
        or equilibrium_confidence
        or equilibrium_api
        or runtime
        or api
        or native
        or native_contracts
        or profile
        or profile_full
        or all_tests
    )
    if all_tests:
        cmd.extend(ALL_TEST_TARGETS)
    elif confidence:
        cmd.extend(CONFIDENCE_TEST_TARGETS)
    elif equilibrium_confidence:
        cmd.extend(EQUILIBRIUM_CONFIDENCE_TEST_TARGETS)
    elif equilibrium_api:
        cmd.extend(EQUILIBRIUM_API_TEST_TARGETS)
    elif generic:
        cmd.extend(GENERIC_TEST_TARGETS)
    elif runtime:
        cmd.extend(RUNTIME_TEST_TARGETS)
    elif api:
        cmd.extend(API_TEST_TARGETS)
    elif native:
        cmd.extend(NATIVE_TEST_TARGETS)
    elif native_contracts:
        cmd.extend(NATIVE_CONTRACT_TEST_TARGETS)
    elif profile:
        cmd.extend(PROFILE_TEST_TARGETS)
    elif profile_full:
        cmd.extend(FULL_PROFILE_TEST_TARGETS)

    if has_predefined_targets:
        cmd.extend(pytest_args)
    else:
        has_positional_target = any(not arg.startswith("-") for arg in pytest_args)
        if has_positional_target:
            cmd.extend(pytest_args)
        else:
            cmd.extend(FAST_TEST_TARGETS)
            cmd.extend(pytest_args)

    if not any(arg == "--basetemp" or arg.startswith("--basetemp=") for arg in cmd):
        cmd.extend(["--basetemp", str(pytest_temp)])
    return cmd


def _reject_unbounded_native_targets(
    pytest_args: list[str],
    *,
    all_tests: bool,
    allow_long_native_tests: bool,
) -> None:
    if all_tests or allow_long_native_tests or os.environ.get("EPCSAFT_ALLOW_LONG_NATIVE_TESTS") == "1":
        return
    for arg in pytest_args:
        if arg.startswith("-"):
            continue
        target_path = arg.split("::", 1)[0].replace("\\", "/").rstrip("/")
        if "::" in arg:
            continue
        if target_path in LONG_NATIVE_TARGETS:
            raise SystemExit(LONG_NATIVE_TARGETS_NOTE)


def _slice_listing_text() -> str:
    lines = [SLICE_SELECTION_NOTE, "", "Available slices:"]
    for name, targets in SLICE_TARGETS.items():
        lines.append(f"{name}:")
        for target in targets:
            lines.append(f"  {target}")
    return "\n".join(lines)


def _patch_windows_pytest_temp_acl() -> None:
    if os.name != "nt":
        return

    original_mkdir = Path.mkdir

    def sandbox_safe_mkdir(self, mode=0o777, parents=False, exist_ok=False):
        return original_mkdir(self, mode=0o777, parents=parents, exist_ok=exist_ok)

    Path.mkdir = sandbox_safe_mkdir


def _patch_pytest_cleanup() -> None:
    # Pytest's Windows cleanup hook can trip over restricted ACLs after
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
    parser = argparse.ArgumentParser(description=SLICE_SELECTION_NOTE)
    predefined = parser.add_mutually_exclusive_group()
    predefined.add_argument("--generic", action="store_true", help="Run the core generic test slice")
    predefined.add_argument(
        "--confidence",
        action="store_true",
        help="Run generic fast contracts plus native runtime contract tests",
    )
    predefined.add_argument(
        "--equilibrium-confidence",
        action="store_true",
        help="Run electrolyte equilibrium confidence contract tests; full reports remain env opt-in",
    )
    predefined.add_argument(
        "--equilibrium-api",
        action="store_true",
        help="Run fast equilibrium/speciation API tests for downstream-agent workflows",
    )
    predefined.add_argument("--runtime", action="store_true", help="Run runtime API and native contract tests")
    predefined.add_argument("--api", action="store_true", help="Run public API and regression API tests")
    predefined.add_argument("--native", action="store_true", help="Run native runtime contract tests")
    predefined.add_argument(
        "--native-contracts",
        action="store_true",
        help="Run fast native route metadata/result contract tests without broad solver route suites",
    )
    predefined.add_argument("--profile", action="store_true", help="Run the opt-in runtime-only profile test")
    predefined.add_argument(
        "--profile-full",
        action="store_true",
        help="Run all opt-in runtime, MIAC, and regression profile tests",
    )
    predefined.add_argument(
        "--all",
        dest="all_tests",
        action="store_true",
        help="Run the exhaustive historical test suite under tests/; this is intentionally opt-in",
    )
    parser.add_argument(
        "--list-slices", action="store_true", help="Print named test slices and exit without running pytest"
    )
    parser.add_argument(
        "--allow-long-native-tests",
        action="store_true",
        help="Allow broad known-slow native route-builder targets such as the full test_route_builders.py file",
    )
    args, pytest_args = parser.parse_known_args()

    if args.list_slices:
        print(_slice_listing_text())
        return 0

    repo_root = _repo_root()
    pytest_temp = _pytest_temp(repo_root)
    env = _pytest_env(pytest_temp, profile=args.profile or args.profile_full)
    src_root = repo_root / "src"
    sys.path.insert(0, str(src_root))
    env["PYTHONPATH"] = str(src_root)

    cmd = _pytest_args(
        pytest_args,
        pytest_temp,
        args.generic,
        confidence=args.confidence,
        equilibrium_confidence=args.equilibrium_confidence,
        equilibrium_api=args.equilibrium_api,
        runtime=args.runtime,
        api=args.api,
        native=args.native,
        native_contracts=args.native_contracts,
        profile=args.profile,
        profile_full=args.profile_full,
        all_tests=args.all_tests,
        allow_long_native_tests=args.allow_long_native_tests,
    )
    print("Running:", f"{sys.executable} -m pytest", " ".join(cmd), flush=True)
    if args.profile_full:
        print(FULL_PROFILE_RUNTIME_NOTE, flush=True)
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
