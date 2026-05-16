from __future__ import annotations

import importlib.util
import shutil
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = REPO_ROOT / "src"
DEV_BUILD_CACHE = REPO_ROOT / "build" / "dev" / "CMakeCache.txt"
STALE_TRACKED_REPORTS: tuple[Path, ...] = ()
REQUIRED_CORE_SYMBOLS = (
    "_fit_pure_neutral_native_ceres",
    "_fit_generic_native_least_squares",
    "_evaluate_generic_native_debug",
    "_solve_equilibrium_native",
    "_native_cppad_smoke",
    "_native_ipopt_smoke",
    "NativeSolutionError",
)


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


def _missing_core_symbols() -> tuple[str, ...]:
    try:
        import epcsaft._core as core
    except Exception:
        return ()
    return tuple(name for name in REQUIRED_CORE_SYMBOLS if not hasattr(core, name))


def _cppad_status() -> str:
    return _runtime_optional_dependency_status("cppad")


def _runtime_optional_dependency_status(name: str) -> str:
    try:
        import epcsaft

        info = epcsaft.runtime_build_info()
    except Exception:
        return "<unknown>"
    optional = info.get("optional_dependencies", {})
    if not isinstance(optional, dict):
        return "<unknown>"
    payload = optional.get(name, {})
    if not isinstance(payload, dict):
        return "<unknown>"
    return str(payload.get("status", "<unknown>"))


def _ceres_status() -> str:
    return _runtime_optional_dependency_status("ceres")


def _ipopt_status() -> str:
    return _runtime_optional_dependency_status("ipopt")


def _cmake_cache_value(name: str, cache_path: Path = DEV_BUILD_CACHE) -> str | None:
    if not cache_path.exists():
        return None
    prefix = f"{name}:"
    for line in cache_path.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith(prefix):
            return line.split("=", 1)[1].strip()
    return None


def _tool_path(name: str) -> str:
    return shutil.which(name) or "<missing>"


def _cmake_generator(cache_path: Path = DEV_BUILD_CACHE) -> str | None:
    if not cache_path.exists():
        return None
    for line in cache_path.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith("CMAKE_GENERATOR:INTERNAL="):
            return line.split("=", 1)[1].strip()
    return None


def _ninja_migration_recommendation(generator: str | None, ninja_path: str) -> str | None:
    if generator == "MinGW Makefiles" and ninja_path != "<missing>":
        return "uv run python scripts\\dev\\build_epcsaft.py --clean --generator ninja"
    return None


def _tracked_generated_count() -> int | None:
    output = _git_output("ls-files")
    if output is None:
        return None
    count = 0
    generated_extensions = {".png", ".svg", ".csv", ".md"}
    for raw_line in output.splitlines():
        path = Path(raw_line)
        parts = path.parts
        if "out" not in parts and "runs" not in parts:
            continue
        if not parts or parts[0] not in {"analyses", "scripts", "tests", "src"}:
            continue
        if path.suffix.lower() in generated_extensions:
            count += 1
    return count


def _stale_report_state() -> str:
    stale = [path.relative_to(REPO_ROOT).as_posix() for path in STALE_TRACKED_REPORTS if path.exists()]
    return ", ".join(stale) if stale else "<none>"


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
    print(f"git_branch: {branch}")
    print(f"git_head: {head}")
    print(f"uv: {_tool_path('uv')}")
    print(f"cmake: {_tool_path('cmake')}")
    ninja_path = _tool_path("ninja")
    print(f"ninja: {ninja_path}")
    cmake_generator = _cmake_generator()
    print(f"cmake_generator: {cmake_generator or '<missing>'}")
    generator_recommendation = _ninja_migration_recommendation(cmake_generator, ninja_path)
    print(f"build_generator_recommendation: {generator_recommendation or '<none>'}")
    print(f"epcsaft_import: {package_path if package_path else '<missing>'}")
    print(f"epcsaft_import_error: {package_error or '<none>'}")
    print(f"epcsaft_core: {core_path if core_path else '<missing>'}")
    print(f"epcsaft_core_error: {core_error or '<none>'}")
    missing_core_symbols = _missing_core_symbols() if core_path is not None else ()
    print(f"epcsaft_core_missing_symbols: {', '.join(missing_core_symbols) if missing_core_symbols else '<none>'}")
    print(f"ceres_configured: {_cmake_cache_value('EPCSAFT_ENABLE_CERES') or '<unconfigured>'}")
    print(f"cppad_status: {_cppad_status()}")
    print(f"ceres_status: {_ceres_status()}")
    print(f"ipopt_configured: {_cmake_cache_value('EPCSAFT_ENABLE_IPOPT') or '<unconfigured>'}")
    print(f"ipopt_status: {_ipopt_status()}")
    print(f"stale_generated_reports: {_stale_report_state()}")
    tracked_generated = _tracked_generated_count()
    print(f"tracked_generated_run_files: {tracked_generated if tracked_generated is not None else '<unknown>'}")

    if package_path is None:
        print("install_state: missing-package")
        print("next_command: uv sync --no-install-project")
        return 1
    if core_path is None:
        print("install_state: missing-core")
        print("next_command: uv run python scripts\\dev\\build_epcsaft.py")
        return 1
    if missing_core_symbols:
        print("install_state: stale-core")
        print("next_command: uv run python scripts\\dev\\build_epcsaft.py")
        return 1
    print("install_state: current")
    print("next_command: none")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
