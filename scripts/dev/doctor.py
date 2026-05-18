from __future__ import annotations

import argparse
import importlib
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
    "_fit_generic_native_ceres",
    "_evaluate_generic_native_debug",
    "_native_cppad_smoke",
    "_native_ipopt_smoke",
    "NativeSolutionError",
)

try:
    from scripts.dev.native_runtime_env import apply_to_current_process
except ModuleNotFoundError:  # pragma: no cover - direct script execution
    from native_runtime_env import apply_to_current_process


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
        module = importlib.import_module(module_name)
    except Exception as exc:
        return None, f"{type(exc).__name__}: {exc}"
    origin = getattr(module, "__file__", None)
    if origin is None:
        return None, None
    return Path(origin).resolve(), None


def _missing_core_symbols() -> tuple[str, ...]:
    try:
        import epcsaft._core as core
    except Exception as exc:
        raise RuntimeError(f"epcsaft._core import failed: {type(exc).__name__}: {exc}") from exc
    return tuple(name for name in REQUIRED_CORE_SYMBOLS if not hasattr(core, name))


def _missing_core_symbols_or_error() -> tuple[tuple[str, ...], str | None]:
    try:
        return _missing_core_symbols(), None
    except RuntimeError as exc:
        return (), str(exc)


def _cppad_status() -> str:
    return _runtime_native_dependency_status("cppad")


def _runtime_native_dependency_status(name: str) -> str:
    try:
        import epcsaft

        info = epcsaft.runtime_build_info()
    except Exception:
        return "<unknown>"
    native_dependencies = info.get("native_dependencies", {})
    if not isinstance(native_dependencies, dict):
        return "<unknown>"
    payload = native_dependencies.get(name, {})
    if not isinstance(payload, dict):
        return "<unknown>"
    return str(payload.get("status", "<unknown>"))


def _ceres_status() -> str:
    return _runtime_native_dependency_status("ceres")


def _ipopt_status() -> str:
    return _runtime_native_dependency_status("ipopt")


def _ipopt_available() -> bool:
    try:
        import epcsaft

        info = epcsaft.runtime_build_info()
    except Exception:
        return False
    native_dependencies = info.get("native_dependencies", {})
    if not isinstance(native_dependencies, dict):
        return False
    payload = native_dependencies.get("ipopt", {})
    if not isinstance(payload, dict):
        return False
    return payload.get("status") == "enabled_available" and bool(payload.get("available"))


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Check local ePC-SAFT development runtime health.")
    parser.add_argument(
        "--require-ipopt",
        action="store_true",
        help="Fail unless native Ipopt is compiled and available in the current runtime.",
    )
    return parser


def _cmake_cache_value(name: str, cache_path: Path = DEV_BUILD_CACHE) -> str | None:
    if not cache_path.exists():
        return None
    prefix = f"{name}:"
    for line in cache_path.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith(prefix):
            return line.split("=", 1)[1].strip()
    return None


def _tool_path(name: str) -> str:
    found = shutil.which(name)
    if found is not None:
        return found
    if name == "uv" and sys.platform.startswith("win"):
        local_uv = Path.home() / ".local" / "bin" / "uv.exe"
        if local_uv.is_file():
            return str(local_uv)
    return "<missing>"


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


def main(argv: list[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    if str(SRC_ROOT) not in sys.path:
        sys.path.insert(0, str(SRC_ROOT))
    runtime_env = apply_to_current_process()

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
    missing_core_symbols, missing_core_error = _missing_core_symbols_or_error() if core_path is not None else ((), None)
    print(f"epcsaft_core_missing_symbols: {', '.join(missing_core_symbols) if missing_core_symbols else '<none>'}")
    print(f"epcsaft_core_symbol_error: {missing_core_error or '<none>'}")
    print(f"ceres_configured: {_cmake_cache_value('EPCSAFT_ENABLE_CERES') or '<unconfigured>'}")
    print(f"cppad_status: {_cppad_status()}")
    print(f"ceres_status: {_ceres_status()}")
    print(f"ipopt_configured: {_cmake_cache_value('EPCSAFT_ENABLE_IPOPT') or '<unconfigured>'}")
    print(f"ipopt_root: {_cmake_cache_value('EPCSAFT_IPOPT_ROOT') or '<unconfigured>'}")
    print(f"ipopt_runtime_dll_dir: {runtime_env.ipopt_runtime_dir if runtime_env.ipopt_runtime_dir else '<none>'}")
    print(f"ipopt_runtime_env_applied: {'true' if runtime_env.applied else 'false'}")
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
    if missing_core_error is not None:
        print("install_state: broken-core-import")
        print("next_command: uv run python scripts\\dev\\build_epcsaft.py")
        return 1
    if missing_core_symbols:
        print("install_state: stale-core")
        print("next_command: uv run python scripts\\dev\\build_epcsaft.py")
        return 1
    if args.require_ipopt and not _ipopt_available():
        print("install_state: missing-ipopt")
        print(
            "next_command: uv run python scripts\\dev\\build_epcsaft.py --profile ipopt --ipopt-root <IpoptRoot>"
        )
        return 1
    print("install_state: current")
    print("next_command: none")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
