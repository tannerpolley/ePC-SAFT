from __future__ import annotations

import argparse
import os
import re
import shutil
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import NamedTuple

REPO_ROOT = Path(__file__).resolve().parents[2]
BUILD_DIR = REPO_ROOT / "build" / "dev"
PACKAGE_DIR = REPO_ROOT / "src" / "epcsaft"
LOG_FILE_NAME = "build_epcsaft.log"
STALE_LOCK_SECONDS = 120


class BuildProfile(NamedTuple):
    enable_ceres: bool
    enable_cppad: bool
    windows_parallel: str
    description: str


class BuildSettings(NamedTuple):
    enable_ceres: bool
    enable_cppad: bool
    parallel: str | None


BUILD_PROFILES: dict[str, BuildProfile] = {
    "fast": BuildProfile(
        enable_ceres=False,
        enable_cppad=True,
        windows_parallel="2",
        description="default local iteration profile: CppAD enabled, Ceres disabled",
    ),
    "cppad": BuildProfile(
        enable_ceres=False,
        enable_cppad=True,
        windows_parallel="2",
        description="explicit CppAD validation profile without the Ceres FetchContent build",
    ),
    "full": BuildProfile(
        enable_ceres=True,
        enable_cppad=True,
        windows_parallel="4",
        description="full optional native backend profile: Ceres and CppAD enabled",
    ),
    "minimal": BuildProfile(
        enable_ceres=False,
        enable_cppad=False,
        windows_parallel="10",
        description="smallest native extension profile: Ceres and CppAD disabled",
    ),
}


def _run(cmd: list[str], *, env: dict[str, str]) -> None:
    print("Running:", " ".join(cmd), flush=True)
    _append_log(f"\n[{_timestamp()}] Running: {' '.join(cmd)}\n")
    completed = subprocess.Popen(
        cmd,
        cwd=str(REPO_ROOT),
        env=env,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        bufsize=1,
    )
    assert completed.stdout is not None
    with _log_path().open("a", encoding="utf-8", errors="replace") as handle:
        for line in completed.stdout:
            print(line, end="", flush=True)
            handle.write(line)
    returncode = completed.wait()
    if returncode != 0:
        raise subprocess.CalledProcessError(returncode, cmd)


def _capture(cmd: list[str], *, env: dict[str, str]) -> str:
    return subprocess.check_output(cmd, cwd=str(REPO_ROOT), env=env, text=True).strip()


def _timestamp() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def _log_path() -> Path:
    BUILD_DIR.mkdir(parents=True, exist_ok=True)
    return BUILD_DIR / LOG_FILE_NAME


def _append_log(text: str) -> None:
    with _log_path().open("a", encoding="utf-8", errors="replace") as handle:
        handle.write(text)


def _pyproject_version() -> str:
    text = (REPO_ROOT / "pyproject.toml").read_text(encoding="utf-8")
    match = re.search(r'(?m)^version\s*=\s*"([^"]+)"', text)
    if not match:
        raise RuntimeError("Could not derive package version from pyproject.toml")
    return match.group(1)


def _env() -> dict[str, str]:
    temp_root = REPO_ROOT / "build" / "temp"
    temp_root.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    env["TMP"] = str(temp_root.resolve())
    env["TEMP"] = str(temp_root.resolve())
    env["TMPDIR"] = str(temp_root.resolve())
    env.setdefault("PIP_DISABLE_PIP_VERSION_CHECK", "1")
    return env


def _clean() -> None:
    # Remove the importable extension first. If Windows has it locked, fail
    # before deleting the reusable CMake build tree.
    for artifact in PACKAGE_DIR.glob("_core*.pyd"):
        _remove_extension_artifact(artifact)
    for artifact in PACKAGE_DIR.glob("_core*.so"):
        _remove_extension_artifact(artifact)
    shutil.rmtree(BUILD_DIR, ignore_errors=True)


def _remove_extension_artifact(artifact: Path) -> None:
    try:
        artifact.unlink()
    except PermissionError as exc:
        message = (
            f"Unable to remove {artifact}. The compiled extension is probably "
            "loaded by an active Python process. Close Python terminals, IDE "
            "test runners, or parallel workers that imported epcsaft._core, "
            "then rerun the build."
        )
        raise PermissionError(message) from exc


def _configured_generator() -> str | None:
    return _cmake_cache_value("CMAKE_GENERATOR")


def _cmake_cache_value(name: str) -> str | None:
    cache = BUILD_DIR / "CMakeCache.txt"
    if not cache.exists():
        return None
    prefix = f"{name}:"
    for line in cache.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith(prefix):
            return line.split("=", 1)[1].strip()
    return None


def _native_artifacts() -> list[Path]:
    return sorted([*PACKAGE_DIR.glob("_core*.pyd"), *PACKAGE_DIR.glob("_core*.so")])


def _last_ninja_target() -> str | None:
    log = BUILD_DIR / ".ninja_log"
    if not log.exists():
        return None
    for line in reversed(log.read_text(encoding="utf-8", errors="replace").splitlines()):
        if not line or line.startswith("#"):
            continue
        parts = line.split("\t")
        if len(parts) >= 4:
            return parts[3]
    return None


def _repo_build_processes() -> list[str]:
    if os.name != "nt":
        return []
    root = str(REPO_ROOT).replace("'", "''")
    current_pid = os.getpid()
    command = (
        "$root = '" + root + "'; "
        f"$currentPid = {current_pid}; "
        "Get-CimInstance Win32_Process | "
        "Where-Object { $_.Name -in @('python.exe','uv.exe','cmake.exe','ninja.exe') "
        "-and $_.ProcessId -ne $currentPid "
        '-and $_.CommandLine -like "*$root*" } | '
        'ForEach-Object { "$($_.ProcessId) $($_.Name) $($_.CommandLine)" }'
    )
    try:
        output = subprocess.check_output(
            ["powershell.exe", "-NoProfile", "-Command", command],
            cwd=str(REPO_ROOT),
            text=True,
            stderr=subprocess.DEVNULL,
            timeout=5,
        )
    except (subprocess.SubprocessError, OSError):
        return []
    processes = []
    for line in output.splitlines():
        process = line.strip()
        normalized = process.lower().replace("\\", "/")
        is_build_script = "scripts/dev/build_epcsaft.py" in normalized and "--status" not in normalized
        is_cmake_build = "cmake" in normalized and "--build" in normalized
        is_ninja_build = "ninja.exe" in normalized
        if is_build_script or is_cmake_build or is_ninja_build:
            processes.append(process)
    return processes


def _status_lines(*, stale_lock_seconds: int = STALE_LOCK_SECONDS) -> list[str]:
    generator = _configured_generator() or "<unconfigured>"
    ceres = _cmake_cache_value("EPCSAFT_ENABLE_CERES") or "<unconfigured>"
    system_ceres = _cmake_cache_value("EPCSAFT_USE_SYSTEM_CERES") or "<unconfigured>"
    ceres_dir = _cmake_cache_value("Ceres_DIR") or "<unconfigured>"
    cppad = _cmake_cache_value("EPCSAFT_ENABLE_CPPAD") or "<unconfigured>"
    artifacts = _native_artifacts()
    lock = BUILD_DIR / ".ninja_lock"
    lock_present = lock.exists()
    lock_age = time.time() - lock.stat().st_mtime if lock_present else None
    stale_lock = bool(lock_present and lock_age is not None and lock_age >= stale_lock_seconds)
    processes = _repo_build_processes()

    lines = [
        f"repo_root: {REPO_ROOT}",
        f"build_dir: {BUILD_DIR}",
        f"configured_generator: {generator}",
        f"ceres_configured: {ceres}",
        f"system_ceres_configured: {system_ceres}",
        f"ceres_dir: {ceres_dir}",
        f"cppad_configured: {cppad}",
        f"profile_hint: {_profile_hint(ceres=ceres, cppad=cppad)}",
        f"native_core: {'present' if artifacts else 'missing'}",
        "native_core_paths: " + (", ".join(str(path) for path in artifacts) if artifacts else "<none>"),
        f"ninja_lock: {'present' if lock_present else 'missing'}",
        f"stale_ninja_lock: {'true' if stale_lock else 'false'}",
        f"last_ninja_target: {_last_ninja_target() or '<none>'}",
        f"live_build_processes: {len(processes)}",
        f"log_file: {_log_path()}",
    ]
    lines.extend(f"live_build_process: {process}" for process in processes)
    if stale_lock:
        lines.append(
            "stale_lock_remediation: inspect live_build_process entries; stop only repo-owned build processes before retrying"
        )
    return lines


def _profile_hint(*, ceres: str, cppad: str) -> str:
    if "<unconfigured>" in {ceres, cppad}:
        return "<unconfigured>"
    normalized = (ceres.upper(), cppad.upper())
    if normalized == ("OFF", "ON"):
        return "fast/cppad"
    if normalized == ("ON", "ON"):
        return "full"
    if normalized == ("OFF", "OFF"):
        return "minimal"
    return "custom"


def _print_status() -> None:
    for line in _status_lines():
        print(line)


def _warn_if_stale_build_lock() -> None:
    lines = _status_lines()
    if "stale_ninja_lock: true" not in lines:
        return
    print("Warning: stale Ninja lock detected for build/dev.", flush=True)
    for line in lines:
        if line.startswith(
            ("ninja_lock:", "stale_ninja_lock:", "last_ninja_target:", "live_build_process", "stale_lock")
        ):
            print(line, flush=True)


def _generator_args(env: dict[str, str], configured_generator: str | None = None) -> list[str]:
    requested = env.get("EPCSAFT_CMAKE_GENERATOR", "").strip().lower()
    if requested == "auto":
        requested = ""

    known = {
        "ninja": "Ninja",
        "mingw": "MinGW Makefiles",
        "mingw makefiles": "MinGW Makefiles",
    }
    if configured_generator:
        if not requested:
            return []
        target = known.get(requested)
        if target and target == configured_generator:
            return []
        raise RuntimeError(
            "build/dev is already configured with "
            f"{configured_generator!r}. Use --clean before switching to {target or requested!r}."
        )

    if requested in {"ninja", ""} and shutil.which("ninja", path=env.get("PATH")):
        return ["-G", "Ninja"]
    if (
        requested in {"mingw", "mingw makefiles"}
        and os.name == "nt"
        and shutil.which("mingw32-make", path=env.get("PATH"))
    ):
        return ["-G", "MinGW Makefiles"]
    if os.name == "nt" and shutil.which("mingw32-make", path=env.get("PATH")):
        return ["-G", "MinGW Makefiles"]
    return []


def _configure(
    env: dict[str, str],
    *,
    enable_ceres: bool,
    use_system_ceres: bool,
    ceres_dir: Path | None,
    enable_cppad: bool,
    use_system_cppad: bool,
) -> None:
    pybind11_dir = _capture([sys.executable, "-m", "pybind11", "--cmakedir"], env=env)
    cmd = [
        "cmake",
        "-S",
        str(REPO_ROOT),
        "-B",
        str(BUILD_DIR),
        "-DEPCSAFT_DEV_INPLACE=ON",
        f"-DEPCSAFT_ENABLE_CERES={'ON' if enable_ceres else 'OFF'}",
        f"-DEPCSAFT_USE_SYSTEM_CERES={'ON' if use_system_ceres else 'OFF'}",
        f"-DEPCSAFT_ENABLE_CPPAD={'ON' if enable_cppad else 'OFF'}",
        f"-DEPCSAFT_USE_SYSTEM_CPPAD={'ON' if use_system_cppad else 'OFF'}",
        "-DCMAKE_TRY_COMPILE_TARGET_TYPE=STATIC_LIBRARY",
        f"-DSKBUILD_PROJECT_VERSION={_pyproject_version()}",
        f"-DPython_EXECUTABLE={sys.executable}",
        f"-Dpybind11_DIR={pybind11_dir}",
    ]
    if ceres_dir is not None:
        cmd.append(f"-DCeres_DIR={ceres_dir}")
    cmd.extend(_generator_args(env, _configured_generator()))
    _run(cmd, env=env)


def _build(env: dict[str, str], parallel: str | None) -> None:
    cmd = ["cmake", "--build", str(BUILD_DIR)]
    if parallel:
        cmd.extend(["--parallel", parallel])
    _warn_if_stale_build_lock()
    _run(cmd, env=env)


def _verify_native_import(env: dict[str, str]) -> None:
    cmd = [
        sys.executable,
        "-c",
        "import epcsaft._core as core; print(core.__file__)",
    ]
    completed = subprocess.run(
        cmd,
        cwd=str(REPO_ROOT),
        env=env,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    _append_log(f"[{_timestamp()}] Running native import check\n{completed.stdout}")
    if completed.returncode != 0:
        print(completed.stdout, end="", flush=True)
        raise subprocess.CalledProcessError(completed.returncode, cmd)
    core_path = completed.stdout.strip()
    print(f"native import OK: {core_path}", flush=True)


def _timed(label: str, action) -> float:
    start = time.perf_counter()
    action()
    elapsed = time.perf_counter() - start
    print(f"Timing: {label} completed in {elapsed:.2f}s", flush=True)
    return elapsed


def _resolve_settings(args: argparse.Namespace) -> BuildSettings:
    profile = BUILD_PROFILES[args.profile]
    enable_ceres = profile.enable_ceres if args.enable_ceres is None else bool(args.enable_ceres)
    enable_cppad = profile.enable_cppad if args.enable_cppad is None else bool(args.enable_cppad)
    if args.use_system_cppad:
        enable_cppad = True
    if args.use_system_ceres or args.ceres_dir is not None:
        enable_ceres = True

    parallel = args.parallel
    if parallel is None:
        if os.name == "nt":
            parallel = BUILD_PROFILES["full"].windows_parallel if enable_ceres else profile.windows_parallel
        elif enable_ceres:
            parallel = BUILD_PROFILES["full"].windows_parallel
    return BuildSettings(enable_ceres=enable_ceres, enable_cppad=enable_cppad, parallel=parallel)


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Build epcsaft._core in-place with direct CMake.")
    parser.add_argument("--clean", action="store_true", help="Remove build/dev and any in-place _core artifact first.")
    parser.add_argument(
        "--configure-only", action="store_true", help="Configure the CMake dev build tree without building."
    )
    parser.add_argument(
        "--build-only", action="store_true", help="Build the existing CMake dev build tree without reconfiguring."
    )
    parser.add_argument(
        "--status", action="store_true", help="Report build/dev status without configuring or building."
    )
    parser.add_argument("--parallel", help="Optional CMake build parallelism value.")
    parser.add_argument(
        "--profile",
        choices=tuple(BUILD_PROFILES),
        default="fast",
        help=(
            "Native dependency profile. fast/cppad enable CppAD without Ceres; full enables Ceres and CppAD; "
            "minimal disables both. Explicit --enable-* or --disable-* flags override the profile."
        ),
    )
    parser.add_argument(
        "--generator",
        choices=("auto", "ninja", "mingw"),
        default="auto",
        help="CMake generator for a new configure. Auto prefers Ninja when available.",
    )
    parser.set_defaults(enable_ceres=None, enable_cppad=None)
    parser.add_argument(
        "--enable-ceres",
        dest="enable_ceres",
        action="store_true",
        help="Enable Ceres Solver support, overriding the selected profile.",
    )
    parser.add_argument(
        "--disable-ceres",
        dest="enable_ceres",
        action="store_false",
        help="Disable Ceres Solver support for this build.",
    )
    parser.add_argument(
        "--use-system-ceres",
        action="store_true",
        help="Use an installed Ceres Solver package. Implies --enable-ceres.",
    )
    parser.add_argument(
        "--ceres-dir",
        type=Path,
        help="Directory containing CeresConfig.cmake for a prebuilt/exported Ceres package. Implies --use-system-ceres.",
    )
    parser.add_argument(
        "--enable-cppad",
        dest="enable_cppad",
        action="store_true",
        help="Enable package-wide CppAD support, overriding the selected profile.",
    )
    parser.add_argument(
        "--disable-cppad",
        dest="enable_cppad",
        action="store_false",
        help="Disable package-wide CppAD support for this build.",
    )
    parser.add_argument(
        "--use-system-cppad",
        action="store_true",
        help="Use an installed CppAD include tree. Implies --enable-cppad.",
    )
    return parser


def main() -> int:
    parser = _parser()
    args = parser.parse_args()
    if args.status:
        _print_status()
        return 0
    if args.clean and args.build_only:
        parser.error("--clean cannot be combined with --build-only")
    if args.configure_only and args.build_only:
        parser.error("--configure-only cannot be combined with --build-only")
    settings = _resolve_settings(args)
    use_system_ceres = bool(args.use_system_ceres or args.ceres_dir is not None)

    total_start = time.perf_counter()
    if args.clean:
        _timed("clean", _clean)

    env = _env()
    if args.generator != "auto":
        env["EPCSAFT_CMAKE_GENERATOR"] = args.generator
    print(
        "Build profile: "
        f"{args.profile} ({BUILD_PROFILES[args.profile].description}); "
        f"Ceres={'ON' if settings.enable_ceres else 'OFF'}, "
        f"CeresSource={('system' if use_system_ceres else 'FetchContent') if settings.enable_ceres else 'disabled'}, "
        f"CppAD={'ON' if settings.enable_cppad else 'OFF'}, "
        f"parallel={settings.parallel or '<generator default>'}",
        flush=True,
    )
    if not args.build_only:
        _timed(
            "configure",
            lambda: _configure(
                env,
                enable_ceres=settings.enable_ceres,
                use_system_ceres=use_system_ceres,
                ceres_dir=args.ceres_dir,
                enable_cppad=settings.enable_cppad,
                use_system_cppad=args.use_system_cppad,
            ),
        )
    else:
        print("Timing: configure skipped (--build-only)", flush=True)
        print("Build profile does not reconfigure an existing CMake tree when --build-only is used.", flush=True)
    if not args.configure_only:
        _timed("build", lambda: _build(env, settings.parallel))
        _timed("native import", lambda: _verify_native_import(env))
    else:
        print("Timing: build skipped (--configure-only)", flush=True)
    print(f"Timing: total completed in {time.perf_counter() - total_start:.2f}s", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
