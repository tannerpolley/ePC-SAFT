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

try:
    from scripts.dev.native_runtime_env import apply_native_runtime_env
except ModuleNotFoundError:  # pragma: no cover - direct script execution
    from native_runtime_env import apply_native_runtime_env

REPO_ROOT = Path(__file__).resolve().parents[2]
BUILD_DIR = REPO_ROOT / "build" / "dev"
PACKAGE_DIR = REPO_ROOT / "src" / "epcsaft"
LOG_FILE_NAME = "build_epcsaft.log"
STALE_LOCK_SECONDS = 120


class BuildProfile(NamedTuple):
    enable_ipopt: bool
    windows_parallel: str
    description: str


class BuildSettings(NamedTuple):
    enable_ipopt: bool
    parallel: str | None


BUILD_PROFILES: dict[str, BuildProfile] = {
    "fast": BuildProfile(
        enable_ipopt=False,
        windows_parallel="4",
        description="default native dependency profile: Ceres and CppAD enabled, Ipopt disabled",
    ),
    "full": BuildProfile(
        enable_ipopt=False,
        windows_parallel="4",
        description="full native dependency profile: Ceres and CppAD enabled, Ipopt disabled",
    ),
    "ipopt": BuildProfile(
        enable_ipopt=True,
        windows_parallel="4",
        description="native solver dependency profile: Ceres, CppAD, and system Ipopt enabled",
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


def _resolve_optional_dir(raw_path: Path | str | None, *, label: str) -> Path | None:
    if raw_path is None:
        return None
    path = Path(raw_path).expanduser().resolve()
    if not path.is_dir():
        raise FileNotFoundError(f"{label} does not exist or is not a directory: {path}")
    return path


def _ipopt_root_prefers_msvc(ipopt_root: Path | None) -> bool:
    if os.name != "nt" or ipopt_root is None:
        return False
    lib_dir = ipopt_root / "lib"
    if not lib_dir.is_dir():
        return False
    has_msvc_import_lib = any((lib_dir / name).is_file() for name in ("ipopt.lib", "ipopt-3.lib"))
    has_mingw_import_lib = any((lib_dir / name).is_file() for name in ("libipopt.dll.a", "ipopt.dll.a"))
    return has_msvc_import_lib and not has_mingw_import_lib


def _configured_cxx_compiler() -> str | None:
    return _cmake_cache_value("CMAKE_CXX_COMPILER")


def _configured_compiler_is_msvc() -> bool:
    compiler = _configured_cxx_compiler()
    if not compiler:
        return False
    return Path(compiler).name.lower() == "cl.exe"


def _find_vsdevcmd() -> Path:
    explicit = os.environ.get("EPCSAFT_VSDEVCMD")
    if explicit:
        path = Path(explicit).expanduser().resolve()
        if path.is_file():
            return path
    vswhere = Path(os.environ.get("ProgramFiles(x86)", r"C:\Program Files (x86)")) / (
        "Microsoft Visual Studio/Installer/vswhere.exe"
    )
    if vswhere.is_file():
        output = subprocess.check_output(
            [
                str(vswhere),
                "-latest",
                "-products",
                "*",
                "-requires",
                "Microsoft.VisualStudio.Component.VC.Tools.x86.x64",
                "-property",
                "installationPath",
            ],
            text=True,
        ).strip()
        if output:
            candidate = Path(output) / "Common7" / "Tools" / "VsDevCmd.bat"
            if candidate.is_file():
                return candidate
    default_path = Path(r"C:\Program Files (x86)\Microsoft Visual Studio\2022\BuildTools\Common7\Tools\VsDevCmd.bat")
    if default_path.is_file():
        return default_path
    raise FileNotFoundError("MSVC toolchain requested, but VsDevCmd.bat was not found.")


def _load_msvc_env(env: dict[str, str]) -> dict[str, str]:
    vsdevcmd = _find_vsdevcmd()
    command = f'call "{vsdevcmd}" -arch=x64 -host_arch=x64 >nul && set'
    output = subprocess.check_output(command, env=env, text=True, shell=True)
    updated = env.copy()
    for line in output.splitlines():
        if "=" not in line:
            continue
        key, value = line.split("=", 1)
        updated[key] = value
    return updated


def _apply_toolchain_env(
    env: dict[str, str],
    *,
    toolchain: str,
    enable_ipopt: bool,
    ipopt_root: Path | None,
) -> dict[str, str]:
    load_msvc = toolchain == "msvc" or (
        toolchain == "auto"
        and os.name == "nt"
        and (
            (enable_ipopt and _ipopt_root_prefers_msvc(ipopt_root))
            or _configured_compiler_is_msvc()
        )
    )
    if not load_msvc:
        return env
    configured = _configured_cxx_compiler()
    if configured and Path(configured).name.lower() != "cl.exe":
        raise RuntimeError(
            "build/dev is already configured with a non-MSVC compiler. "
            "Use --clean before switching to the MSVC Ipopt toolchain."
        )
    return _load_msvc_env(env)


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
    compiler = _configured_cxx_compiler() or "<unconfigured>"
    build_type = _cmake_cache_value("CMAKE_BUILD_TYPE") or "<unconfigured>"
    ceres = _cmake_cache_value("EPCSAFT_ENABLE_CERES") or "<unconfigured>"
    system_ceres = _cmake_cache_value("EPCSAFT_USE_SYSTEM_CERES") or "<unconfigured>"
    ceres_dir = _cmake_cache_value("Ceres_DIR") or "<unconfigured>"
    cppad = _cmake_cache_value("EPCSAFT_ENABLE_CPPAD") or "<unconfigured>"
    ipopt = _cmake_cache_value("EPCSAFT_ENABLE_IPOPT") or "<unconfigured>"
    system_ipopt = _cmake_cache_value("EPCSAFT_USE_SYSTEM_IPOPT") or "<unconfigured>"
    ipopt_dir = _cmake_cache_value("Ipopt_DIR") or "<unconfigured>"
    ipopt_root = _cmake_cache_value("EPCSAFT_IPOPT_ROOT") or "<unconfigured>"
    runtime_env = apply_native_runtime_env(
        {},
        cache_path=BUILD_DIR / "CMakeCache.txt",
        ipopt_enabled=ipopt.upper() == "ON",
    )
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
        f"configured_cxx_compiler: {compiler}",
        f"build_type: {build_type}",
        f"ceres_configured: {ceres}",
        f"system_ceres_configured: {system_ceres}",
        f"ceres_dir: {ceres_dir}",
        f"cppad_configured: {cppad}",
        f"ipopt_configured: {ipopt}",
        f"system_ipopt_configured: {system_ipopt}",
        f"ipopt_dir: {ipopt_dir}",
        f"ipopt_root: {ipopt_root}",
        f"ipopt_runtime_dll_dir: {runtime_env.ipopt_runtime_dir if runtime_env.ipopt_runtime_dir else '<none>'}",
        f"ipopt_runtime_env_applied: {'true' if runtime_env.applied else 'false'}",
        f"profile_hint: {_profile_hint(ceres=ceres, cppad=cppad, ipopt=ipopt)}",
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


def _profile_hint(*, ceres: str, cppad: str, ipopt: str) -> str:
    if "<unconfigured>" in {ceres, cppad, ipopt}:
        return "<unconfigured>"
    normalized = (ceres.upper(), cppad.upper(), ipopt.upper())
    if normalized == ("ON", "ON", "OFF"):
        return "fast/full"
    if normalized == ("ON", "ON", "ON"):
        return "ipopt"
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
    use_system_ceres: bool,
    ceres_dir: Path | None,
    use_system_cppad: bool,
    enable_ipopt: bool,
    use_system_ipopt: bool,
    ipopt_dir: Path | None,
    ipopt_root: Path | None,
    build_type: str,
) -> None:
    pybind11_dir = _capture([sys.executable, "-m", "pybind11", "--cmakedir"], env=env)
    cmd = [
        "cmake",
        "-S",
        str(REPO_ROOT),
        "-B",
        str(BUILD_DIR),
        "-DEPCSAFT_DEV_INPLACE=ON",
        "-DEPCSAFT_ENABLE_CERES=ON",
        f"-DEPCSAFT_USE_SYSTEM_CERES={'ON' if use_system_ceres else 'OFF'}",
        "-DEPCSAFT_ENABLE_CPPAD=ON",
        f"-DEPCSAFT_USE_SYSTEM_CPPAD={'ON' if use_system_cppad else 'OFF'}",
        f"-DEPCSAFT_ENABLE_IPOPT={'ON' if enable_ipopt else 'OFF'}",
        f"-DEPCSAFT_USE_SYSTEM_IPOPT={'ON' if use_system_ipopt else 'OFF'}",
        "-DCMAKE_TRY_COMPILE_TARGET_TYPE=STATIC_LIBRARY",
        f"-DCMAKE_BUILD_TYPE={build_type}",
        f"-DSKBUILD_PROJECT_VERSION={_pyproject_version()}",
        f"-DPython_EXECUTABLE={sys.executable}",
        f"-Dpybind11_DIR={pybind11_dir}",
    ]
    if ceres_dir is not None:
        cmd.append(f"-DCeres_DIR={ceres_dir}")
    if ipopt_dir is not None:
        cmd.append(f"-DIpopt_DIR={ipopt_dir}")
    else:
        cmd.extend(["-U", "Ipopt_DIR"])
    if ipopt_root is not None:
        cmd.append(f"-DEPCSAFT_IPOPT_ROOT={ipopt_root}")
    else:
        cmd.append("-DEPCSAFT_IPOPT_ROOT=")
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
    enable_ipopt = profile.enable_ipopt if args.enable_ipopt is None else bool(args.enable_ipopt)
    if (
        args.use_system_ipopt
        or args.ipopt_dir is not None
        or args.ipopt_root is not None
        or os.environ.get("EPCSAFT_IPOPT_ROOT")
        or os.environ.get("EPCSAFT_PEP517_IPOPT_ROOT")
    ):
        enable_ipopt = True

    parallel = args.parallel
    if parallel is None:
        if enable_ipopt:
            parallel = BUILD_PROFILES["ipopt"].windows_parallel
        elif os.name == "nt":
            parallel = profile.windows_parallel
        else:
            parallel = BUILD_PROFILES["full"].windows_parallel
    return BuildSettings(
        enable_ipopt=enable_ipopt,
        parallel=parallel,
    )


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
        "--build-type",
        default=os.environ.get("EPCSAFT_CMAKE_BUILD_TYPE", "Release"),
        help="CMake build type for single-config generators. Defaults to Release.",
    )
    parser.add_argument(
        "--profile",
        choices=tuple(BUILD_PROFILES),
        default="fast",
        help=(
            "Native dependency profile. fast/full require Ceres and CppAD; "
            "ipopt additionally enables the system Ipopt dependency."
        ),
    )
    parser.add_argument(
        "--generator",
        choices=("auto", "ninja", "mingw"),
        default="auto",
        help="CMake generator for a new configure. Auto prefers Ninja when available.",
    )
    parser.add_argument(
        "--toolchain",
        choices=("auto", "current", "msvc"),
        default="auto",
        help=(
            "Compiler environment for a new configure or build-only run. Auto loads MSVC on Windows "
            "when --ipopt-root points at an MSVC-style Ipopt tree."
        ),
    )
    parser.set_defaults(enable_ipopt=None)
    parser.add_argument(
        "--use-system-ceres",
        action="store_true",
        help="Use an installed Ceres Solver package instead of FetchContent.",
    )
    parser.add_argument(
        "--ceres-dir",
        type=Path,
        help="Directory containing CeresConfig.cmake for a prebuilt/exported Ceres package. Implies --use-system-ceres.",
    )
    parser.add_argument(
        "--use-system-cppad",
        action="store_true",
        help="Use an installed CppAD include tree instead of FetchContent.",
    )
    parser.add_argument(
        "--enable-ipopt",
        dest="enable_ipopt",
        action="store_true",
        help="Enable native Ipopt support. Requires a system Ipopt package.",
    )
    parser.add_argument(
        "--disable-ipopt",
        dest="enable_ipopt",
        action="store_false",
        help="Disable native Ipopt support for this build.",
    )
    parser.add_argument(
        "--use-system-ipopt",
        action="store_true",
        help="Use an installed Ipopt package. Implies --enable-ipopt.",
    )
    parser.add_argument(
        "--ipopt-dir",
        type=Path,
        help="Directory containing IpoptConfig.cmake for a native Ipopt package. Implies --use-system-ipopt.",
    )
    parser.add_argument(
        "--ipopt-root",
        type=Path,
        help=(
            "Root directory for an installed Ipopt tree with include/, lib/, and optional bin/. "
            "Implies --use-system-ipopt."
        ),
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
    if args.ipopt_dir is not None and args.ipopt_root is not None:
        parser.error("--ipopt-dir and --ipopt-root are mutually exclusive")
    settings = _resolve_settings(args)
    use_system_ceres = bool(args.use_system_ceres or args.ceres_dir is not None)
    ipopt_root_env = os.environ.get("EPCSAFT_IPOPT_ROOT") or os.environ.get("EPCSAFT_PEP517_IPOPT_ROOT")
    ipopt_root = _resolve_optional_dir(
        args.ipopt_root if args.ipopt_root is not None else ipopt_root_env,
        label="Ipopt root",
    )
    use_system_ipopt = bool(
        settings.enable_ipopt or args.use_system_ipopt or args.ipopt_dir is not None or ipopt_root is not None
    )

    total_start = time.perf_counter()
    if args.clean:
        _timed("clean", _clean)

    env = _env()
    if args.generator != "auto":
        env["EPCSAFT_CMAKE_GENERATOR"] = args.generator
    env = _apply_toolchain_env(
        env,
        toolchain=args.toolchain,
        enable_ipopt=settings.enable_ipopt,
        ipopt_root=ipopt_root,
    )
    runtime_env = apply_native_runtime_env(env, ipopt_root=ipopt_root, ipopt_enabled=settings.enable_ipopt)
    print(
        "Build profile: "
        f"{args.profile} ({BUILD_PROFILES[args.profile].description}); "
        "Ceres=ON, "
        f"CeresSource={'system' if use_system_ceres else 'FetchContent'}, "
        "CppAD=ON, "
        f"Ipopt={'ON' if settings.enable_ipopt else 'OFF'}, "
        f"IpoptSource={('system' if use_system_ipopt else 'disabled') if settings.enable_ipopt else 'disabled'}, "
        f"IpoptRoot={ipopt_root if ipopt_root is not None else '<none>'}, "
        f"IpoptRuntimeDllDir={runtime_env.ipopt_runtime_dir if runtime_env.ipopt_runtime_dir else '<none>'}, "
        f"Toolchain={args.toolchain}, "
        f"BuildType={args.build_type}, "
        f"parallel={settings.parallel or '<generator default>'}",
        flush=True,
    )
    if not args.build_only:
        _timed(
            "configure",
            lambda: _configure(
                env,
                use_system_ceres=use_system_ceres,
                ceres_dir=args.ceres_dir,
                use_system_cppad=args.use_system_cppad,
                enable_ipopt=settings.enable_ipopt,
                use_system_ipopt=use_system_ipopt,
                ipopt_dir=args.ipopt_dir,
                ipopt_root=ipopt_root,
                build_type=args.build_type,
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
