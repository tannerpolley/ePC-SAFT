from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
BUILD_DIR = REPO_ROOT / "build" / "dev"
PACKAGE_DIR = REPO_ROOT / "src" / "epcsaft"


def _run(cmd: list[str], *, env: dict[str, str]) -> None:
    print("Running:", " ".join(cmd), flush=True)
    subprocess.run(cmd, cwd=str(REPO_ROOT), env=env, check=True)


def _capture(cmd: list[str], *, env: dict[str, str]) -> str:
    return subprocess.check_output(cmd, cwd=str(REPO_ROOT), env=env, text=True).strip()


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
    cache = BUILD_DIR / "CMakeCache.txt"
    if not cache.exists():
        return None
    for line in cache.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith("CMAKE_GENERATOR:INTERNAL="):
            return line.split("=", 1)[1].strip()
    return None


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


def _configure(env: dict[str, str]) -> None:
    pybind11_dir = _capture([sys.executable, "-m", "pybind11", "--cmakedir"], env=env)
    cmd = [
        "cmake",
        "-S",
        str(REPO_ROOT),
        "-B",
        str(BUILD_DIR),
        "-DEPCSAFT_DEV_INPLACE=ON",
        "-DCMAKE_TRY_COMPILE_TARGET_TYPE=STATIC_LIBRARY",
        f"-DPython_EXECUTABLE={sys.executable}",
        f"-Dpybind11_DIR={pybind11_dir}",
    ]
    cmd.extend(_generator_args(env, _configured_generator()))
    _run(cmd, env=env)


def _build(env: dict[str, str], parallel: str | None) -> None:
    cmd = ["cmake", "--build", str(BUILD_DIR)]
    if parallel:
        cmd.extend(["--parallel", parallel])
    _run(cmd, env=env)


def _timed(label: str, action) -> float:
    start = time.perf_counter()
    action()
    elapsed = time.perf_counter() - start
    print(f"Timing: {label} completed in {elapsed:.2f}s", flush=True)
    return elapsed


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Build epcsaft._core in-place with direct CMake.")
    parser.add_argument("--clean", action="store_true", help="Remove build/dev and any in-place _core artifact first.")
    parser.add_argument(
        "--configure-only", action="store_true", help="Configure the CMake dev build tree without building."
    )
    parser.add_argument(
        "--build-only", action="store_true", help="Build the existing CMake dev build tree without reconfiguring."
    )
    parser.add_argument("--parallel", help="Optional CMake build parallelism value.")
    parser.add_argument(
        "--generator",
        choices=("auto", "ninja", "mingw"),
        default="auto",
        help="CMake generator for a new configure. Auto prefers Ninja when available.",
    )
    return parser


def main() -> int:
    parser = _parser()
    args = parser.parse_args()
    if args.clean and args.build_only:
        parser.error("--clean cannot be combined with --build-only")
    if args.configure_only and args.build_only:
        parser.error("--configure-only cannot be combined with --build-only")

    total_start = time.perf_counter()
    if args.clean:
        _timed("clean", _clean)

    env = _env()
    if args.generator != "auto":
        env["EPCSAFT_CMAKE_GENERATOR"] = args.generator
    if not args.build_only:
        _timed("configure", lambda: _configure(env))
    else:
        print("Timing: configure skipped (--build-only)", flush=True)
    if not args.configure_only:
        _timed("build", lambda: _build(env, args.parallel))
    else:
        print("Timing: build skipped (--configure-only)", flush=True)
    print(f"Timing: total completed in {time.perf_counter() - total_start:.2f}s", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
