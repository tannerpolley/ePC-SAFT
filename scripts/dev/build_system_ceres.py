from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
CERES_VERSION = "2.2.0"
CERES_REPOSITORY = "https://github.com/ceres-solver/ceres-solver.git"
DEFAULT_ROOT = REPO_ROOT / "build" / "system-ceres" / CERES_VERSION


def _run(cmd: list[str], *, cwd: Path | None = None, env: dict[str, str] | None = None) -> None:
    print("Running:", " ".join(cmd), flush=True)
    subprocess.run(cmd, cwd=str(cwd or REPO_ROOT), env=env, check=True)


def _capture(cmd: list[str], *, cwd: Path | None = None, env: dict[str, str] | None = None) -> str:
    return subprocess.check_output(cmd, cwd=str(cwd or REPO_ROOT), env=env, text=True).strip()


def _env() -> dict[str, str]:
    env = os.environ.copy()
    temp_root = REPO_ROOT / "build" / "temp"
    temp_root.mkdir(parents=True, exist_ok=True)
    env["TMP"] = str(temp_root.resolve())
    env["TEMP"] = str(temp_root.resolve())
    env["TMPDIR"] = str(temp_root.resolve())
    return env


def _generator_args(env: dict[str, str], requested: str) -> list[str]:
    if requested == "ninja":
        return ["-G", "Ninja"]
    if requested == "mingw":
        return ["-G", "MinGW Makefiles"]
    if shutil.which("ninja", path=env.get("PATH")):
        return ["-G", "Ninja"]
    if os.name == "nt" and shutil.which("mingw32-make", path=env.get("PATH")):
        return ["-G", "MinGW Makefiles"]
    return []


def _write_eigen_config(config_dir: Path, env: dict[str, str]) -> None:
    include_dir = _capture([sys.executable, "-c", "import includeigen; print(includeigen.get_include())"], env=env)
    include_cmake = include_dir.replace("\\", "/")
    config_dir.mkdir(parents=True, exist_ok=True)
    (config_dir / "Eigen3Config.cmake").write_text(
        "\n".join(
            [
                "if(NOT TARGET Eigen3::Eigen)",
                "    add_library(Eigen3::Eigen INTERFACE IMPORTED)",
                f'    set_target_properties(Eigen3::Eigen PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "{include_cmake}")',
                "endif()",
                "set(Eigen3_VERSION 3.4.0)",
                "set(Eigen3_FOUND TRUE)",
                "",
            ]
        ),
        encoding="utf-8",
    )
    (config_dir / "Eigen3ConfigVersion.cmake").write_text(
        "\n".join(
            [
                "set(PACKAGE_VERSION 3.4.0)",
                "set(PACKAGE_VERSION_COMPATIBLE TRUE)",
                "set(PACKAGE_VERSION_EXACT FALSE)",
                "",
            ]
        ),
        encoding="utf-8",
    )


def _clone_ceres(source_dir: Path, env: dict[str, str]) -> None:
    if source_dir.exists():
        return
    source_dir.parent.mkdir(parents=True, exist_ok=True)
    _run(
        [
            "git",
            "clone",
            "--depth",
            "1",
            "--branch",
            CERES_VERSION,
            CERES_REPOSITORY,
            str(source_dir),
        ],
        env=env,
    )


def _ceres_config_dir(install_dir: Path) -> Path:
    candidates = [
        install_dir / "lib" / "cmake" / "Ceres",
        install_dir / "lib64" / "cmake" / "Ceres",
    ]
    for candidate in candidates:
        if (candidate / "CeresConfig.cmake").is_file():
            return candidate
    return candidates[0]


def _configure_ceres(root: Path, generator: str, env: dict[str, str]) -> None:
    source_dir = root / "src"
    build_dir = root / "build"
    install_dir = root / "install"
    eigen_config_dir = root / "generated_eigen3_config"
    _write_eigen_config(eigen_config_dir, env)
    cmd = [
        "cmake",
        "-S",
        str(source_dir),
        "-B",
        str(build_dir),
        "-DCMAKE_BUILD_TYPE=Release",
        "-DCMAKE_CXX_STANDARD=17",
        "-DCMAKE_CXX_STANDARD_REQUIRED=ON",
        "-DCMAKE_CXX_EXTENSIONS=OFF",
        "-DCMAKE_POSITION_INDEPENDENT_CODE=ON",
        f"-DCMAKE_INSTALL_PREFIX={install_dir}",
        f"-DEigen3_DIR={eigen_config_dir}",
        "-DBUILD_TESTING=OFF",
        "-DBUILD_EXAMPLES=OFF",
        "-DBUILD_BENCHMARKS=OFF",
        "-DMINIGLOG=ON",
        "-DGFLAGS=OFF",
        "-DSUITESPARSE=OFF",
        "-DEIGENSPARSE=OFF",
        "-DLAPACK=OFF",
        "-DUSE_CUDA=OFF",
        "-DSCHUR_SPECIALIZATIONS=OFF",
        "-DBUILD_SHARED_LIBS=OFF",
        "-DEXPORT_BUILD_DIR=ON",
    ]
    if os.name == "nt":
        cmd.append("-DCMAKE_CXX_FLAGS=-include cstdint")
    cmd.extend(_generator_args(env, generator))
    _run(cmd, env=env)


def _build_and_install(root: Path, parallel: str, env: dict[str, str]) -> None:
    _run(["cmake", "--build", str(root / "build"), "--target", "install", "--parallel", parallel], env=env)


def _print_usage(root: Path) -> None:
    config_dir = _ceres_config_dir(root / "install")
    print(f"CeresConfigDir: {config_dir}")
    print("PowerShell:")
    print(f'  $env:EPCSAFT_PEP517_CERES_DIR = "{config_dir}"')
    print('  $env:EPCSAFT_PEP517_USE_SYSTEM_CERES = "1"')
    print('  $env:EPCSAFT_PEP517_BUILD_DIR = "$PWD\\.uv-cache\\epcsaft-build"')
    print("  uv sync --reinstall-package epcsaft")


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Build and install a reusable local Ceres package for ePC-SAFT.")
    parser.add_argument("--root", type=Path, default=DEFAULT_ROOT, help="Reusable Ceres build/install root.")
    parser.add_argument("--parallel", default="4", help="CMake build parallelism.")
    parser.add_argument("--generator", choices=("auto", "ninja", "mingw"), default="auto")
    parser.add_argument("--configure-only", action="store_true")
    parser.add_argument("--build-only", action="store_true")
    parser.add_argument("--print-env", action="store_true", help="Print package-install environment variables.")
    return parser


def main() -> int:
    args = _parser().parse_args()
    if args.configure_only and args.build_only:
        raise SystemExit("--configure-only cannot be combined with --build-only")
    root = args.root.expanduser().resolve()
    env = _env()
    if args.print_env:
        _print_usage(root)
        return 0
    if not args.build_only:
        _clone_ceres(root / "src", env)
        _configure_ceres(root, args.generator, env)
    if not args.configure_only:
        _build_and_install(root, args.parallel, env)
    _print_usage(root)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
