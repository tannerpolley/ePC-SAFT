import argparse
import importlib.util
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any, Literal


REPO_ROOT = Path(__file__).resolve().parents[1]
PACKAGE_ROOT = REPO_ROOT / "src" / "epcsaft"
NATIVE_ROOT = PACKAGE_ROOT / "native"
EXPECTED_PACKAGE_INIT = (PACKAGE_ROOT / "__init__.py").resolve()
BUILD_ROOT = REPO_ROOT / "build"
PIP_TEMP_ROOT = BUILD_ROOT / "pip-temp"
PIP_CACHE_ROOT = BUILD_ROOT / "pip-cache"
BUILD_STAMP = BUILD_ROOT / "epcsaft-editable-build.json"
COMPILED_INPUTS = [
    REPO_ROOT / "pyproject.toml",
    REPO_ROOT / "setup.py",
    REPO_ROOT / "MANIFEST.in",
    PACKAGE_ROOT / "epcsaft.pyx",
    PACKAGE_ROOT / "epcsaft.pxd",
    NATIVE_ROOT / "epcsaft_electrolyte.h",
    NATIVE_ROOT / "epcsaft_core_internal.h",
    NATIVE_ROOT / "epcsaft_autodiff_internal.h",
    NATIVE_ROOT / "third_party" / "autodiff" / "dual.hpp",
    NATIVE_ROOT / "contributions" / "epcsaft_contrib_internal.h",
    NATIVE_ROOT / "epcsaft_mixture_state.cpp",
    NATIVE_ROOT / "epcsaft_dielectric.cpp",
    NATIVE_ROOT / "contributions" / "epcsaft_contrib_hc.cpp",
    NATIVE_ROOT / "contributions" / "epcsaft_contrib_disp.cpp",
    NATIVE_ROOT / "contributions" / "epcsaft_contrib_polar.cpp",
    NATIVE_ROOT / "contributions" / "epcsaft_contrib_assoc.cpp",
    NATIVE_ROOT / "contributions" / "epcsaft_contrib_ion.cpp",
    NATIVE_ROOT / "contributions" / "epcsaft_contrib_born.cpp",
    NATIVE_ROOT / "epcsaft_ares.cpp",
    NATIVE_ROOT / "epcsaft_dadrho.cpp",
    NATIVE_ROOT / "epcsaft_dadx.cpp",
    NATIVE_ROOT / "epcsaft_dadt.cpp",
    NATIVE_ROOT / "epcsaft_hres.cpp",
    NATIVE_ROOT / "epcsaft_gres.cpp",
    NATIVE_ROOT / "epcsaft_sres.cpp",
    NATIVE_ROOT / "epcsaft_Z.cpp",
    NATIVE_ROOT / "epcsaft_mu.cpp",
    NATIVE_ROOT / "epcsaft_fugcoef.cpp",
    NATIVE_ROOT / "epcsaft_activity.cpp",
    NATIVE_ROOT / "epcsaft_density_solver.cpp",
    NATIVE_ROOT / "epcsaft_state.cpp",
]


def _cleanup_generated_cpp() -> None:
    generated_cpp = PACKAGE_ROOT / "epcsaft.cpp"
    if generated_cpp.exists():
        generated_cpp.unlink()


def _cleanup_transient_build_dirs() -> None:
    if not BUILD_ROOT.exists():
        return
    for child in BUILD_ROOT.iterdir():
        if not child.is_dir():
            continue
        if child.name.startswith(("bdist.", "lib.", "temp.")):
            shutil.rmtree(child, ignore_errors=True)


def _installed_module_path() -> Path | None:
    spec = importlib.util.find_spec("epcsaft")
    if spec is None or spec.origin is None:
        return None
    return Path(spec.origin).resolve()


def _build_signature() -> dict[str, Any]:
    signature: dict[str, Any] = {}
    head = _git_head(REPO_ROOT)
    if head is not None:
        signature["git_head"] = head

    inputs: dict[str, dict[str, int]] = {}
    for path in COMPILED_INPUTS:
        if path.exists():
            stat = path.stat()
            inputs[str(path.relative_to(REPO_ROOT))] = {
                "mtime_ns": stat.st_mtime_ns,
                "size": stat.st_size,
            }
    signature["inputs"] = inputs
    return signature


def _load_build_stamp() -> dict[str, Any] | None:
    if not BUILD_STAMP.exists():
        return None
    try:
        data = json.loads(BUILD_STAMP.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None
    return data if isinstance(data, dict) else None


def _stamp_matches() -> bool:
    stamp = _load_build_stamp()
    return stamp == _build_signature()


def _git_output(repo_root: Path, *args: str) -> str | None:
    try:
        result = subprocess.run(
            ["git", "-C", str(repo_root), *args],
            check=True,
            capture_output=True,
            text=True,
        )
    except (OSError, subprocess.CalledProcessError):
        return None
    return result.stdout.strip()


def _git_head(repo_root: Path) -> str | None:
    return _git_output(repo_root, "rev-parse", "HEAD")


def _editable_install_state() -> tuple[Literal["current", "stale", "other"], str]:
    module_path = _installed_module_path()
    if module_path is None:
        return "other", "epcsaft is not importable in the active environment"

    if module_path != EXPECTED_PACKAGE_INIT:
        return "other", f"epcsaft resolves to {module_path}, not this repo checkout"

    if _stamp_matches():
        return "current", "editable install is current for this checkout"

    return "stale", "editable build stamp is missing or out of date"


def _pip_environment() -> dict[str, str]:
    PIP_TEMP_ROOT.mkdir(parents=True, exist_ok=True)
    PIP_CACHE_ROOT.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    env["TMP"] = str(PIP_TEMP_ROOT.resolve())
    env["TEMP"] = str(PIP_TEMP_ROOT.resolve())
    env["PIP_CACHE_DIR"] = str(PIP_CACHE_ROOT.resolve())
    env.setdefault("PIP_DISABLE_PIP_VERSION_CHECK", "1")
    return env


def _run_command(cmd: list[str], env: dict[str, str] | None = None) -> None:
    print("Running:", " ".join(cmd))
    subprocess.run(cmd, cwd=str(REPO_ROOT), check=True, env=env)


def _run_editable_install() -> None:
    _run_command(
        [
            sys.executable,
            "-m",
            "pip",
            "install",
            "-e",
            ".",
            "--config-settings",
            "editable_mode=compat",
        ],
        env=_pip_environment(),
    )


def _write_build_stamp() -> None:
    BUILD_ROOT.mkdir(parents=True, exist_ok=True)
    BUILD_STAMP.write_text(
        json.dumps(_build_signature(), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--force",
        action="store_true",
        help="Reinstall the editable package even when the existing build looks current.",
    )
    args = parser.parse_args()

    pyproject = REPO_ROOT / "pyproject.toml"
    if not pyproject.exists():
        print(f"error: expected {pyproject} to exist")
        return 1

    install_state, install_reason = _editable_install_state()
    if not args.force and install_state == "current":
        print(f"No rebuild required: {install_reason}.")
        _cleanup_generated_cpp()
        _cleanup_transient_build_dirs()
        return 0

    if args.force:
        print("Reinstalling editable package because forced rebuild was requested.")
    else:
        print(f"Reinstalling editable package because {install_reason}.")

    _run_editable_install()
    _write_build_stamp()
    _cleanup_generated_cpp()
    _cleanup_transient_build_dirs()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
