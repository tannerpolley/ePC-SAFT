import argparse
import importlib.machinery
import importlib.util
import os
import subprocess
import sys
from pathlib import Path
from typing import Literal


REPO_ROOT = Path(__file__).resolve().parents[1]
PACKAGE_ROOT = REPO_ROOT / "src" / "pcsaft"
EXPECTED_PACKAGE_INIT = (PACKAGE_ROOT / "__init__.py").resolve()
BUILD_ROOT = REPO_ROOT / "build"
PIP_TEMP_ROOT = BUILD_ROOT / "codex-pip-temp"
PIP_CACHE_ROOT = BUILD_ROOT / "codex-pip-cache"
COMPILED_INPUTS = [
    REPO_ROOT / "pyproject.toml",
    REPO_ROOT / "setup.py",
    REPO_ROOT / "MANIFEST.in",
    PACKAGE_ROOT / "pcsaft.pyx",
    PACKAGE_ROOT / "pcsaft.pxd",
    PACKAGE_ROOT / "pcsaft_electrolyte.cpp",
    PACKAGE_ROOT / "pcsaft_electrolyte.h",
]
INPLACE_BUILD_INPUTS = {
    REPO_ROOT / "setup.py",
    PACKAGE_ROOT / "pcsaft.pyx",
    PACKAGE_ROOT / "pcsaft.pxd",
    PACKAGE_ROOT / "pcsaft_electrolyte.cpp",
    PACKAGE_ROOT / "pcsaft_electrolyte.h",
}


def _cleanup_generated_cpp() -> None:
    generated_cpp = PACKAGE_ROOT / "pcsaft.cpp"
    if generated_cpp.exists():
        generated_cpp.unlink()


def _installed_module_path() -> Path | None:
    spec = importlib.util.find_spec("pcsaft")
    if spec is None or spec.origin is None:
        return None
    return Path(spec.origin).resolve()


def _repo_root_from_module_path(module_path: Path) -> Path | None:
    if module_path.name != "__init__.py":
        return None
    package_root = module_path.parent
    if package_root.name != "pcsaft":
        return None
    src_root = package_root.parent
    if src_root.name != "src":
        return None
    return src_root.parent


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


def _editable_install_state() -> tuple[Literal["exact", "compatible_checkout", "mismatch"], str]:
    module_path = _installed_module_path()
    if module_path is None:
        return "mismatch", "pcsaft is not importable in the active environment"

    if module_path == EXPECTED_PACKAGE_INIT:
        return "exact", "pcsaft is installed editable from this repo"

    installed_repo_root = _repo_root_from_module_path(module_path)
    if installed_repo_root is None:
        return "mismatch", f"pcsaft resolves to {module_path}, not to a src/pcsaft checkout"

    current_head = _git_head(REPO_ROOT)
    installed_head = _git_head(installed_repo_root)
    if current_head is not None and installed_head is not None and current_head == installed_head:
        short_head = current_head[:12]
        return (
            "compatible_checkout",
            f"pcsaft is installed from {installed_repo_root}, but both checkouts share git HEAD {short_head}",
        )

    return "mismatch", f"pcsaft is installed from {installed_repo_root}, not this repo"


def _extension_candidates() -> list[Path]:
    candidates: list[Path] = []
    for suffix in importlib.machinery.EXTENSION_SUFFIXES:
        candidates.extend(PACKAGE_ROOT.glob(f"pcsaft*{suffix}"))
    return sorted({path.resolve() for path in candidates if path.is_file()})


def _latest_extension_mtime() -> float | None:
    candidates = _extension_candidates()
    if not candidates:
        return None
    return max(path.stat().st_mtime for path in candidates)


def _rebuild_plan() -> tuple[Literal["build_ext", "editable_install"], str] | None:
    install_state, install_reason = _editable_install_state()
    if install_state == "mismatch":
        return "editable_install", install_reason
    if install_state == "compatible_checkout":
        return None

    extension_mtime = _latest_extension_mtime()
    if extension_mtime is None:
        return "build_ext", "compiled extension artifact is missing under src/pcsaft"

    latest_input = max(path.stat().st_mtime for path in COMPILED_INPUTS if path.exists())
    if latest_input > extension_mtime:
        newest_input = max(
            (path for path in COMPILED_INPUTS if path.exists()),
            key=lambda path: path.stat().st_mtime,
        )
        reason = f"compiled build inputs are newer than the extension ({newest_input.name})"
        if newest_input in INPLACE_BUILD_INPUTS and install_state == "exact":
            return "build_ext", reason
        return "editable_install", reason

    return None


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


def _run_build_ext() -> None:
    _run_command([sys.executable, "setup.py", "build_ext", "--inplace"])


def _run_editable_install() -> None:
    _run_command(
        [sys.executable, "-m", "pip", "install", "-e", ".", "--no-build-isolation", "--no-deps"],
        env=_pip_environment(),
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

    rebuild_plan = ("build_ext", "forced rebuild requested") if args.force else _rebuild_plan()
    if rebuild_plan is None:
        _, install_reason = _editable_install_state()
        print(f"No rebuild required: {install_reason}.")
        _cleanup_generated_cpp()
        return 0

    action, rebuild_reason = rebuild_plan
    if action == "build_ext":
        print(f"Rebuilding compiled extension in place because {rebuild_reason}.")
        _run_build_ext()
    else:
        print(f"Rebuilding editable install because {rebuild_reason}.")
        _run_editable_install()

    _cleanup_generated_cpp()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
