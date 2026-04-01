import argparse
import importlib.machinery
import importlib.util
import subprocess
import sys
from pathlib import Path
from typing import Literal


REPO_ROOT = Path(__file__).resolve().parents[1]
PACKAGE_ROOT = REPO_ROOT / "src" / "pcsaft"
EXPECTED_PACKAGE_INIT = (PACKAGE_ROOT / "__init__.py").resolve()
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


def _editable_install_matches_repo() -> bool:
    spec = importlib.util.find_spec("pcsaft")
    if spec is None or spec.origin is None:
        return False

    module_path = Path(spec.origin).resolve()
    return module_path == EXPECTED_PACKAGE_INIT


def _rebuild_plan() -> tuple[Literal["build_ext", "editable_install"], str] | None:
    editable_matches_repo = _editable_install_matches_repo()
    if not editable_matches_repo:
        return "editable_install", "pcsaft is not installed editable from this repo in the active environment"

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
        if newest_input in INPLACE_BUILD_INPUTS and editable_matches_repo:
            return "build_ext", reason
        return "editable_install", reason

    return None


def _run_command(cmd: list[str]) -> None:
    print("Running:", " ".join(cmd))
    subprocess.run(cmd, cwd=str(REPO_ROOT), check=True)


def _run_build_ext() -> None:
    _run_command([sys.executable, "setup.py", "build_ext", "--inplace"])


def _run_editable_install() -> None:
    _run_command([sys.executable, "-m", "pip", "install", "-e", ".", "--no-build-isolation", "--no-deps"])


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
        print("Editable install is up to date; skipping rebuild.")
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
