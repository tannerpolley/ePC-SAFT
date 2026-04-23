import argparse
import importlib.util
import json
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any, Literal


REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

import build_config  # noqa: E402

EXPECTED_PACKAGE_INIT = (build_config.PACKAGE_ROOT / "__init__.py").resolve()


def _cleanup_generated_cpp() -> None:
    generated_cpp = build_config.generated_cython_cpp()
    if generated_cpp.exists():
        generated_cpp.unlink()


def _cleanup_transient_build_dirs() -> None:
    if not build_config.BUILD_ROOT.exists():
        return
    for child in build_config.BUILD_ROOT.iterdir():
        if not child.is_dir():
            continue
        if child.name in {"pip-temp", "pip-cache"} or child.name.startswith(("bdist.", "lib.", "temp.")):
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
    for path in build_config.rebuild_input_paths():
        if path.exists():
            stat = path.stat()
            inputs[str(path.relative_to(REPO_ROOT))] = {
                "mtime_ns": stat.st_mtime_ns,
                "size": stat.st_size,
            }
    signature["inputs"] = inputs
    return signature


def _load_build_stamp() -> dict[str, Any] | None:
    if not build_config.BUILD_STAMP.exists():
        return None
    try:
        data = json.loads(build_config.BUILD_STAMP.read_text(encoding="utf-8"))
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


def _rebuild_plan() -> tuple[str, str] | None:
    install_state, install_reason = _editable_install_state()
    if install_state == "current":
        return None
    if install_state == "other":
        return "install-dev", install_reason
    return "build-ext-inplace", install_reason


def _run_command(cmd: list[str]) -> None:
    print("Running:", " ".join(cmd))
    subprocess.run(cmd, cwd=str(REPO_ROOT), check=True)


def _run_inplace_extension_build() -> None:
    _run_command([sys.executable, "setup.py", "build_ext", "--inplace"])


def _write_build_stamp() -> None:
    build_config.BUILD_ROOT.mkdir(parents=True, exist_ok=True)
    build_config.BUILD_STAMP.write_text(
        json.dumps(_build_signature(), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--force",
        action="store_true",
        help="Rebuild the native extension even when the existing build looks current.",
    )
    args = parser.parse_args()

    pyproject = REPO_ROOT / "pyproject.toml"
    if not pyproject.exists():
        print(f"error: expected {pyproject} to exist")
        return 1

    install_state, install_reason = _editable_install_state()
    if install_state == "other":
        print(f"error: {install_reason}.")
        print("Run: python scripts/install_dev.py")
        return 2

    if not args.force and install_state == "current":
        print(f"No rebuild required: {install_reason}.")
        _cleanup_generated_cpp()
        _cleanup_transient_build_dirs()
        return 0

    if args.force:
        print("Rebuilding in-place native extension because forced rebuild was requested.")
    else:
        print(f"Rebuilding in-place native extension because {install_reason}.")
        print("Using in-place native extension rebuild for the existing editable checkout.")

    _run_inplace_extension_build()
    _write_build_stamp()
    _cleanup_generated_cpp()
    _cleanup_transient_build_dirs()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
