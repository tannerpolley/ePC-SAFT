from __future__ import annotations

import importlib.util
import json
import shutil
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = REPO_ROOT / "src"
PLOT_MANIFEST = REPO_ROOT / "docs" / "plots" / "manifest.json"
STALE_TRACKED_REPORTS = (
    REPO_ROOT / "scripts" / "paper_validation" / "tools" / "out" / "plot_asset_report.csv",
    REPO_ROOT / "tests" / "plots" / "out" / "plot_asset_report.csv",
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


def _tool_path(name: str) -> str:
    return shutil.which(name) or "<missing>"


def _tracked_generated_count() -> int | None:
    output = _git_output("ls-files")
    if output is None:
        return None
    count = 0
    generated_extensions = {".png", ".svg", ".csv", ".md"}
    for raw_line in output.splitlines():
        path = Path(raw_line)
        parts = path.parts
        if "out" not in parts:
            continue
        if not parts or parts[0] not in {"scripts", "tests", "src"}:
            continue
        if path.suffix.lower() in generated_extensions:
            count += 1
    return count


def _manifest_state() -> tuple[str, str]:
    if not PLOT_MANIFEST.exists():
        return "missing", "uv run python scripts\\build_plot_manifest.py --refresh"
    try:
        payload = json.loads(PLOT_MANIFEST.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        return "invalid-json", "uv run python scripts\\build_plot_manifest.py --refresh"
    items = payload.get("items")
    if not isinstance(items, list):
        return "invalid-schema", "uv run python scripts\\build_plot_manifest.py --refresh"
    return f"current ({len(items)} item(s))", "none"


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
    print(f"ninja: {_tool_path('ninja')}")
    print(f"epcsaft_import: {package_path if package_path else '<missing>'}")
    print(f"epcsaft_import_error: {package_error or '<none>'}")
    print(f"epcsaft_core: {core_path if core_path else '<missing>'}")
    print(f"epcsaft_core_error: {core_error or '<none>'}")
    manifest_state, manifest_next = _manifest_state()
    print(f"plot_manifest: {manifest_state}")
    print(f"stale_generated_reports: {_stale_report_state()}")
    tracked_generated = _tracked_generated_count()
    print(f"tracked_source_out_generated_files: {tracked_generated if tracked_generated is not None else '<unknown>'}")

    if package_path is None:
        print("install_state: missing-package")
        print("next_command: uv sync --no-install-project")
        return 1
    if core_path is None:
        print("install_state: missing-core")
        print("next_command: uv run python scripts\\build_epcsaft.py")
        return 1
    if manifest_next != "none":
        print("install_state: missing-plot-manifest")
        print(f"next_command: {manifest_next}")
        return 1

    print("install_state: current")
    print("next_command: none")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
