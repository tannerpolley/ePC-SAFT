from __future__ import annotations

import argparse
import importlib
import sys
from pathlib import Path
from types import ModuleType


def _load_modules() -> tuple[ModuleType, ModuleType]:
    package = importlib.import_module("epcsaft")
    core = importlib.import_module("epcsaft._core")
    return package, core


def _failure_message(exc: BaseException) -> str:
    return (
        "status: error\n"
        f"epcsaft._core import failed: {exc}\n"
        "source checkout hint: run `uv run python scripts/build_epcsaft.py`\n"
        "installed package hint: reinstall from a wheel or rebuild from source with a working C++ toolchain"
    )


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Report epcsaft package and native extension status.")
    parser.parse_args(argv)
    try:
        package, core = _load_modules()
    except Exception as exc:
        print(_failure_message(exc))
        return 1

    package_path = Path(package.__file__).resolve()
    core_path = Path(core.__file__).resolve()
    build_info = package.runtime_build_info()
    print(f"epcsaft package: {package_path}")
    print(f"epcsaft._core: {core_path}")
    print(f"version: {build_info['package_version']}")
    print(f"source_git_commit: {build_info['source_git_commit']}")
    print("status: ok")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
