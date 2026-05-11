from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest
import tomllib

REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPT_PATH = REPO_ROOT / "scripts" / "build_epcsaft.py"


def _load_script():
    spec = importlib.util.spec_from_file_location("build_epcsaft_for_test", SCRIPT_PATH)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_build_script_preserves_existing_generator_for_auto() -> None:
    build = _load_script()

    args = build._generator_args({"EPCSAFT_CMAKE_GENERATOR": ""}, configured_generator="MinGW Makefiles")

    assert args == []


def test_build_script_rejects_explicit_generator_switch_without_clean() -> None:
    build = _load_script()

    with pytest.raises(RuntimeError, match="Use --clean before switching"):
        build._generator_args({"EPCSAFT_CMAKE_GENERATOR": "ninja"}, configured_generator="MinGW Makefiles")


def test_build_script_auto_prefers_ninja_for_new_build_tree(monkeypatch) -> None:
    build = _load_script()

    monkeypatch.setattr(build.shutil, "which", lambda name, path=None: "ninja.exe" if name == "ninja" else None)

    args = build._generator_args({"EPCSAFT_CMAKE_GENERATOR": ""}, configured_generator=None)

    assert args == ["-G", "Ninja"]


def test_build_script_reads_version_from_pyproject() -> None:
    build = _load_script()
    pyproject = tomllib.loads((REPO_ROOT / "pyproject.toml").read_text(encoding="utf-8"))

    assert build._pyproject_version() == pyproject["project"]["version"]
