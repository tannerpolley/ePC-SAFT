from __future__ import annotations

import importlib.util
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_PATH = REPO_ROOT / "scripts" / "dev" / "build_system_ceres.py"


def _load_script():
    spec = importlib.util.spec_from_file_location("build_system_ceres_for_test", SCRIPT_PATH)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_system_ceres_helper_uses_build_scoped_default_root() -> None:
    script = _load_script()
    script_text = SCRIPT_PATH.read_text(encoding="utf-8")

    assert script.DEFAULT_ROOT == REPO_ROOT / "build" / "system-ceres" / script.CERES_VERSION
    assert script.CERES_VERSION == "2.2.0"
    assert "-DCMAKE_CXX_STANDARD=17" in script_text
    assert "-DCMAKE_CXX_EXTENSIONS=OFF" in script_text


def test_system_ceres_config_dir_prefers_installed_cmake_package(tmp_path) -> None:
    script = _load_script()
    config_dir = tmp_path / "install" / "lib" / "cmake" / "Ceres"
    config_dir.mkdir(parents=True)
    (config_dir / "CeresConfig.cmake").write_text("# test\n", encoding="utf-8")

    assert script._ceres_config_dir(tmp_path / "install") == config_dir


def test_system_ceres_print_env_mentions_pep517_ceres_dir(tmp_path, capsys) -> None:
    script = _load_script()
    config_dir = tmp_path / "install" / "lib" / "cmake" / "Ceres"
    config_dir.mkdir(parents=True)
    (config_dir / "CeresConfig.cmake").write_text("# test\n", encoding="utf-8")

    script._print_usage(tmp_path)

    output = capsys.readouterr().out
    assert "EPCSAFT_PEP517_CERES_DIR" in output
    assert "EPCSAFT_PEP517_BUILD_DIR" in output
    assert str(config_dir) in output
