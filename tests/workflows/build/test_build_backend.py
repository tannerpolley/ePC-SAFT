from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[3]
BACKEND_PATH = REPO_ROOT / "build_backend" / "epcsaft_build_backend.py"


def _load_backend():
    spec = importlib.util.spec_from_file_location("epcsaft_build_backend_for_test", BACKEND_PATH)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_pep517_build_backend_uses_isolated_build_dir_by_default(tmp_path, monkeypatch) -> None:
    backend = _load_backend()
    monkeypatch.setenv("TMP", str(tmp_path))
    monkeypatch.setenv("TEMP", str(tmp_path))
    monkeypatch.setenv("TMPDIR", str(tmp_path))

    config = backend._isolated_build_config(None)

    build_dir = Path(config["build-dir"])
    assert build_dir.exists()
    assert REPO_ROOT not in build_dir.parents
    assert build_dir.name.startswith("epcsaft-pep517-build-")


def test_pep517_build_backend_preserves_explicit_build_dir(tmp_path) -> None:
    backend = _load_backend()
    requested = tmp_path / "build"

    config = backend._isolated_build_config({"build-dir": str(requested)})

    assert config["build-dir"] == str(requested)


def test_pep517_build_backend_honors_persistent_build_dir_env(tmp_path, monkeypatch) -> None:
    backend = _load_backend()
    requested = tmp_path / "persistent"
    monkeypatch.setenv("EPCSAFT_PEP517_BUILD_DIR", str(requested))

    config = backend._isolated_build_config(None)

    assert Path(config["build-dir"]) == requested.resolve()
    assert requested.exists()


def test_pep517_build_backend_uses_prebuilt_ceres_env(tmp_path, monkeypatch) -> None:
    backend = _load_backend()
    ceres_dir = tmp_path / "ceres" / "lib" / "cmake" / "Ceres"
    ceres_dir.mkdir(parents=True)
    (ceres_dir / "CeresConfig.cmake").write_text("# test config\n", encoding="utf-8")
    monkeypatch.setenv("EPCSAFT_PEP517_CERES_DIR", str(ceres_dir))

    config = backend._isolated_build_config(None)

    assert config["cmake.define.EPCSAFT_ENABLE_CERES"] == "ON"
    assert config["cmake.define.EPCSAFT_USE_SYSTEM_CERES"] == "ON"
    assert Path(config["cmake.define.Ceres_DIR"]) == ceres_dir.resolve()


def test_pep517_build_backend_does_not_disable_ceres_by_default(monkeypatch) -> None:
    backend = _load_backend()
    monkeypatch.delenv("EPCSAFT_PEP517_CERES_DIR", raising=False)
    monkeypatch.delenv("EPCSAFT_PEP517_USE_SYSTEM_CERES", raising=False)
    monkeypatch.delenv("Ceres_DIR", raising=False)

    config = backend._isolated_build_config(None)

    assert "cmake.define.EPCSAFT_ENABLE_CERES" not in config
    assert "cmake.define.EPCSAFT_USE_SYSTEM_CERES" not in config
    assert "cmake.define.Ceres_DIR" not in config


def test_pep517_build_backend_rejects_bad_ceres_dir_env(tmp_path, monkeypatch) -> None:
    backend = _load_backend()
    bad_dir = tmp_path / "missing-config"
    bad_dir.mkdir()
    monkeypatch.setenv("EPCSAFT_PEP517_CERES_DIR", str(bad_dir))

    with pytest.raises(FileNotFoundError) as excinfo:
        backend._isolated_build_config(None)
    assert "CeresConfig.cmake" in str(excinfo.value)


def test_pep660_editable_hook_uses_isolated_build_dir(tmp_path, monkeypatch) -> None:
    backend = _load_backend()
    wheel_dir = tmp_path / "wheelhouse"
    wheel_dir.mkdir()
    captured = {}

    def fake_build_editable(wheel_directory, config_settings=None, metadata_directory=None):
        captured["wheel_directory"] = wheel_directory
        captured["config_settings"] = config_settings
        captured["metadata_directory"] = metadata_directory
        return "epcsaft-1.5.2-0.editable-py3-none-any.whl"

    monkeypatch.setattr(backend._scikit_build, "build_editable", fake_build_editable)

    result = backend.build_editable(str(wheel_dir), config_settings={"editable.mode": "redirect"})

    assert result == "epcsaft-1.5.2-0.editable-py3-none-any.whl"
    assert captured["wheel_directory"] == str(wheel_dir)
    assert captured["metadata_directory"] is None
    assert captured["config_settings"]["editable.mode"] == "redirect"
    build_dir = Path(captured["config_settings"]["build-dir"])
    assert build_dir.exists()
    assert REPO_ROOT not in build_dir.parents
