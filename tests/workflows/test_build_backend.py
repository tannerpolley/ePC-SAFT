from __future__ import annotations

import importlib.util
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
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
