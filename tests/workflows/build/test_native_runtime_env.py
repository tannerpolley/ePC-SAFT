from __future__ import annotations

import os
import types
from pathlib import Path

from scripts.dev import build_epcsaft
from scripts.dev import doctor
from scripts.dev.native_runtime_env import apply_native_runtime_env
from scripts.dev.native_runtime_env import resolve_default_windows_ipopt_sdk_root


def test_native_runtime_env_injects_ipopt_bin_from_cmake_cache(tmp_path: Path) -> None:
    ipopt_root = tmp_path / "Ipopt"
    ipopt_bin = ipopt_root / "bin"
    ipopt_bin.mkdir(parents=True)
    cache = tmp_path / "CMakeCache.txt"
    cache.write_text(
        "\n".join(
            [
                "EPCSAFT_ENABLE_IPOPT:BOOL=ON",
                f"EPCSAFT_IPOPT_ROOT:PATH={ipopt_root}",
            ]
        ),
        encoding="utf-8",
    )
    env = {"PATH": f"C:\\existing{os.pathsep}{ipopt_bin}", "EPCSAFT_RUNTIME_DLL_DIRS": "C:\\runtime"}

    runtime = apply_native_runtime_env(env, cache_path=cache)

    assert runtime.ipopt_configured is True
    assert runtime.ipopt_root == ipopt_root.resolve()
    assert runtime.ipopt_runtime_dir == ipopt_bin.resolve()
    assert runtime.applied is True
    assert env["PATH"].split(os.pathsep)[0] == str(ipopt_bin.resolve())
    assert env["PATH"].split(os.pathsep).count(str(ipopt_bin.resolve())) == 1
    assert env["EPCSAFT_RUNTIME_DLL_DIRS"].split(os.pathsep)[0] == str(ipopt_bin.resolve())


def test_native_runtime_env_leaves_non_ipopt_builds_untouched(tmp_path: Path) -> None:
    ipopt_root = tmp_path / "Ipopt"
    (ipopt_root / "bin").mkdir(parents=True)
    cache = tmp_path / "CMakeCache.txt"
    cache.write_text(
        "\n".join(
            [
                "EPCSAFT_ENABLE_IPOPT:BOOL=OFF",
                f"EPCSAFT_IPOPT_ROOT:PATH={ipopt_root}",
            ]
        ),
        encoding="utf-8",
    )
    env = {"PATH": "C:\\existing"}

    runtime = apply_native_runtime_env(env, cache_path=cache)

    assert runtime.ipopt_configured is False
    assert runtime.applied is False
    assert env == {"PATH": "C:\\existing"}


def test_native_runtime_env_resolves_local_windows_ipopt_sdk_default(tmp_path: Path) -> None:
    ipopt_root = tmp_path / "Documents" / "deps" / "ipopt-msvc"
    ipopt_root.mkdir(parents=True)

    assert resolve_default_windows_ipopt_sdk_root(home=tmp_path, platform_name="nt") == ipopt_root.resolve()
    assert resolve_default_windows_ipopt_sdk_root(home=tmp_path, platform_name="posix") is None


def test_doctor_resolves_module_paths_by_importing(monkeypatch, tmp_path: Path) -> None:
    core_path = tmp_path / "_core.pyd"

    def fake_import(name: str):
        if name == "epcsaft._core":
            raise ImportError("DLL load failed while importing _core")
        return types.SimpleNamespace(__file__=str(core_path))

    monkeypatch.setattr(doctor.importlib, "import_module", fake_import)

    path, error = doctor._module_path("epcsaft._core")

    assert path is None
    assert error == "ImportError: DLL load failed while importing _core"


def test_doctor_reports_user_local_uv_when_not_on_path(monkeypatch, tmp_path: Path) -> None:
    local_uv = tmp_path / ".local" / "bin" / "uv.exe"
    local_uv.parent.mkdir(parents=True)
    local_uv.write_text("", encoding="utf-8")
    monkeypatch.setattr(doctor.shutil, "which", lambda name: None)
    monkeypatch.setattr(doctor.sys, "platform", "win32")
    monkeypatch.setattr(doctor.Path, "home", staticmethod(lambda: tmp_path))

    assert doctor._tool_path("uv") == str(local_uv)


def test_build_status_reports_ipopt_runtime_dll_dir(monkeypatch, tmp_path: Path) -> None:
    build_dir = tmp_path / "build" / "dev"
    package_dir = tmp_path / "src" / "epcsaft"
    ipopt_root = tmp_path / "Ipopt"
    (ipopt_root / "bin").mkdir(parents=True)
    build_dir.mkdir(parents=True)
    package_dir.mkdir(parents=True)
    (build_dir / "CMakeCache.txt").write_text(
        "\n".join(
            [
                "CMAKE_GENERATOR:INTERNAL=Ninja",
                "EPCSAFT_ENABLE_CERES:BOOL=ON",
                "EPCSAFT_USE_SYSTEM_CERES:BOOL=OFF",
                "EPCSAFT_ENABLE_CPPAD:BOOL=ON",
                "EPCSAFT_ENABLE_IPOPT:BOOL=ON",
                "EPCSAFT_USE_SYSTEM_IPOPT:BOOL=ON",
                f"EPCSAFT_IPOPT_ROOT:PATH={ipopt_root}",
            ]
        ),
        encoding="utf-8",
    )

    monkeypatch.setattr(build_epcsaft, "BUILD_DIR", build_dir)
    monkeypatch.setattr(build_epcsaft, "PACKAGE_DIR", package_dir)
    monkeypatch.setattr(build_epcsaft, "_repo_build_processes", lambda: [])

    lines = build_epcsaft._status_lines()

    assert f"ipopt_runtime_dll_dir: {(ipopt_root / 'bin').resolve()}" in lines
    assert "ipopt_runtime_env_applied: true" in lines
