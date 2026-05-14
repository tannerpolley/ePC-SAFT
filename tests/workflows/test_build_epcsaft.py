from __future__ import annotations

import json

import pytest

from scripts import build_epcsaft


def test_build_script_rejects_clean_build_only_combination(monkeypatch, capsys) -> None:
    monkeypatch.setattr("sys.argv", ["build_epcsaft.py", "--clean", "--build-only"])

    with pytest.raises(SystemExit) as excinfo:
        build_epcsaft.main()

    assert excinfo.value.code == 2
    assert "--clean cannot be combined with --build-only" in capsys.readouterr().err


def test_build_script_help_lists_incremental_workflow_flags(capsys) -> None:
    parser = build_epcsaft._parser()

    parser.print_help()
    help_text = capsys.readouterr().out

    assert "--configure-only" in help_text
    assert "--build-only" in help_text
    assert "--parallel" in help_text
    assert "--status" in help_text


def test_build_script_defaults_enable_ceres_and_cppad() -> None:
    args = build_epcsaft._parser().parse_args([])

    assert args.enable_ceres is True
    assert args.enable_cppad is True


def test_build_script_disable_flags_turn_off_optional_native_dependencies() -> None:
    args = build_epcsaft._parser().parse_args(["--disable-ceres", "--disable-cppad"])

    assert args.enable_ceres is False
    assert args.enable_cppad is False


def test_cmake_and_presets_match_cli_default_dependency_state() -> None:
    cmake_text = (build_epcsaft.REPO_ROOT / "CMakeLists.txt").read_text(encoding="utf-8")
    presets = json.loads((build_epcsaft.REPO_ROOT / "CMakePresets.json").read_text(encoding="utf-8"))

    assert 'option(EPCSAFT_ENABLE_CERES "Enable Ceres Solver support for native equilibrium nonlinear solves" ON)' in cmake_text
    assert 'option(EPCSAFT_ENABLE_CPPAD "Enable package-wide CppAD support" ON)' in cmake_text

    native_default = next(p for p in presets["configurePresets"] if p["name"] == "native-default")
    assert native_default["cacheVariables"]["EPCSAFT_ENABLE_CERES"] == "ON"
    assert native_default["cacheVariables"]["EPCSAFT_ENABLE_CPPAD"] == "ON"

    dev_mingw = next(p for p in presets["configurePresets"] if p["name"] == "dev-mingw")
    assert dev_mingw["cacheVariables"]["EPCSAFT_ENABLE_CERES"] == "OFF"
    assert dev_mingw["cacheVariables"]["EPCSAFT_ENABLE_CPPAD"] == "OFF"


def test_build_status_reports_generator_core_optional_flags_and_stale_lock(tmp_path, monkeypatch) -> None:
    build_dir = tmp_path / "build" / "dev"
    package_dir = tmp_path / "src" / "epcsaft"
    build_dir.mkdir(parents=True)
    package_dir.mkdir(parents=True)
    (package_dir / "_core.cp313-win_amd64.pyd").write_bytes(b"native")
    (build_dir / ".ninja_lock").write_text("", encoding="utf-8")
    (build_dir / ".ninja_log").write_text(
        "# ninja log v7\n1\t2\t3\tCMakeFiles/example.cpp.obj\tabc\n",
        encoding="utf-8",
    )
    (build_dir / "CMakeCache.txt").write_text(
        "\n".join(
            [
                "CMAKE_GENERATOR:INTERNAL=Ninja",
                "EPCSAFT_ENABLE_CERES:BOOL=OFF",
                "EPCSAFT_ENABLE_CPPAD:BOOL=ON",
            ]
        ),
        encoding="utf-8",
    )

    monkeypatch.setattr(build_epcsaft, "BUILD_DIR", build_dir)
    monkeypatch.setattr(build_epcsaft, "PACKAGE_DIR", package_dir)
    monkeypatch.setattr(build_epcsaft, "_repo_build_processes", lambda: [])

    lines = build_epcsaft._status_lines(stale_lock_seconds=0)

    assert "configured_generator: Ninja" in lines
    assert "native_core: present" in lines
    assert "ceres_configured: OFF" in lines
    assert "cppad_configured: ON" in lines
    assert "ninja_lock: present" in lines
    assert "stale_ninja_lock: true" in lines
    assert "last_ninja_target: CMakeFiles/example.cpp.obj" in lines
