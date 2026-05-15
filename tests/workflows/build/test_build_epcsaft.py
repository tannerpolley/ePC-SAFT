from __future__ import annotations

import json

import pytest

from scripts.dev import build_epcsaft


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
    assert "--profile" in help_text
    assert "--status" in help_text


def test_build_script_default_profile_skips_ceres_and_keeps_cppad() -> None:
    args = build_epcsaft._parser().parse_args([])
    settings = build_epcsaft._resolve_settings(args)

    assert args.profile == "fast"
    assert settings.enable_ceres is False
    assert settings.enable_cppad is True
    assert settings.parallel == "2"


def test_build_script_disable_flags_turn_off_optional_native_dependencies() -> None:
    args = build_epcsaft._parser().parse_args(["--disable-ceres", "--disable-cppad"])
    settings = build_epcsaft._resolve_settings(args)

    assert settings.enable_ceres is False
    assert settings.enable_cppad is False


def test_build_script_profiles_resolve_optional_native_dependency_state() -> None:
    full = build_epcsaft._resolve_settings(build_epcsaft._parser().parse_args(["--profile", "full"]))
    cppad = build_epcsaft._resolve_settings(build_epcsaft._parser().parse_args(["--profile", "cppad"]))
    minimal = build_epcsaft._resolve_settings(build_epcsaft._parser().parse_args(["--profile", "minimal"]))
    explicit_ceres = build_epcsaft._resolve_settings(build_epcsaft._parser().parse_args(["--enable-ceres"]))
    system_ceres = build_epcsaft._resolve_settings(
        build_epcsaft._parser().parse_args(["--ceres-dir", "C:/ceres/lib/cmake/Ceres"])
    )

    assert full.enable_ceres is True
    assert full.enable_cppad is True
    assert full.parallel == "4"
    assert cppad.enable_ceres is False
    assert cppad.enable_cppad is True
    assert cppad.parallel == "2"
    assert minimal.enable_ceres is False
    assert minimal.enable_cppad is False
    assert explicit_ceres.enable_ceres is True
    assert explicit_ceres.enable_cppad is True
    assert system_ceres.enable_ceres is True
    assert system_ceres.enable_cppad is True


def test_package_defaults_keep_ceres_while_dev_script_offers_fast_profile() -> None:
    cmake_text = (build_epcsaft.REPO_ROOT / "CMakeLists.txt").read_text(encoding="utf-8")
    pyproject_text = (build_epcsaft.REPO_ROOT / "pyproject.toml").read_text(encoding="utf-8")
    presets = json.loads((build_epcsaft.REPO_ROOT / "CMakePresets.json").read_text(encoding="utf-8"))

    assert 'option(EPCSAFT_ENABLE_CERES "Enable Ceres Solver support for native equilibrium nonlinear solves" ON)' in cmake_text
    assert 'option(EPCSAFT_ENABLE_CPPAD "Enable package-wide CppAD support" ON)' in cmake_text
    assert "GIT_SHALLOW TRUE" in cmake_text
    assert "unset(Ceres_BINARY_DIR CACHE)" in cmake_text
    assert "unset(Ceres_SOURCE_DIR CACHE)" in cmake_text
    assert "EPCSAFT_ENABLE_CERES" not in pyproject_text

    native_default = next(p for p in presets["configurePresets"] if p["name"] == "native-default")
    assert native_default["cacheVariables"]["EPCSAFT_ENABLE_CERES"] == "ON"
    assert native_default["cacheVariables"]["EPCSAFT_ENABLE_CPPAD"] == "ON"

    native_cppad = next(p for p in presets["configurePresets"] if p["name"] == "native-cppad")
    assert native_cppad["cacheVariables"]["EPCSAFT_ENABLE_CERES"] == "OFF"
    assert native_cppad["cacheVariables"]["EPCSAFT_ENABLE_CPPAD"] == "ON"

    native_ceres_cppad = next(p for p in presets["configurePresets"] if p["name"] == "native-ceres-cppad")
    assert native_ceres_cppad["cacheVariables"]["EPCSAFT_ENABLE_CERES"] == "ON"
    assert native_ceres_cppad["cacheVariables"]["EPCSAFT_ENABLE_CPPAD"] == "ON"

    native_system_ceres = next(p for p in presets["configurePresets"] if p["name"] == "native-system-ceres")
    assert native_system_ceres["cacheVariables"]["EPCSAFT_ENABLE_CERES"] == "ON"
    assert native_system_ceres["cacheVariables"]["EPCSAFT_USE_SYSTEM_CERES"] == "ON"
    assert native_system_ceres["cacheVariables"]["EPCSAFT_ENABLE_CPPAD"] == "ON"

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
                "EPCSAFT_USE_SYSTEM_CERES:BOOL=OFF",
                "EPCSAFT_ENABLE_CPPAD:BOOL=ON",
                "Ceres_DIR:PATH=C:/ceres/lib/cmake/Ceres",
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
    assert "system_ceres_configured: OFF" in lines
    assert "ceres_dir: C:/ceres/lib/cmake/Ceres" in lines
    assert "cppad_configured: ON" in lines
    assert "profile_hint: fast/cppad" in lines
    assert "ninja_lock: present" in lines
    assert "stale_ninja_lock: true" in lines
    assert "last_ninja_target: CMakeFiles/example.cpp.obj" in lines
