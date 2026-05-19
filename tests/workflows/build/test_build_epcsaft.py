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


def test_build_script_default_profile_requires_ceres_and_cppad() -> None:
    args = build_epcsaft._parser().parse_args([])
    settings = build_epcsaft._resolve_settings(args)

    assert args.profile == "fast"
    assert settings.enable_ipopt is True
    assert settings.parallel == "4"


def test_build_script_can_disable_only_ipopt() -> None:
    args = build_epcsaft._parser().parse_args(["--disable-ipopt"])
    settings = build_epcsaft._resolve_settings(args)

    assert settings.enable_ipopt is False


def test_build_script_uses_local_windows_ipopt_sdk_default(monkeypatch, tmp_path) -> None:
    ipopt_root = tmp_path / "Documents" / "deps" / "ipopt-msvc"
    ipopt_root.mkdir(parents=True)
    monkeypatch.delenv("EPCSAFT_IPOPT_ROOT", raising=False)
    monkeypatch.delenv("EPCSAFT_PEP517_IPOPT_ROOT", raising=False)

    args = build_epcsaft._parser().parse_args([])
    settings = build_epcsaft._resolve_settings(args)

    assert build_epcsaft.resolve_ipopt_root_for_build(
        None,
        enable_ipopt=settings.enable_ipopt,
        ipopt_dir=args.ipopt_dir,
        default_root=ipopt_root,
        label="Ipopt root",
    ) == ipopt_root.resolve()


def test_build_script_profiles_resolve_optional_native_dependency_state() -> None:
    full = build_epcsaft._resolve_settings(build_epcsaft._parser().parse_args(["--profile", "full"]))
    fast = build_epcsaft._resolve_settings(build_epcsaft._parser().parse_args(["--profile", "fast"]))
    ipopt = build_epcsaft._resolve_settings(build_epcsaft._parser().parse_args(["--profile", "ipopt"]))
    system_ceres = build_epcsaft._resolve_settings(
        build_epcsaft._parser().parse_args(["--ceres-dir", "C:/ceres/lib/cmake/Ceres"])
    )
    system_ipopt = build_epcsaft._resolve_settings(
        build_epcsaft._parser().parse_args(["--ipopt-dir", "C:/ipopt/lib/cmake/Ipopt"])
    )

    assert full.enable_ipopt is True
    assert full.parallel == "4"
    assert fast.enable_ipopt is True
    assert fast.parallel == "4"
    assert ipopt.enable_ipopt is True
    assert ipopt.parallel == "4"
    assert system_ceres.enable_ipopt is True
    assert system_ipopt.enable_ipopt is True


def test_package_and_dev_defaults_require_ceres_and_cppad() -> None:
    cmake_text = (build_epcsaft.REPO_ROOT / "CMakeLists.txt").read_text(encoding="utf-8")
    pyproject_text = (build_epcsaft.REPO_ROOT / "pyproject.toml").read_text(encoding="utf-8")
    presets = json.loads((build_epcsaft.REPO_ROOT / "CMakePresets.json").read_text(encoding="utf-8"))

    assert 'option(EPCSAFT_ENABLE_CERES "Enable Ceres Solver support required for native regression solves" ON)' in cmake_text
    assert 'option(EPCSAFT_ENABLE_CPPAD "Enable package-wide CppAD support" ON)' in cmake_text
    assert 'option(EPCSAFT_ENABLE_IPOPT "Enable native Ipopt support for production equilibrium NLP solves" ON)' in cmake_text
    assert "native regression builds require Ceres" in cmake_text
    assert "derivative-capable package builds require CppAD" in cmake_text
    assert "GIT_SHALLOW TRUE" in cmake_text
    assert 'add_subdirectory("${ceres_solver_SOURCE_DIR}" "${ceres_solver_BINARY_DIR}" EXCLUDE_FROM_ALL)' in cmake_text
    assert 'src/epcsaft/native/equilibrium_nlp/*.cpp' in cmake_text
    assert 'src/epcsaft/native/equilibrium_nlp' in cmake_text
    assert "unset(Ceres_BINARY_DIR CACHE)" in cmake_text
    assert "unset(Ceres_SOURCE_DIR CACHE)" in cmake_text
    assert "EPCSAFT_ENABLE_CERES" not in pyproject_text

    native_default = next(p for p in presets["configurePresets"] if p["name"] == "native-default")
    assert native_default["cacheVariables"]["EPCSAFT_ENABLE_CERES"] == "ON"
    assert native_default["cacheVariables"]["EPCSAFT_ENABLE_CPPAD"] == "ON"
    assert native_default["cacheVariables"]["EPCSAFT_ENABLE_IPOPT"] == "ON"

    native_required = next(p for p in presets["configurePresets"] if p["name"] == "native-required")
    assert native_required["cacheVariables"]["EPCSAFT_ENABLE_CERES"] == "ON"
    assert native_required["cacheVariables"]["EPCSAFT_ENABLE_CPPAD"] == "ON"

    native_ceres_cppad = next(p for p in presets["configurePresets"] if p["name"] == "native-ceres-cppad")
    assert native_ceres_cppad["cacheVariables"]["EPCSAFT_ENABLE_CERES"] == "ON"
    assert native_ceres_cppad["cacheVariables"]["EPCSAFT_ENABLE_CPPAD"] == "ON"

    native_system_ceres = next(p for p in presets["configurePresets"] if p["name"] == "native-system-ceres")
    assert native_system_ceres["cacheVariables"]["EPCSAFT_ENABLE_CERES"] == "ON"
    assert native_system_ceres["cacheVariables"]["EPCSAFT_USE_SYSTEM_CERES"] == "ON"
    assert native_system_ceres["cacheVariables"]["EPCSAFT_ENABLE_CPPAD"] == "ON"

    dev_native = next(p for p in presets["configurePresets"] if p["name"] == "dev-native")
    assert dev_native["generator"] == "Ninja"
    assert dev_native["binaryDir"] == "${sourceDir}/build/dev"
    assert dev_native["cacheVariables"]["CMAKE_BUILD_TYPE"] == "Release"
    assert dev_native["cacheVariables"]["EPCSAFT_ENABLE_CERES"] == "ON"
    assert dev_native["cacheVariables"]["EPCSAFT_ENABLE_CPPAD"] == "ON"


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
                "EPCSAFT_ENABLE_CERES:BOOL=ON",
                "EPCSAFT_USE_SYSTEM_CERES:BOOL=OFF",
                "EPCSAFT_ENABLE_CPPAD:BOOL=ON",
                "EPCSAFT_ENABLE_IPOPT:BOOL=OFF",
                "EPCSAFT_USE_SYSTEM_IPOPT:BOOL=ON",
                "Ipopt_DIR:PATH=C:/ipopt/lib/cmake/Ipopt",
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
    assert "ceres_configured: ON" in lines
    assert "system_ceres_configured: OFF" in lines
    assert "ceres_dir: C:/ceres/lib/cmake/Ceres" in lines
    assert "cppad_configured: ON" in lines
    assert "ipopt_configured: OFF" in lines
    assert "system_ipopt_configured: ON" in lines
    assert "ipopt_dir: C:/ipopt/lib/cmake/Ipopt" in lines
    assert "profile_hint: fast/full" in lines
    assert "ninja_lock: present" in lines
    assert "stale_ninja_lock: true" in lines
    assert "last_ninja_target: CMakeFiles/example.cpp.obj" in lines
