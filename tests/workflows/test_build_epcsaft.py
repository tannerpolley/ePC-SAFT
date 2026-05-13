from __future__ import annotations

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
