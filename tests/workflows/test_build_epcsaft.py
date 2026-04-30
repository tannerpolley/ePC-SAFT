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
