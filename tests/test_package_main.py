from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

from epcsaft.__main__ import main


REPO_ROOT = Path(__file__).resolve().parents[1]


def test_python_m_epcsaft_reports_package_and_core_status(capsys) -> None:
    exit_code = main([])

    output = capsys.readouterr().out
    assert exit_code == 0
    assert "epcsaft package:" in output
    assert "epcsaft._core:" in output
    assert "status: ok" in output


def test_python_m_epcsaft_works_from_source_checkout_without_pythonpath() -> None:
    env = os.environ.copy()
    env.pop("PYTHONPATH", None)

    result = subprocess.run(
        [sys.executable, "-m", "epcsaft"],
        cwd=REPO_ROOT,
        env=env,
        capture_output=True,
        text=True,
        check=False,
    )

    assert result.returncode == 0, result.stdout + result.stderr
    assert "epcsaft package:" in result.stdout
    assert "epcsaft._core:" in result.stdout
    assert "status: ok" in result.stdout
