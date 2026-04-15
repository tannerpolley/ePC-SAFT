from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest


REPO_ROOT = Path(__file__).resolve().parents[1]
TEX_PATH = REPO_ROOT / "docs" / "latex" / "equations.tex"


def test_equation_registry_outputs_are_synced() -> None:
    if not TEX_PATH.exists():
        pytest.skip("docs/latex/equations.tex is unavailable because the docs submodule is not checked out")

    result = subprocess.run(
        [sys.executable, "scripts/sync_equation_registry.py", "--check"],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, result.stdout + result.stderr
