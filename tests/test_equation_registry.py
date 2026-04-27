from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest
from scripts import sync_equation_registry


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


def test_traceability_report_flags_missing_cpp_refs_except_documentation_only() -> None:
    entries = [
        {"eqid": "implemented_without_owner", "status": "Implemented", "cpp_refs": []},
        {"eqid": "documentation_helper", "status": "Documentation-only", "cpp_refs": []},
        {"eqid": "implemented_with_owner", "status": "Implemented", "cpp_refs": [{"file": "x.cpp"}]},
    ]

    missing = sync_equation_registry.missing_cpp_ref_entries(entries)
    report = sync_equation_registry.render_traceability_report(missing)

    assert [entry["eqid"] for entry in missing] == ["implemented_without_owner"]
    assert "implemented_without_owner" in report
    assert "documentation_helper" not in report
