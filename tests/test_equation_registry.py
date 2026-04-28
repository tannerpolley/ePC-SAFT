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


def test_equation_registry_strict_traceability_passes_current_registry() -> None:
    if not TEX_PATH.exists():
        pytest.skip("docs/latex/equations.tex is unavailable because the docs submodule is not checked out")

    result = subprocess.run(
        [sys.executable, "scripts/sync_equation_registry.py", "--check", "--strict-traceability"],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, result.stdout + result.stderr


def test_traceability_report_flags_missing_cpp_refs_except_documentation_only() -> None:
    entries = [
        {"eqid": "implemented_without_owner", "section": "Runtime", "status": "Implemented", "cpp_refs": []},
        {"eqid": "other_missing_owner", "section": "Parameter Setup", "status": "Implemented", "cpp_refs": []},
        {"eqid": "documentation_helper", "section": "Runtime", "status": "Documentation-only", "cpp_refs": []},
        {"eqid": "implemented_with_owner", "section": "Runtime", "status": "Implemented", "cpp_refs": [{"file": "x.cpp"}]},
    ]

    missing = sync_equation_registry.missing_cpp_ref_entries(entries)
    report = sync_equation_registry.render_traceability_report(missing)

    assert [entry["eqid"] for entry in missing] == ["implemented_without_owner", "other_missing_owner"]
    assert "2 EqIDs without C++ owner comments" in report
    assert "Parameter Setup: 1" in report
    assert "Runtime: 1" in report
    assert "implemented_without_owner" in report
    assert "other_missing_owner" in report
    assert "documentation_helper" not in report


def test_strict_traceability_allows_documentation_only_entries() -> None:
    entries = [
        {"eqid": "documentation_helper", "status": "Documentation-only", "cpp_refs": []},
        {"eqid": "docs_helper", "status": "docs_only", "cpp_refs": []},
        {"eqid": "reference_helper", "status": "reference-only", "cpp_refs": []},
        {"eqid": "implemented_with_owner", "status": "Implemented", "cpp_refs": [{"file": "x.cpp"}]},
    ]

    sync_equation_registry.enforce_traceability(entries)


def test_strict_traceability_fails_missing_implemented_owner() -> None:
    entries = [
        {"eqid": "implemented_without_owner", "section": "Runtime", "status": "Implemented", "cpp_refs": []},
    ]

    with pytest.raises(SystemExit) as excinfo:
        sync_equation_registry.enforce_traceability(entries)

    assert "implemented_without_owner" in str(excinfo.value)


def test_validate_links_fails_unknown_cpp_eqid() -> None:
    entries = [{"eqid": "known_equation"}]
    code_refs = {"known_equation": [], "unknown_equation": [{"file": "x.cpp"}]}

    with pytest.raises(ValueError, match="unknown_equation"):
        sync_equation_registry.validate_links(entries, code_refs)


def test_sync_equation_registry_help_includes_strict_traceability_flag() -> None:
    result = subprocess.run(
        [sys.executable, "scripts/sync_equation_registry.py", "--help"],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
        check=False,
    )

    assert result.returncode == 0, result.stdout + result.stderr
    assert "--strict-traceability" in result.stdout
