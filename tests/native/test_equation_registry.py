from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest

from scripts import sync_equation_registry

REPO_ROOT = Path(__file__).resolve().parents[2]
TEX_PATH = REPO_ROOT / "docs" / "latex" / "equations.tex"

EQUATION_FAMILY_COVERAGE = {
    "density closure / root solving": [
        "test_pressure_density_edge_cases_cover_vapor_and_liquid_extremes",
        "test_pressure_density_phase_branches_do_not_cross_at_two_phase_like_state",
        "test_ionic_high_pressure_liquid_density_branch_remains_stable",
    ],
    "residual Helmholtz energy": [
        "test_native_residual_helmholtz_and_compressibility_contributions_match_neutral_contract"
    ],
    "compressibility factor": ["test_native_residual_helmholtz_and_compressibility_contributions_match_ionic_contract"],
    "temperature derivative": ["test_temperature_derivative_matches_neutral_unsupported_derivative_across_density_branches"],
    "composition derivative": ["test_composition_derivative_matches_constrained_composition_unsupported_derivative"],
    "residual chemical potential": ["test_public_methods_expose_eqid_owned_contribution_groups"],
    "fugacity coefficient": ["test_public_methods_expose_eqid_owned_contribution_groups"],
    "Debye-Huckel / Born / ionic activity coefficient": [
        "test_miac_electrolyte_variants_cover_water_nonaqueous_and_mixed_solvents",
        "test_nonionic_state_rejects_electrolyte_only_activity_methods",
    ],
    "reference-state and density cache behavior": [
        "test_runtime_cache_stats_track_density_and_reference_state_reuse",
        "test_activity_coefficient_cache_behavior_distinguishes_aux_cache_from_solvent_override",
        "test_miac_activity_cache_reuse_keeps_results_stable",
    ],
}


def test_equation_registry_outputs_are_synced() -> None:
    assert TEX_PATH.exists(), "docs/latex/equations.tex must be present as a tracked repo file"

    result = subprocess.run(
        [sys.executable, "scripts/sync_equation_registry.py", "--check"],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, result.stdout + result.stderr


def test_equation_registry_strict_traceability_passes_current_registry() -> None:
    assert TEX_PATH.exists(), "docs/latex/equations.tex must be present as a tracked repo file"

    result = subprocess.run(
        [sys.executable, "scripts/sync_equation_registry.py", "--check", "--strict-traceability"],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, result.stdout + result.stderr


def test_sync_equation_registry_missing_tracked_latex_reports_actionable_message(tmp_path, capsys) -> None:
    missing_tex = tmp_path / "docs" / "latex" / "equations.tex"

    try:
        sync_equation_registry.require_equations_tex(missing_tex)
    except SystemExit as exc:
        assert exc.code == 1
    else:
        raise AssertionError("missing equations.tex must stop registry sync")

    combined = capsys.readouterr().err
    assert "docs/latex/equations.tex" in combined
    assert "tracked repo content" in combined
    removed_workflow = "git " + "sub" + "module"
    assert removed_workflow not in combined
    assert "Traceback" not in combined


def test_traceability_report_flags_missing_cpp_refs_except_documentation_only() -> None:
    entries = [
        {"eqid": "implemented_without_owner", "section": "Runtime", "status": "Implemented", "cpp_refs": []},
        {"eqid": "other_missing_owner", "section": "Parameter Setup", "status": "Implemented", "cpp_refs": []},
        {"eqid": "documentation_helper", "section": "Runtime", "status": "Documentation-only", "cpp_refs": []},
        {
            "eqid": "implemented_with_owner",
            "section": "Runtime",
            "status": "Implemented",
            "cpp_refs": [{"file": "x.cpp"}],
        },
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


def test_docs_only_audit_lists_exempt_entries_by_section() -> None:
    entries = [
        {
            "eqid": "documentation_helper",
            "section": "Runtime",
            "status": "Documentation-only",
            "tex_file": "docs/latex/equations.tex",
            "tex_line": 10,
            "description": "notation helper",
            "cpp_refs": [],
        },
        {
            "eqid": "implemented_with_owner",
            "section": "Runtime",
            "status": "Implemented",
            "cpp_refs": [{"file": "x.cpp"}],
        },
    ]

    audit = sync_equation_registry.render_docs_only_audit(entries)

    assert "1 EqIDs are exempt from strict C++ owner enforcement" in audit
    assert "Runtime: 1" in audit
    assert "documentation_helper" in audit
    assert "notation helper" in audit
    assert "implemented_with_owner" not in audit


def test_generated_markdown_explains_documentation_only_cpp_owner_absence() -> None:
    entries = [
        {
            "eqid": "documentation_helper",
            "section": "Runtime",
            "subsection": "",
            "subsubsection": "",
            "label": "",
            "status": "Documentation-only",
            "description": "notation helper",
            "tex_file": "docs/latex/equations.tex",
            "tex_line": 10,
            "latex": "a=b",
            "cpp_refs": [],
        },
    ]

    markdown = sync_equation_registry.render_markdown(entries)

    assert "Documentation-only: no direct native owner expected" in markdown
    assert "No `EqID` owner comment has been attached yet" not in markdown


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
    assert "--docs-only-audit" in result.stdout


def test_native_equation_family_coverage_matrix_is_documented_and_complete() -> None:
    docs = (REPO_ROOT / "docs" / "pages" / "equation_traceability.rst").read_text(encoding="utf-8")
    native_contracts = (REPO_ROOT / "tests" / "native" / "test_runtime_contracts.py").read_text(encoding="utf-8")

    assert "Equation-Family Coverage Matrix" in docs
    for family, tests in EQUATION_FAMILY_COVERAGE.items():
        assert family in docs
        for test_name in tests:
            assert test_name in native_contracts



