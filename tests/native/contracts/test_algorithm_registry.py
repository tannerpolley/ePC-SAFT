from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest

from scripts.docs import sync_algorithm_registry

REPO_ROOT = Path(__file__).resolve().parents[3]
TEX_PATH = REPO_ROOT / "docs" / "latex" / "algorithms.tex"


def test_algorithm_registry_outputs_are_synced() -> None:
    assert TEX_PATH.exists(), "docs/latex/algorithms.tex must be present as a tracked repo file"

    result = subprocess.run(
        [sys.executable, "scripts/docs/sync_algorithm_registry.py", "--check"],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
        check=False,
    )

    assert result.returncode == 0, result.stdout + result.stderr


def test_algorithm_registry_strict_traceability_passes_current_registry() -> None:
    assert TEX_PATH.exists(), "docs/latex/algorithms.tex must be present as a tracked repo file"

    result = subprocess.run(
        [sys.executable, "scripts/docs/sync_algorithm_registry.py", "--check", "--strict-traceability"],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
        check=False,
    )

    assert result.returncode == 0, result.stdout + result.stderr


def test_sync_algorithm_registry_missing_tracked_latex_reports_actionable_message(tmp_path, capsys) -> None:
    missing_tex = tmp_path / "docs" / "latex" / "algorithms.tex"

    with pytest.raises(SystemExit) as excinfo:
        sync_algorithm_registry.require_algorithms_tex(missing_tex)

    assert excinfo.value.code == 1
    combined = capsys.readouterr().err
    assert "docs/latex/algorithms.tex" in combined
    assert "tracked algorithm registry source" in combined
    assert "Traceback" not in combined


def test_parse_algorithms_rejects_duplicate_algids(tmp_path) -> None:
    tex = tmp_path / "algorithms.tex"
    tex.write_text(
        "\n".join(
            [
                r"\section{Algorithms}",
                "% AlgID: duplicate_route",
                "% Status: Implemented",
                "First route.",
                "% AlgID: duplicate_route",
                "% Status: Implemented",
                "Second route.",
            ]
        ),
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="duplicate_route"):
        sync_algorithm_registry.parse_algorithms(tex)


def test_validate_links_fails_unknown_source_algid(tmp_path) -> None:
    source = tmp_path / "src"
    source.mkdir()
    owner = source / "owner.py"
    owner.write_text("# " + "AlgID" + ": unknown_route\n" "def route():\n" "    return None\n", encoding="utf-8")

    refs = sync_algorithm_registry.parse_code_refs([source])

    with pytest.raises(ValueError, match="unknown_route"):
        sync_algorithm_registry.validate_links([{"algid": "known_route"}], refs)


def test_strict_traceability_allows_documentation_only_and_planned_entries() -> None:
    entries = [
        {"algid": "docs_only", "status": "Documentation-only", "owner_refs": []},
        {"algid": "planned_route", "status": "planned", "owner_refs": []},
        {"algid": "implemented_route", "status": "Implemented", "owner_refs": [{"file": "x.py"}]},
    ]

    sync_algorithm_registry.enforce_traceability(entries)


def test_strict_traceability_fails_missing_implemented_owner() -> None:
    entries = [{"algid": "implemented_without_owner", "section": "Routes", "status": "Implemented", "owner_refs": []}]

    with pytest.raises(SystemExit) as excinfo:
        sync_algorithm_registry.enforce_traceability(entries)

    assert "implemented_without_owner" in str(excinfo.value)


def test_python_and_cpp_owner_comments_attach_to_next_code_line(tmp_path) -> None:
    source = tmp_path / "src"
    source.mkdir()
    py_owner = source / "owner.py"
    py_owner.write_text(
        "# " + "AlgID" + ": python_route\n\n" "def python_route():\n" "    return None\n",
        encoding="utf-8",
    )
    cpp_owner = source / "owner.cpp"
    cpp_owner.write_text(
        "// " + "AlgID" + ": cpp_route\n\n" "int cpp_route() {\n" "    return 0;\n" "}\n",
        encoding="utf-8",
    )

    refs = sync_algorithm_registry.parse_code_refs([source])

    assert refs["python_route"][0]["line"] == 3
    assert refs["python_route"][0]["context"] == "def python_route():"
    assert refs["cpp_route"][0]["line"] == 3
    assert refs["cpp_route"][0]["context"] == "int cpp_route() {"


def test_ceres_and_ipopt_algorithm_ids_are_present_and_distinct() -> None:
    entries = sync_algorithm_registry.parse_algorithms(TEX_PATH)
    algids = {entry["algid"] for entry in entries}

    assert "ipopt_tnlp_adapter" in algids
    assert "native_ceres_regression_adapter" in algids
    assert "neutral_tp_flash_ipopt" in algids
    assert "pure_neutral_ceres_regression" in algids
    for entry in entries:
        if entry["algid"].endswith("_ipopt"):
            assert "Ceres" not in str(entry.get("backend", "")) + str(entry.get("solver_role", ""))
        if entry["algid"].endswith("_ceres_regression") or entry["algid"] == "native_ceres_regression_adapter":
            assert "Ipopt" not in str(entry.get("backend", "")) + str(entry.get("solver_role", ""))


def test_generated_markdown_names_public_api_backend_and_dependencies() -> None:
    entries = sync_algorithm_registry.parse_algorithms(TEX_PATH)
    math_map = sync_algorithm_registry.build_algorithm_math_map(entries)
    markdown = sync_algorithm_registry.render_markdown(entries, math_map)

    assert 'mix.equilibrium(kind="tp_flash", ...)' in markdown
    assert "mix.solve_equilibrium(TPFlash(...))" in markdown
    assert "Native C++ Ipopt equilibrium NLP" in markdown
    assert "Dependency: Ipopt" in markdown
    assert "Native C++ Ceres regression" in markdown
    assert "Dependency: Ceres" in markdown
    assert "not a production Ceres optimizer" in markdown
    assert "**LaTeX source**" in markdown
    assert "```tex" in markdown
    assert "**Rendered formulae**" in markdown
    assert "$$\nl_{ij}\n$$" in markdown
    assert "$$\nk^{hb}_{ij}\n$$" in markdown
    assert "`l_{ij}`" in markdown
    assert "`k^{hb}_{ij}`" in markdown
    assert "`NlpProblem`" in markdown
    assert r"\end{document}" not in markdown
