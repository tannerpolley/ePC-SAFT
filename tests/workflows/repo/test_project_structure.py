from __future__ import annotations

import subprocess
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[3]
ANALYSIS_IDS = {
    "2012_held",
    "2014_held",
    "2015_baygi",
    "2019_bulow",
    "2020_bulow",
    "2025_figiel",
    "2026_khudaida",
    "dielectric_fits",
    "miac_fits",
    "osmotic_validation",
    "package_plot_smokes",
}
MIGRATED_ANALYSIS_IDS = ANALYSIS_IDS - {"2025_figiel"}
TEST_SUBGROUP_ROOTS = {
    "tests/api/package",
    "tests/api/parameters",
    "tests/api/reactive",
    "tests/api/regression",
    "tests/api/runtime",
    "tests/equilibrium/core",
    "tests/equilibrium/electrolyte",
    "tests/equilibrium/reactive",
    "tests/native/ceres",
    "tests/native/contracts",
    "tests/native/cppad",
    "tests/native/equilibrium",
    "tests/native/runtime",
    "tests/regression/core",
    "tests/regression/electrolyte",
    "tests/regression/literature",
    "tests/workflows/benchmarks",
    "tests/workflows/build",
    "tests/workflows/repo",
}
REPLACED_FLAT_TEST_FILES = {
    "tests/api/test_runtime.py",
    "tests/api/test_regression_api.py",
    "tests/api/test_reactive_speciation.py",
    "tests/api/test_reactive_regression.py",
    "tests/api/test_reactive_electrolyte_bubble.py",
    "tests/equilibrium/test_electrolyte_lle.py",
    "tests/native/test_runtime_contracts.py",
    "tests/native/test_chemical_equilibrium_native.py",
}


def _tracked_files(*paths: str) -> list[str]:
    result = subprocess.run(
        ["git", "ls-files", *paths],
        cwd=REPO_ROOT,
        text=True,
        capture_output=True,
        check=True,
    )
    return [line.strip() for line in result.stdout.splitlines() if line.strip()]


def test_reference_data_root_is_canonical() -> None:
    assert (REPO_ROOT / "data" / "reference" / "epcsaft_parameters").is_dir()
    assert not (REPO_ROOT / "data" / "epcsaft_parameters").exists()


def test_migrated_analyses_have_local_contract_files() -> None:
    for analysis_id in sorted(ANALYSIS_IDS):
        root = REPO_ROOT / "analyses" / analysis_id
        assert (root / "README.md").is_file(), analysis_id
        assert (root / "analysis.yaml").is_file(), analysis_id
        assert (root / "scripts").is_dir(), analysis_id


def test_old_gallery_and_script_roots_are_not_tracked() -> None:
    assert _tracked_files("docs/plots") == []
    assert _tracked_files("scripts/paper_validation") == []
    assert _tracked_files("scripts/fits") == []
    assert _tracked_files("tests/plots") == []


def test_reorganized_test_subgroup_roots_exist() -> None:
    for relpath in sorted(TEST_SUBGROUP_ROOTS):
        assert (REPO_ROOT / relpath).is_dir(), relpath
        assert (REPO_ROOT / relpath / "__init__.py").is_file(), relpath


def test_replaced_flat_test_modules_are_absent_from_the_working_tree() -> None:
    for relpath in sorted(REPLACED_FLAT_TEST_FILES):
        assert not (REPO_ROOT / relpath).exists(), relpath


def test_generated_output_roots_are_not_tracked_in_analyses() -> None:
    tracked = _tracked_files("analyses")
    stale = [
        path
        for path in tracked
        if "/out/" in path.replace("\\", "/")
        or "/results/runs/" in path.replace("\\", "/")
        or "/results/final/" in path.replace("\\", "/")
        or "/figures/" in path.replace("\\", "/") and "/output/runs/" in path.replace("\\", "/")
    ]
    assert stale == []


def test_analysis_metadata_does_not_reference_removed_final_results_layout() -> None:
    tracked = _tracked_files("analyses")
    metadata_files = [REPO_ROOT / path for path in tracked if path.endswith("/analysis.yaml")]
    assert metadata_files
    for path in metadata_files:
        text = path.read_text(encoding="utf-8")
        assert "results/final" not in text, path


def test_migrated_analysis_metadata_uses_figure_owned_outputs() -> None:
    for analysis_id in sorted(MIGRATED_ANALYSIS_IDS):
        path = REPO_ROOT / "analyses" / analysis_id / "analysis.yaml"
        text = path.read_text(encoding="utf-8")
        assert "figures: figures/<figure_id>/output" in text, path
        assert "runs: figures/<figure_id>/output/runs" in text, path


def test_migrated_analyses_use_complete_figure_owned_roots() -> None:
    for analysis_id in sorted(MIGRATED_ANALYSIS_IDS):
        figures_root = REPO_ROOT / "analyses" / analysis_id / "figures"
        assert figures_root.is_dir(), analysis_id
        figure_roots = sorted(path for path in figures_root.iterdir() if path.is_dir())
        assert figure_roots, analysis_id
        for figure_root in figure_roots:
            assert (figure_root / "input").is_dir(), figure_root
            assert (figure_root / "output").is_dir(), figure_root
            assert (figure_root / "scripts").is_dir(), figure_root


def test_selected_figure_scripts_do_not_read_canonical_data_root_directly() -> None:
    analysis_ids = ("2012_held", "2014_held", "2019_bulow", "2020_bulow")
    forbidden_snippets = (
        'REPO_ROOT / "data"',
        "REPO_ROOT / 'data'",
        'common.REPO_ROOT / "data"',
        "common.REPO_ROOT / 'data'",
        'data/reference',
    )
    for analysis_id in analysis_ids:
        figures_root = REPO_ROOT / "analyses" / analysis_id / "figures"
        for path in sorted(figures_root.rglob("*.py")):
            if path.name == "generate_data.py":
                continue
            text = path.read_text(encoding="utf-8")
            for snippet in forbidden_snippets:
                assert snippet not in text, f"{path} still references canonical data root via {snippet!r}"


def test_analysis_template_uses_figure_owned_outputs() -> None:
    text = (REPO_ROOT / "analyses" / "_template" / "analysis.yaml").read_text(encoding="utf-8")
    assert "figures: figures/<figure_id>/output" in text
    assert "runs: figures/<figure_id>/output/runs" in text


def test_figiel_analysis_is_migrated_to_figure_owned_layout() -> None:
    root = REPO_ROOT / "analyses" / "2025_figiel"
    text = (root / "analysis.yaml").read_text(encoding="utf-8")
    assert "figures: figures/<figure_id>/output" in text
    assert "runs: figures/<figure_id>/output/runs" in text
    assert not (root / "data").exists()
    assert (root / "figures").is_dir()
    for figure_id in ("figure_4", "figure_5", "figure_6", "figure_7", "figure_8", "figure_9"):
        figure_root = root / "figures" / figure_id
        assert (figure_root / "input").is_dir(), figure_id
        assert (figure_root / "output").is_dir(), figure_id
        assert (figure_root / "scripts").is_dir(), figure_id
