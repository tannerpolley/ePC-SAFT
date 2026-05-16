from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[3]
IMPORT_BOUNDARY_WATCHLIST = {
    "epcsaft.electrolyte_bubble",
    "epcsaft.epcsaft",
    "epcsaft.eos",
    "epcsaft.equilibrium",
    "epcsaft.properties",
    "epcsaft.reactive_electrolyte",
    "epcsaft.reactive_regression",
    "epcsaft.reactive_speciation",
    "epcsaft.regression",
    "epcsaft.runtime",
}
ANALYSIS_ROOTS = {
    "2012_held": REPO_ROOT / "analyses" / "paper_validation" / "native" / "2012_held",
    "2014_held": REPO_ROOT / "analyses" / "paper_validation" / "native" / "2014_held",
    "2015_baygi": REPO_ROOT / "analyses" / "paper_validation" / "application" / "2015_baygi",
    "2019_bulow": REPO_ROOT / "analyses" / "paper_validation" / "native" / "2019_bulow",
    "2020_bulow": REPO_ROOT / "analyses" / "paper_validation" / "native" / "2020_bulow",
    "2025_figiel": REPO_ROOT / "analyses" / "paper_validation" / "native" / "2025_figiel",
    "2026_khudaida": REPO_ROOT / "analyses" / "paper_validation" / "application" / "2026_khudaida",
    "dielectric_fits": REPO_ROOT / "analyses" / "data_validation" / "dielectric_fits",
    "miac_fits": REPO_ROOT / "analyses" / "data_validation" / "miac_fits",
    "osmotic_validation": REPO_ROOT / "analyses" / "data_validation" / "osmotic_validation",
    "package_plot_smokes": REPO_ROOT / "analyses" / "package_validation" / "package_plot_smokes",
}
MIGRATED_ANALYSIS_IDS = set(ANALYSIS_ROOTS) - {"2025_figiel"}
CATEGORY_ROOTS = {
    REPO_ROOT / "analyses" / "paper_validation",
    REPO_ROOT / "analyses" / "paper_validation" / "native",
    REPO_ROOT / "analyses" / "paper_validation" / "application",
    REPO_ROOT / "analyses" / "data_validation",
    REPO_ROOT / "analyses" / "package_validation",
}
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


def _probe_epcsaft_import_modules(source: str) -> set[str]:
    probe = f"""
import json
import sys

{source}

watchlist = {sorted(IMPORT_BOUNDARY_WATCHLIST)!r}
print(json.dumps(sorted(name for name in watchlist if name in sys.modules)))
"""
    env = os.environ.copy()
    src_root = str(REPO_ROOT / "src")
    env["PYTHONPATH"] = src_root + os.pathsep + env["PYTHONPATH"] if env.get("PYTHONPATH") else src_root
    result = subprocess.run(
        [sys.executable, "-c", probe],
        cwd=REPO_ROOT,
        text=True,
        capture_output=True,
        check=True,
        env=env,
    )
    return set(json.loads(result.stdout))


def _tracked_files(*paths: str) -> list[str]:
    result = subprocess.run(
        ["git", "ls-files", *paths],
        cwd=REPO_ROOT,
        text=True,
        capture_output=True,
        check=True,
    )
    return [line.strip() for line in result.stdout.splitlines() if line.strip()]


def test_strict_solver_derivative_text_gate_passes() -> None:
    subprocess.run(
        [sys.executable, "scripts/dev/check_text_gates.py"],
        cwd=REPO_ROOT,
        text=True,
        check=True,
    )


def test_removed_numerics_stack_is_not_a_package_dev_test_or_analysis_runtime_dependency() -> None:
    removed_dependency_name = "sci" + "py"
    removed_python_ipopt_wrapper = "cy" + "ipopt"
    pyproject_text = (REPO_ROOT / "pyproject.toml").read_text(encoding="utf-8").lower()
    lock_text = (REPO_ROOT / "uv.lock").read_text(encoding="utf-8").lower()
    assert removed_dependency_name not in pyproject_text
    assert removed_python_ipopt_wrapper not in pyproject_text
    assert removed_python_ipopt_wrapper not in lock_text

    tracked = _tracked_files("src", "tests", "scripts", "analyses", "docs", "README.md", "CHANGELOG.md")
    import_offenders: list[str] = []
    import_snippets = (f"import {removed_dependency_name}", f"from {removed_dependency_name}")
    for relpath in tracked:
        if relpath.replace("\\", "/").startswith("docs/papers/"):
            continue
        if not relpath.endswith(".py"):
            if Path(relpath).suffix.lower() not in {".md", ".rst", ".toml", ".yaml", ".yml", ".txt", ".ps1"}:
                continue
        path = REPO_ROOT / relpath
        if not path.exists():
            continue
        text = path.read_text(encoding="utf-8", errors="ignore").lower()
        assert removed_python_ipopt_wrapper not in text, relpath
        if not relpath.endswith(".py"):
            continue
        path = REPO_ROOT / relpath
        if not path.exists():
            continue
        if any(snippet in text for snippet in import_snippets):
            import_offenders.append(relpath)
    assert import_offenders == []
    removed_fit_script = "/".join(
        (
            "analyses",
            "paper_validation",
            "application",
            "2026_rezaee",
            "scripts",
            "rezaee_reactive_" + "equilibrium_fit.py",
        )
    )
    assert not (REPO_ROOT / removed_fit_script).exists()
    assert not (REPO_ROOT / "src" / "epcsaft" / ("_optional" + "_backends")).exists()
    removed_ipopt_helper = REPO_ROOT / "scripts" / "dev" / ("setup_windows_" + removed_python_ipopt_wrapper + "_uv.ps1")
    assert not removed_ipopt_helper.exists()


def test_package_import_is_lazy_across_equilibrium_and_regression_extensions() -> None:
    loaded = _probe_epcsaft_import_modules("import epcsaft")
    assert loaded == set()


def test_eos_and_property_imports_do_not_load_equilibrium_or_regression_extensions() -> None:
    loaded = _probe_epcsaft_import_modules(
        "from epcsaft import ePCSAFTMixture\n"
        "import epcsaft.eos\n"
        "import epcsaft.properties\n"
        "_ = ePCSAFTMixture"
    )
    assert {"epcsaft.epcsaft", "epcsaft.eos", "epcsaft.properties"} <= loaded
    assert loaded.isdisjoint(
        {
            "epcsaft.electrolyte_bubble",
            "epcsaft.equilibrium",
            "epcsaft.reactive_electrolyte",
            "epcsaft.reactive_regression",
            "epcsaft.reactive_speciation",
            "epcsaft.regression",
        }
    )


def test_top_level_public_exports_load_only_the_requested_extension() -> None:
    loaded = _probe_epcsaft_import_modules("import epcsaft\n_ = epcsaft.EquilibriumOptions")
    assert "epcsaft.equilibrium" in loaded
    assert "epcsaft.regression" not in loaded
    assert "epcsaft.reactive_regression" not in loaded


def test_reference_data_root_is_canonical() -> None:
    assert (REPO_ROOT / "data" / "reference" / "epcsaft_parameters").is_dir()
    assert not (REPO_ROOT / "data" / "epcsaft_parameters").exists()


def test_analysis_category_roots_exist() -> None:
    for root in sorted(CATEGORY_ROOTS):
        assert root.is_dir(), root


def test_migrated_analyses_have_local_contract_files() -> None:
    for analysis_id, root in sorted(ANALYSIS_ROOTS.items()):
        assert (root / "README.md").is_file(), analysis_id
        assert (root / "analysis.yaml").is_file(), analysis_id
        scripts_root = root / "scripts"
        figure_scripts = sorted((root / "figures").glob("*/scripts"))
        assert scripts_root.is_dir() or figure_scripts, analysis_id


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
        or ("/figures/" in path.replace("\\", "/") and "/output/runs/" in path.replace("\\", "/"))
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
        path = ANALYSIS_ROOTS[analysis_id] / "analysis.yaml"
        text = path.read_text(encoding="utf-8")
        assert "figures: figures/<figure_id>/output" in text, path
        assert "runs: figures/<figure_id>/output/runs" in text, path


def test_migrated_analyses_use_complete_figure_owned_roots() -> None:
    for analysis_id in sorted(MIGRATED_ANALYSIS_IDS):
        figures_root = ANALYSIS_ROOTS[analysis_id] / "figures"
        assert figures_root.is_dir(), analysis_id
        figure_roots = sorted(path for path in figures_root.iterdir() if path.is_dir())
        assert figure_roots, analysis_id
        for figure_root in figure_roots:
            assert (figure_root / "scripts").is_dir(), figure_root
            optional_roots = (figure_root / "input", figure_root / "output")
            assert any(path.is_dir() for path in optional_roots) or any(
                (figure_root / "scripts").glob("*.py")
            ), figure_root


def test_selected_figure_scripts_do_not_read_canonical_data_root_directly() -> None:
    analysis_ids = ("2012_held", "2014_held", "2019_bulow", "2020_bulow")
    forbidden_snippets = (
        'REPO_ROOT / "data"',
        "REPO_ROOT / 'data'",
        'common.REPO_ROOT / "data"',
        "common.REPO_ROOT / 'data'",
        "data/reference",
    )
    for analysis_id in analysis_ids:
        figures_root = ANALYSIS_ROOTS[analysis_id] / "figures"
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
    root = ANALYSIS_ROOTS["2025_figiel"]
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
