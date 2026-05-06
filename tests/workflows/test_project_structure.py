from __future__ import annotations

import subprocess
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
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
        assert (root / "scripts").is_dir() or analysis_id == "package_plot_smokes", analysis_id


def test_old_gallery_and_script_roots_are_not_tracked() -> None:
    assert _tracked_files("docs/plots") == []
    assert _tracked_files("scripts/paper_validation") == []
    assert _tracked_files("scripts/fits") == []
    assert _tracked_files("tests/plots") == []


def test_generated_output_roots_are_not_tracked_in_analyses() -> None:
    tracked = _tracked_files("analyses")
    stale = [
        path for path in tracked if "/out/" in path.replace("\\", "/") or "/results/runs/" in path.replace("\\", "/")
    ]
    assert stale == []
