from __future__ import annotations

import csv
import subprocess
import sys
from pathlib import Path

import pytest


REPO_ROOT = Path(__file__).resolve().parents[2]
RUN_ALL = REPO_ROOT / "scripts" / "paper_validation" / "2015_Baygi_analysis" / "run_all.py"
PLOT_ROOT = REPO_ROOT / "docs" / "plots" / "paper_validation" / "2015_Baygi"
EXPECTED_LABELS = {"DIPPR", "MEA 2B", "MEA 3B", "MEA 4C"}
EXPECTED_FIGURES = ("figure_2", "figure_3", "figure_2_regressed", "figure_3_regressed")


@pytest.fixture(scope="module", autouse=True)
def _generate_baygi_outputs() -> None:
    subprocess.run([sys.executable, str(RUN_ALL)], cwd=REPO_ROOT, check=True)


def _artist_labels(path: Path) -> set[str]:
    with path.open("r", newline="", encoding="utf-8-sig") as handle:
        return {
            str(row.get("artist_label", "")).strip()
            for row in csv.DictReader(handle)
            if str(row.get("artist_label", "")).strip()
        }


def _diagnostic_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def _metric_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def _aard_by_series(rows: list[dict[str, str]], field: str) -> dict[str, float]:
    reference = {
        round(float(row["T_K"]), 6): float(row[field])
        for row in rows
        if row["series"] == "DIPPR"
    }
    out: dict[str, float] = {}
    for series in ("MEA 2B", "MEA 3B", "MEA 4C"):
        residuals = []
        for row in rows:
            if row["series"] != series or row["status"] != "solved":
                continue
            expected = reference[round(float(row["T_K"]), 6)]
            residuals.append(abs(float(row[field]) / expected - 1.0) * 100.0)
        assert len(residuals) >= 20
        out[series] = sum(residuals) / len(residuals)
    return out


def test_2015_baygi_figure_2_and_3_workflow_outputs_expected_series():
    assert (PLOT_ROOT / "data" / "regressed_parameters.csv").exists()
    for figure in EXPECTED_FIGURES:
        figure_dir = PLOT_ROOT / figure
        image = figure_dir / f"{figure}.png"
        svg = figure_dir / f"{figure}.svg"
        plot_data = figure_dir / "data" / f"{figure}_plot_data.csv"
        diagnostics = figure_dir / "data" / f"{figure}_diagnostics.csv"
        metrics = figure_dir / "data" / f"{figure}_metrics.csv"

        assert image.exists()
        assert svg.exists()
        assert plot_data.exists()
        assert diagnostics.exists()
        assert metrics.exists()
        assert EXPECTED_LABELS <= _artist_labels(plot_data)

        rows = _diagnostic_rows(diagnostics)
        assert rows
        assert {row["contribution_terms"] for row in rows} <= {"hc;disp;assoc", "reference"}


def test_2015_baygi_figures_report_numeric_fit_quality():
    figure_2_rows = _diagnostic_rows(PLOT_ROOT / "figure_2" / "data" / "figure_2_diagnostics.csv")
    figure_3_rows = _diagnostic_rows(PLOT_ROOT / "figure_3" / "data" / "figure_3_diagnostics.csv")
    figure_2_metrics = _metric_rows(PLOT_ROOT / "figure_2" / "data" / "figure_2_metrics.csv")
    figure_3_metrics = _metric_rows(PLOT_ROOT / "figure_3" / "data" / "figure_3_metrics.csv")

    psat_aard = _aard_by_series(figure_2_rows, "P_Pa")
    rho_aard = _aard_by_series(figure_3_rows, "rho_mol_m3")

    assert all(value < 10.0 for value in psat_aard.values()), psat_aard
    assert all(value < 2.0 for value in rho_aard.values()), rho_aard
    assert {row["field"] for row in figure_2_metrics} == {"P_Pa"}
    assert {row["field"] for row in figure_3_metrics} == {"rho_mol_m3"}
    assert all(row["reported_table2_aad_percent"] for row in [*figure_2_metrics, *figure_3_metrics])
    assert any(row["within_reported_aad"] == "False" for row in figure_2_metrics)
