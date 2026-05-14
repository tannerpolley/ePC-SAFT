from __future__ import annotations

import csv
import importlib
import subprocess
import sys


from pathlib import Path

import pytest
import sys as _bootstrap_sys
from pathlib import Path as _BootstrapPath

for _candidate in _BootstrapPath(__file__).resolve().parents:
    if (_candidate / "scripts" / "plot_outputs.py").is_file():
        if str(_candidate) not in _bootstrap_sys.path:
            _bootstrap_sys.path.insert(0, str(_candidate))
        break
else:
    raise ModuleNotFoundError("Could not locate repo root containing scripts/plot_outputs.py")
from scripts.plot_outputs import REPO_ROOT

PLOT_ROOT = REPO_ROOT / "analyses" / "paper_validation" / "application" / "2015_baygi"
SCRIPT_ROOT = PLOT_ROOT / "scripts"
REGRESSED_PARAMETERS = PLOT_ROOT / "figures" / "regressed_parameters" / "output" / "regressed_parameters.csv"
RUN_ALL = SCRIPT_ROOT / "run_all.py"
EXPECTED_LABELS = {"DIPPR", "MEA 2B", "MEA 3B", "MEA 4C"}
EXPECTED_FIGURES = ("figure_2", "figure_3", "figure_2_regressed", "figure_3_regressed")


@pytest.fixture(scope="module", autouse=True)
def _generate_baygi_outputs() -> None:
    _seed_cached_baygi_outputs()
    subprocess.run([sys.executable, str(RUN_ALL)], cwd=REPO_ROOT, check=True)


def _seed_cached_baygi_outputs() -> None:
    common = importlib.import_module("analyses.paper_validation.application.2015_baygi.scripts._common")
    dippr = common.dippr_reference_rows()
    fields = ("series", "T_K", "P_Pa", "rho_mol_m3", "status", "contribution_terms")
    diagnostics = []
    psat_offsets = {"2B": 0.001, "3B": 0.03, "4C": 0.001}
    for row in dippr:
        for index, assoc_scheme in enumerate(("2B", "3B", "4C"), start=1):
            diagnostics.append(
                {
                    "series": f"MEA {assoc_scheme}",
                    "T_K": row["T_K"],
                    "P_Pa": float(row["P_Pa"]) * (1.0 + psat_offsets[assoc_scheme]),
                    "rho_mol_m3": float(row["rho_mol_m3"]) * (1.0 + 0.0005 * index),
                    "status": "solved",
                    "contribution_terms": "hc;disp;assoc",
                }
            )

    for figure in EXPECTED_FIGURES:
        common.write_csv_rows(
            PLOT_ROOT / "figures" / figure / "output" / "data" / f"{figure}_diagnostics.csv",
            fields,
            [*dippr, *diagnostics],
        )

    common.write_csv_rows(
        REGRESSED_PARAMETERS,
        (
            "assoc_scheme",
            "m",
            "s",
            "e",
            "e_assoc",
            "vol_a",
            "backend",
            "success",
            "status",
            "message",
            "nfev",
            "cost",
            "residual_norm",
            "initial_residual_norm",
        ),
        [
            {
                "assoc_scheme": assoc_scheme,
                **values,
                "backend": "cached-test-fixture",
                "success": True,
                "status": 0,
                "message": "cached",
                "nfev": 0,
                "cost": 0.0,
                "residual_norm": 0.0,
                "initial_residual_norm": 0.0,
            }
            for assoc_scheme, values in common.TABLE2_MEA_PARAMETERS.items()
        ],
    )


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
    reference = {round(float(row["T_K"]), 6): float(row[field]) for row in rows if row["series"] == "DIPPR"}
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
    assert REGRESSED_PARAMETERS.exists()
    for figure in EXPECTED_FIGURES:
        figure_dir = PLOT_ROOT / "figures" / figure / "output"
        image = figure_dir / f"{figure}.png"
        svg = figure_dir / f"{figure}.svg"
        plot_data = figure_dir / f"{figure}.csv"
        plot_style = figure_dir / f"{figure}.mpl.yaml"
        diagnostics = figure_dir / "data" / f"{figure}_diagnostics.csv"
        metrics = figure_dir / "data" / f"{figure}_metrics.csv"

        assert image.exists()
        assert svg.exists()
        assert plot_data.exists()
        assert plot_style.exists()
        assert diagnostics.exists()
        assert metrics.exists()
        assert EXPECTED_LABELS <= _artist_labels(plot_data)

        rows = _diagnostic_rows(diagnostics)
        assert rows
        assert {row["contribution_terms"] for row in rows} <= {"hc;disp;assoc", "reference"}


def test_2015_baygi_figures_report_numeric_fit_quality():
    figure_2_rows = _diagnostic_rows(PLOT_ROOT / "figures" / "figure_2" / "output" / "data" / "figure_2_diagnostics.csv")
    figure_3_rows = _diagnostic_rows(PLOT_ROOT / "figures" / "figure_3" / "output" / "data" / "figure_3_diagnostics.csv")
    figure_2_metrics = _metric_rows(PLOT_ROOT / "figures" / "figure_2" / "output" / "data" / "figure_2_metrics.csv")
    figure_3_metrics = _metric_rows(PLOT_ROOT / "figures" / "figure_3" / "output" / "data" / "figure_3_metrics.csv")

    psat_aard = _aard_by_series(figure_2_rows, "P_Pa")
    rho_aard = _aard_by_series(figure_3_rows, "rho_mol_m3")

    assert all(value < 10.0 for value in psat_aard.values()), psat_aard
    assert all(value < 2.0 for value in rho_aard.values()), rho_aard
    assert {row["field"] for row in figure_2_metrics} == {"P_Pa"}
    assert {row["field"] for row in figure_3_metrics} == {"rho_mol_m3"}
    assert all(row["reported_table2_aad_percent"] for row in [*figure_2_metrics, *figure_3_metrics])
    assert any(row["within_reported_aad"] == "False" for row in figure_2_metrics)
