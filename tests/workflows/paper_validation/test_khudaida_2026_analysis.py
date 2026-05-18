from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[3]
KHUDAIDA = REPO_ROOT / "analyses" / "paper_validation" / "application" / "2026_khudaida"
INPUT = KHUDAIDA / "data" / "input"
PARAMS = REPO_ROOT / "data" / "reference" / "epcsaft_parameters" / "2026_Khudaida"


def _rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def test_khudaida_source_tables_cover_parameters_and_born_options() -> None:
    pure_rows = {row["component"]: row for row in _rows(INPUT / "table_5_pure_component_parameters.csv")}
    dielectric_rows = {row["component"]: row for row in _rows(INPUT / "table_6_relative_dielectric_constants.csv")}
    kij_rows = _rows(INPUT / "table_7_epcsaft_kij.csv")
    source_options = json.loads((INPUT / "born_ssm_ds_options.json").read_text(encoding="utf-8"))
    runtime_options = json.loads((PARAMS / "user_options.json").read_text(encoding="utf-8"))

    assert set(pure_rows) == {"H2O", "Ethanol", "Butanol", "Na+", "Cl-"}
    assert float(pure_rows["Na+"]["d_born_A"]) == pytest.approx(3.445)
    assert float(pure_rows["Cl-"]["d_born_A"]) == pytest.approx(4.100)
    assert float(pure_rows["Ethanol"]["f_solv"]) == pytest.approx(1.6)
    assert set(dielectric_rows) == {"water", "ethanol", "isobutanol", "NaCl ions"}
    assert len(kij_rows) == 10
    assert source_options == runtime_options
    born_model = runtime_options["elec_model"]["born_model"]
    assert born_model["d_Born_mode"] == 3
    assert born_model["solvation_shell_model"] is True
    assert born_model["dielectric_saturation"] is True
    assert born_model["mu_born_model"]["differential_mode"] == "analytical"
    assert born_model["mu_born_model"]["comp_dep_delta_d"] is True


def test_khudaida_source_tables_cover_paper_and_si_lle_rows() -> None:
    tielines = _rows(INPUT / "table_3_4_experimental_tielines.csv")
    metrics = _rows(INPUT / "table_s1_s2_distribution_separation.csv")
    enrtl_bips = _rows(INPUT / "table_s3_enrtl_bips.csv")
    enrtl_aad = _rows(INPUT / "table_s4_enrtl_aad.csv")

    assert len(tielines) == 39
    assert len(metrics) == 39
    assert len(enrtl_bips) == 6
    assert len(enrtl_aad) == 12
    assert {row["source_table"] for row in tielines} == {"Table 3", "Table 4"}
    assert {row["source_table"] for row in metrics} == {"Table S1", "Table S2"}
    for row in tielines:
        organic = [float(row[key]) for key in ("x_water_org", "x_ethanol_org", "x_isobutanol_org", "x_nacl_org")]
        aqueous = [float(row[key]) for key in ("x_water_aq", "x_ethanol_aq", "x_isobutanol_aq", "x_nacl_aq")]
        # Retain the paper-table values as source data; one Table 4 organic row sums to 0.99 in the local source text.
        assert sum(organic) == pytest.approx(1.0, abs=1.1e-2)
        assert sum(aqueous) == pytest.approx(1.0, abs=1.1e-2)


def test_khudaida_supporting_subfigure_workflows_are_declared() -> None:
    manifest = _rows(INPUT / "figure_manifest.csv")
    panels = {(row["figure_id"], row["panel"]) for row in manifest}

    assert {("figure_s2", panel) for panel in ("a", "b", "c")} <= panels
    assert {("figure_s3", panel) for panel in ("a", "b", "c")} <= panels
    for figure_id in ("figure_s2", "figure_s3"):
        scripts = KHUDAIDA / "figures" / figure_id / "scripts"
        assert (scripts / "generate_data.py").is_file()
        assert (scripts / f"plot_{figure_id}.py").is_file()
