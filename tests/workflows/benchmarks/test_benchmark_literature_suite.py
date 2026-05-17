from __future__ import annotations

import json
import subprocess
import sys
from dataclasses import replace


def test_literature_benchmark_module_exposes_issue_119_case_order():
    from scripts.benchmarks.helpers.literature import LITERATURE_CASES

    assert tuple(LITERATURE_CASES) == (
        "gross_sadowski_pure_nonassociating",
        "gross_sadowski_associating_systems",
        "baygi_mea_association_and_mea_water_binary",
        "cameretti_held_aqueous_electrolyte_density_miac",
        "held_mixed_solvent_density_osmotic_miac",
        "bulow_ascani_dielectric_born",
        "figiel_2025_ssm_ds_born",
        "ascani_2022_distributed_ion_lle",
        "ascani_2023_reactive_phase_equilibrium",
        "khudaida_2026_salting_out_lle",
        "rezaee_lithium_extraction_inputs",
        "mea_true_species_pressure_speciation",
    )


def test_literature_benchmark_payload_tracks_executable_and_blocked_issue_119_cases():
    from scripts.benchmarks.helpers.literature import BLOCKED, EXECUTABLE, run_literature_benchmarks

    payload = run_literature_benchmarks()

    assert payload["issue"] == 119
    assert payload["title"] == "Executable literature benchmark and downstream gate registry"
    assert payload["selected_cases"][0] == "gross_sadowski_pure_nonassociating"
    assert payload["status_counts"] == {EXECUTABLE: 7, BLOCKED: 5}
    assert len(payload["executable_cases"]) == 7
    assert len(payload["blocked_cases"]) == 5

    by_case = {row["id"]: row for row in payload["cases"]}
    assert by_case["gross_sadowski_pure_nonassociating"]["status"] == EXECUTABLE
    assert by_case["gross_sadowski_pure_nonassociating"]["expected"] is not None
    assert by_case["gross_sadowski_pure_nonassociating"]["tolerances"] is not None
    assert by_case["figiel_2025_ssm_ds_born"]["status"] == EXECUTABLE
    assert by_case["figiel_2025_ssm_ds_born"]["expected"]["miac_probe"] == 0.7732309439080085
    assert by_case["khudaida_2026_salting_out_lle"]["status"] == EXECUTABLE
    assert by_case["baygi_mea_association_and_mea_water_binary"]["status"] == BLOCKED
    assert by_case["baygi_mea_association_and_mea_water_binary"]["expected"] is None
    assert by_case["baygi_mea_association_and_mea_water_binary"]["tolerances"] is None
    assert by_case["ascani_2023_reactive_phase_equilibrium"]["status"] == BLOCKED
    assert by_case["ascani_2023_reactive_phase_equilibrium"]["blocked_by_issue"] == 119
    assert by_case["rezaee_lithium_extraction_inputs"]["status"] == EXECUTABLE
    assert by_case["rezaee_lithium_extraction_inputs"]["expected"]["direct_published_constant_closure_supported"] is False
    assert by_case["rezaee_lithium_extraction_inputs"]["validation_paths"]


def test_literature_benchmark_payload_executes_only_executable_cases_with_injected_runner():
    from scripts.benchmarks.helpers.literature import BLOCKED, EXECUTABLE, run_literature_benchmarks

    calls: list[str] = []

    def fake_runner(command: str) -> dict[str, object]:
        calls.append(command)
        return {
            "returncode": 0,
            "duration_seconds": 0.01,
            "stdout": f"ran {command}",
            "stderr": "",
        }

    payload = run_literature_benchmarks(execute_commands=True, command_runner=fake_runner)

    assert payload["run_mode"] == "execute_executable_cases"
    assert payload["execution_summary"] == {
        "attempted": 7,
        "passed": 7,
        "failed": 0,
        "blocked": 5,
    }
    assert len(calls) == 7

    by_case = {row["id"]: row for row in payload["cases"]}
    assert by_case["gross_sadowski_pure_nonassociating"]["status"] == EXECUTABLE
    assert by_case["gross_sadowski_pure_nonassociating"]["execution"]["attempted"] is True
    assert by_case["gross_sadowski_pure_nonassociating"]["execution"]["passed"] is True
    assert by_case["baygi_mea_association_and_mea_water_binary"]["status"] == BLOCKED
    assert by_case["baygi_mea_association_and_mea_water_binary"]["execution"]["attempted"] is False
    assert by_case["baygi_mea_association_and_mea_water_binary"]["execution"]["passed"] is None


def test_literature_benchmark_payload_rejects_executable_case_without_tolerance_contract(monkeypatch):
    import pytest

    from scripts.benchmarks.helpers import literature as bench

    broken = replace(bench.LITERATURE_CASES["gross_sadowski_pure_nonassociating"], tolerances=None)
    monkeypatch.setitem(bench.LITERATURE_CASES, "gross_sadowski_pure_nonassociating", broken)

    with pytest.raises(ValueError, match="missing expected values or tolerances"):
        bench.run_literature_benchmarks()


def test_literature_benchmark_script_runs_and_writes_json(tmp_path):
    from scripts.benchmarks.helpers.literature import EXECUTABLE

    output_path = tmp_path / "literature_suite.json"

    result = subprocess.run(
        [
            sys.executable,
            "scripts/benchmarks/benchmark_literature_suite.py",
            "--case",
            "figiel_2025_ssm_ds_born",
            "--json",
            str(output_path),
        ],
        capture_output=True,
        text=True,
        check=False,
    )

    assert result.returncode == 0, result.stdout + result.stderr
    assert "figiel_2025_ssm_ds_born" in result.stdout
    assert "status" in result.stdout
    assert output_path.exists()

    payload = json.loads(output_path.read_text(encoding="utf-8"))
    assert payload["issue"] == 119
    assert payload["selected_cases"] == ["figiel_2025_ssm_ds_born"]
    assert payload["run_mode"] == "execute_executable_cases"
    assert payload["execution_summary"] == {
        "attempted": 1,
        "passed": 1,
        "failed": 0,
        "blocked": 0,
    }
    assert payload["status_counts"] == {EXECUTABLE: 1}
    assert payload["cases"][0]["status"] == EXECUTABLE
    assert payload["cases"][0]["expected"] is not None
    assert payload["cases"][0]["tolerances"] is not None
    assert payload["cases"][0]["execution"]["attempted"] is True
    assert payload["cases"][0]["execution"]["passed"] is True
