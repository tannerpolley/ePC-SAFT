from __future__ import annotations

import json
import subprocess
import sys


def test_literature_benchmark_module_exposes_issue_scope_order():
    from scripts.benchmarks.helpers.literature import LITERATURE_CASES

    assert tuple(LITERATURE_CASES) == (
        "mea_simple_workflow",
        "mdea_epcsaft",
        "figiel_2025_ssm_ds_born",
        "held_2014_revised_epcsaft",
        "non_electrolyte_lle",
        "ascani_2022_electrolyte_lle",
        "ascani_2023_reactive_lle",
        "khudaida_2026_salting_out_lle",
        "hubach_yu_lithium_equilibrium",
    )


def test_literature_benchmark_payload_classifies_supported_and_blocked_cases():
    from scripts.benchmarks.helpers.literature import run_literature_benchmarks

    payload = run_literature_benchmarks()

    assert payload["issue"] == 95
    assert payload["selected_cases"][0] == "mea_simple_workflow"
    assert payload["classification_counts"]["already_supported_with_tests"] == 4
    assert payload["classification_counts"]["blocker_requires_followup"] == 5

    by_case = {row["case"]: row for row in payload["cases"]}
    assert by_case["mea_simple_workflow"]["classification"] == "already_supported_with_tests"
    assert "application-specific" in by_case["mea_simple_workflow"]["notes"]
    assert by_case["figiel_2025_ssm_ds_born"]["classification"] == "already_supported_with_tests"
    assert by_case["khudaida_2026_salting_out_lle"]["classification"] == "already_supported_with_tests"
    assert by_case["ascani_2023_reactive_lle"]["classification"] == "blocker_requires_followup"
    assert by_case["mdea_epcsaft"]["validation_paths"] == []


def test_literature_benchmark_script_runs_and_writes_json(tmp_path):
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
    assert output_path.exists()

    payload = json.loads(output_path.read_text(encoding="utf-8"))
    assert payload["selected_cases"] == ["figiel_2025_ssm_ds_born"]
    assert payload["cases"][0]["classification"] == "already_supported_with_tests"
    assert payload["cases"][0]["coverage_kind"] == "smoke_regression"
