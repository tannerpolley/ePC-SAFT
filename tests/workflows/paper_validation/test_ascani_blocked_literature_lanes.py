from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[3]
ASCANI_2022 = REPO_ROOT / "analyses" / "paper_validation" / "native" / "2022_ascani"
ASCANI_2023 = REPO_ROOT / "analyses" / "paper_validation" / "native" / "2023_ascani"


def _run(script: Path) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, str(script)],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
        check=False,
    )


def test_ascani_2022_lle_lane_records_accepted_ipopt_phase_split() -> None:
    result = _run(ASCANI_2022 / "scripts" / "run_all.py")

    assert result.returncode == 0, result.stdout + result.stderr
    summary_path = ASCANI_2022 / "results" / "electrolyte_lle" / "summary.json"
    payload = json.loads(summary_path.read_text(encoding="utf-8"))
    assert payload["status"] == "accepted"
    assert payload["solve"]["runtime_ipopt"]["status"] == "enabled_available"
    assert payload["solve"]["accepted"] is True
    assert payload["solve"]["solver_backend"] == "ipopt"

    diagnostics = payload["solve"]["diagnostics"]
    assert diagnostics["route_status"] == "accepted"
    assert diagnostics["solver_status"] in {"success", "acceptable_point"}
    assert diagnostics["application_status"] in {"solve_succeeded", "solved_to_acceptable_level"}
    assert diagnostics["material_balance_norm"] <= payload["expected"]["material_balance_abs"]
    assert diagnostics["charge_balance_norm"] <= payload["expected"]["charge_balance_abs"]
    assert diagnostics["phase_distance"] >= payload["expected"]["phase_distance_min"]

    org_phase, aq_phase = diagnostics["phase_compositions"]
    org_butanol = org_phase[1]
    aq_butanol = aq_phase[1]
    org_ion_fraction = org_phase[2] + org_phase[3]
    aq_ion_fraction = aq_phase[2] + aq_phase[3]
    assert org_butanol > aq_butanol
    assert aq_ion_fraction > org_ion_fraction

    expected_outputs = (
        ASCANI_2022 / "figures" / "figure_4b" / "output" / "figure_4b.svg",
        ASCANI_2022 / "figures" / "figure_4b" / "output" / "data" / "current_phase_compositions.csv",
        ASCANI_2022 / "figures" / "table_5" / "output" / "table_5_fugacity.svg",
        ASCANI_2022 / "figures" / "table_5" / "output" / "data" / "table_5_fugacity_comparison.csv",
        ASCANI_2022 / "figures" / "gibbs_summary" / "output" / "gibbs_summary.svg",
        ASCANI_2022 / "figures" / "gibbs_summary" / "output" / "data" / "gibbs_comparison.csv",
    )
    for output in expected_outputs:
        assert output.exists(), output


def test_ascani_2023_reactive_lane_records_missing_source_targets_without_passing() -> None:
    result = _run(ASCANI_2023 / "scripts" / "run_all.py")

    assert result.returncode == 1, result.stdout + result.stderr
    summary_path = ASCANI_2023 / "results" / "reactive_phase_equilibrium" / "summary.json"
    payload = json.loads(summary_path.read_text(encoding="utf-8"))
    assert payload["status"] == "blocked"
    assert payload["blocker"]["kind"] == "missing_source_target_rows"
    assert "Table 4. Obtained" in payload["source_markers_present"]
    assert "toy reactive LLE fixtures" in payload["not_substituted"]
