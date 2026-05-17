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


def test_ascani_2022_lle_lane_records_solver_rejection_without_passing() -> None:
    result = _run(ASCANI_2022 / "scripts" / "run_all.py")

    assert result.returncode == 1, result.stdout + result.stderr
    summary_path = ASCANI_2022 / "results" / "electrolyte_lle" / "summary.json"
    payload = json.loads(summary_path.read_text(encoding="utf-8"))
    assert payload["status"] == "blocked"
    assert payload["solve"]["runtime_ipopt"]["status"] == "enabled_available"
    assert payload["solve"]["blocker"] == {
        "kind": "native_ipopt_solver_rejected",
        "route_status": "solver_rejected",
        "solver_status": "ipopt_status_16",
    }
    assert payload["solve"]["accepted"] is False


def test_ascani_2023_reactive_lane_records_missing_source_targets_without_passing() -> None:
    result = _run(ASCANI_2023 / "scripts" / "run_all.py")

    assert result.returncode == 1, result.stdout + result.stderr
    summary_path = ASCANI_2023 / "results" / "reactive_phase_equilibrium" / "summary.json"
    payload = json.loads(summary_path.read_text(encoding="utf-8"))
    assert payload["status"] == "blocked"
    assert payload["blocker"]["kind"] == "missing_source_target_rows"
    assert "Table 4. Obtained" in payload["source_markers_present"]
    assert "toy reactive LLE fixtures" in payload["not_substituted"]
