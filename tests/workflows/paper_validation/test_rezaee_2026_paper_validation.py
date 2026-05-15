from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
ANALYSIS = REPO_ROOT / "analyses" / "paper_validation" / "application" / "2026_rezaee"
SCRIPTS = ANALYSIS / "scripts"


def _run_script(name: str) -> None:
    result = subprocess.run(
        [sys.executable, str(SCRIPTS / name)],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, result.stdout + result.stderr


def test_rezaee_source_backed_paper_validation_generates_pre_surrogate_rows() -> None:
    _run_script("rezaee_2025_target_summary.py")
    _run_script("rezaee_reactive_equilibrium_replay.py")
    _run_script("rezaee_section32_basis_inference.py")
    _run_script("rezaee_section32_equilibrium_replication.py")

    target_summary = ANALYSIS / "data" / "processed" / "rezaee_2025_extraction_target_summary.csv"
    replay_summary = ANALYSIS / "results" / "reaction_equilibrium" / "rezaee_2026_reactive_equilibrium_replay_summary.json"
    basis_summary = ANALYSIS / "results" / "reaction_equilibrium" / "rezaee_2026_section32_basis_inference_summary.json"
    section32_rows = ANALYSIS / "data" / "processed" / "rezaee_2026_section32_equilibrium_replication_rows.csv"
    section32_summary = (
        ANALYSIS / "results" / "reaction_equilibrium" / "rezaee_2026_section32_equilibrium_replication_summary.json"
    )

    assert target_summary.exists()
    assert section32_rows.exists()

    replay_payload = json.loads(replay_summary.read_text(encoding="utf-8"))
    package_cross_phase = replay_payload["package_phase_tagged_cross_phase"]
    assert replay_payload["row_count"] == 26
    assert package_cross_phase["evaluated_rows"] == 26
    assert package_cross_phase["reaction_phase_scope"] == "phase_tagged_cross_phase"
    assert package_cross_phase["native_reaction_residual_size"] == 2
    assert package_cross_phase["max_element_balance_norm"] <= 1.0e-10

    basis_payload = json.loads(basis_summary.read_text(encoding="utf-8"))
    assert basis_payload["row_count"] == 26
    assert basis_payload["status"] == "section32_basis_inference_complete"

    section32_payload = json.loads(section32_summary.read_text(encoding="utf-8"))
    assert section32_payload["row_count"] == 26
    assert section32_payload["status"] == "section32_equation_replication_ran"
    assert section32_payload["direct_held2014_table9_pH_stoich"]["converged_rows"] == 26
