from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
RUN_ALL = REPO_ROOT / "analyses" / "data_validation" / "mea_co2_pressure_speciation" / "scripts" / "run_all.py"
SUMMARY_JSON = (
    REPO_ROOT
    / "analyses"
    / "data_validation"
    / "mea_co2_pressure_speciation"
    / "results"
    / "pressure_speciation"
    / "summary.json"
)


def test_mea_co2_pressure_speciation_lane_records_public_api_gate_status() -> None:
    result = subprocess.run(
        [sys.executable, str(RUN_ALL)],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
        check=False,
    )

    assert SUMMARY_JSON.exists(), result.stdout + result.stderr
    summary = json.loads(SUMMARY_JSON.read_text(encoding="utf-8"))

    assert summary["lane"] == "mea_co2_pressure_speciation"
    assert summary["public_api"] == "epcsaft.solve_reactive_speciation"
    assert summary["required_solver_backend"] == "ipopt"
    assert summary["required_derivative_backend"] == "cppad_implicit"
    assert all("MEA-Thermodynamics" not in path for path in summary["source_paths"])

    if summary["status"] == "accepted_public_native_ipopt":
        assert result.returncode == 0
        assert summary["accepted_native_ipopt_speciation"] is True
        diagnostics = summary["solve"]["diagnostics"]
        assert diagnostics["selected_solver_backend"] == "native_ipopt"
        assert diagnostics["derivative_backend"] == "cppad_implicit"
        assert summary["solve"]["pressure_comparison"]["pressure_model"] == "liquid_fugacity_with_ideal_vapor_side"
        assert summary["pressure_model"] == "liquid_co2_fugacity_with_ideal_vapor_side"
        return

    assert result.returncode != 0
    assert summary["status"] == "blocked_capability"
    assert summary["accepted_native_ipopt_speciation"] is False
    assert summary["required_activity_source"] == "liquid ePC-SAFT state at true-species composition"
    assert summary["pressure_model"] == "liquid_co2_fugacity_with_ideal_vapor_side"
    blocker_text = json.dumps(summary, sort_keys=True)
    assert (
        "SSM/DS Born composition sensitivity" in blocker_text
        or "reactive_speciation requires a native Ipopt" in blocker_text
    )
