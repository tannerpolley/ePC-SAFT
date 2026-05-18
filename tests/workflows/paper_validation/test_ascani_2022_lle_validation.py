from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import yaml

REPO_ROOT = Path(__file__).resolve().parents[3]
ASCANI_2022 = REPO_ROOT / "analyses" / "paper_validation" / "native" / "2022_ascani"


def _run(script: Path) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, str(script)],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
        check=False,
    )


def test_ascani_2022_full_case2_public_native_ipopt_lle_validation() -> None:
    result = _run(ASCANI_2022 / "scripts" / "run_all.py")

    assert result.returncode == 0, result.stdout + result.stderr
    summary_path = ASCANI_2022 / "results" / "electrolyte_lle" / "summary.json"
    payload = json.loads(summary_path.read_text(encoding="utf-8"))

    assert payload["schema_version"] == 1
    assert payload["stage"] == "A-C"
    assert payload["lane_id"] == "ascani_2022_distributed_ion_lle"
    assert payload["status"] == "accepted_public_native_ipopt"
    assert payload["feed"]["species"] == ["H2O", "Butanol", "Na+", "K+", "Cl-"]
    assert payload["ipopt_runtime_diagnostics"]["solver_backend"] == "ipopt"
    assert payload["derivative_diagnostics"]["derivative_backend"] == "cppad_implicit"
    assert payload["hessian_approximation_diagnostics"]["hessian_approximation"] == "limited-memory"
    assert payload["density_diagnostics"]["density_backend"] == "liquid_pressure_root"

    tolerances = payload["resolved_tolerances"]
    residuals = payload["material_charge_fugacity_residuals"]
    assert residuals["material_balance_norm"] <= tolerances["material_balance_abs"]
    assert residuals["charge_balance_norm"] <= tolerances["charge_balance_abs"]
    assert residuals["neutral_fugacity_residual_norm"] <= tolerances["neutral_fugacity_abs"]
    assert residuals["salt_pair_fugacity_residual_norm"] <= tolerances["salt_pair_fugacity_abs"]
    assert payload["ghat_delta"] < tolerances["ghat_delta_max"]
    assert payload["tpdf_stability_results"]["accepted"] is True
    assert payload["tpdf_stability_results"]["feed_unstable"] is True
    assert payload["tpdf_stability_results"]["final_phases_stable"] is True

    for output in payload["retained_outputs"]:
        assert (REPO_ROOT / output).exists(), output


def test_ascani_2022_tolerance_contracts_match_yaml_summary_and_registry() -> None:
    from scripts.benchmarks.helpers.literature import LITERATURE_CASES

    analysis = yaml.safe_load((ASCANI_2022 / "analysis.yaml").read_text(encoding="utf-8"))
    payload = json.loads((ASCANI_2022 / "results" / "electrolyte_lle" / "summary.json").read_text(encoding="utf-8"))
    registry = LITERATURE_CASES["ascani_2022_distributed_ion_lle"]

    assert payload["resolved_tolerances"] == analysis["tolerances"]
    assert registry.tolerances == analysis["tolerances"]
    assert payload["expected"]["status"] == analysis["expected"]["status"]
    assert registry.expected["status"] == analysis["expected"]["status"]
