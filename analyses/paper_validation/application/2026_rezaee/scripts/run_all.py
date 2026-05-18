from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path
from typing import Any

from _paths import ANALYSIS_DIR, REPO_ROOT
from scripts.dev.native_runtime_env import apply_to_current_process

SCRIPT_DIR = Path(__file__).resolve().parent
INPUT_DIR = ANALYSIS_DIR / "data" / "input"
PROCESSED_DIR = ANALYSIS_DIR / "data" / "processed"
RESULTS_DIR = ANALYSIS_DIR / "results"
REACTION_RESULTS_DIR = RESULTS_DIR / "reaction_equilibrium"
SUMMARY_JSON = REACTION_RESULTS_DIR / "summary.json"

SCRIPT_SEQUENCE = (
    "rezaee_des_epcsaft_parameter_smoke.py",
    "rezaee_2025_target_summary.py",
    "rezaee_reactive_equilibrium_replay.py",
    "rezaee_reactive_convention_scan.py",
    "rezaee_reactive_epcsaft_option_scan.py",
    "rezaee_paper_basis_reaction_coordinate.py",
    "rezaee_section32_basis_inference.py",
    "rezaee_section32_equilibrium_replication.py",
    "generate_rezaee_2026_figure_comparison_data.py",
    "render_rezaee_2026_figure_comparisons.py",
)

apply_to_current_process()


def _rel(path: Path) -> str:
    return path.relative_to(REPO_ROOT).as_posix()


def _read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _run_script(script_name: str) -> dict[str, Any]:
    command = [sys.executable, str(SCRIPT_DIR / script_name)]
    print("Running:", " ".join(command), flush=True)
    completed = subprocess.run(
        command,
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
        check=False,
    )
    if completed.stdout:
        print(completed.stdout, end="" if completed.stdout.endswith("\n") else "\n", flush=True)
    if completed.stderr:
        print(completed.stderr, end="" if completed.stderr.endswith("\n") else "\n", file=sys.stderr, flush=True)
    if completed.returncode != 0:
        raise subprocess.CalledProcessError(completed.returncode, command, completed.stdout, completed.stderr)
    return {
        "script": script_name,
        "returncode": completed.returncode,
    }


def _direct_closure_supported(section32: dict[str, Any], convention: dict[str, Any]) -> bool:
    direct = section32["direct_held2014_table9_pH_stoich"]
    paper = section32["paper_reference_AARD_pct"]["after_table9_kij"]
    return (
        float(direct["Li_extraction_AARD_pct"]) <= float(paper["Li_extraction"])
        and float(direct["selectivity_AARD_pct"]) <= float(paper["selectivity"])
        and bool(convention["acceptance"]["closed_by_simple_convention_scan"])
    )


def _build_summary(command_results: list[dict[str, Any]]) -> dict[str, Any]:
    replay = _read_json(REACTION_RESULTS_DIR / "rezaee_2026_reactive_equilibrium_replay_summary.json")
    convention = _read_json(REACTION_RESULTS_DIR / "rezaee_2026_reactive_convention_scan_summary.json")
    option_scan = _read_json(REACTION_RESULTS_DIR / "rezaee_2026_reactive_epcsaft_option_scan_summary.json")
    paper_basis = _read_json(REACTION_RESULTS_DIR / "rezaee_2026_paper_basis_reaction_coordinate_summary.json")
    basis = _read_json(REACTION_RESULTS_DIR / "rezaee_2026_section32_basis_inference_summary.json")
    section32 = _read_json(REACTION_RESULTS_DIR / "rezaee_2026_section32_equilibrium_replication_summary.json")
    figure_comparison = _read_json(REACTION_RESULTS_DIR / "rezaee_2026_package_figure_comparison_summary.json")

    closure_supported = _direct_closure_supported(section32, convention)
    direct = section32["direct_held2014_table9_pH_stoich"]
    cross_phase = replay["package_phase_tagged_cross_phase"]
    source_files = sorted(_rel(path) for path in INPUT_DIR.glob("*") if path.is_file())
    retained_outputs = (
        PROCESSED_DIR / "rezaee_2025_extraction_target_summary.csv",
        PROCESSED_DIR / "rezaee_2025_extraction_equilibrium_summary.csv",
        PROCESSED_DIR / "rezaee_2026_reactive_equilibrium_replay.csv",
        PROCESSED_DIR / "rezaee_2026_reactive_convention_scan_rows.csv",
        PROCESSED_DIR / "rezaee_2026_reactive_epcsaft_option_scan_rows.csv",
        PROCESSED_DIR / "rezaee_2026_paper_basis_reaction_coordinate_rows.csv",
        PROCESSED_DIR / "rezaee_2026_section32_basis_inference_rows.csv",
        PROCESSED_DIR / "rezaee_2026_section32_equilibrium_replication_rows.csv",
        REACTION_RESULTS_DIR / "rezaee_2026_reactive_equilibrium_replay_summary.json",
        REACTION_RESULTS_DIR / "rezaee_2026_reactive_convention_scan_summary.json",
        REACTION_RESULTS_DIR / "rezaee_2026_reactive_epcsaft_option_scan_summary.json",
        REACTION_RESULTS_DIR / "rezaee_2026_paper_basis_reaction_coordinate_summary.json",
        REACTION_RESULTS_DIR / "rezaee_2026_section32_basis_inference_summary.json",
        REACTION_RESULTS_DIR / "rezaee_2026_section32_equilibrium_replication_summary.json",
        REACTION_RESULTS_DIR / "rezaee_2026_package_figure_comparison_data_summary.json",
        REACTION_RESULTS_DIR / "rezaee_2026_package_figure_comparison_summary.json",
        REACTION_RESULTS_DIR / "rezaee_2026_package_figure_comparison.md",
    )
    figure_outputs = [
        ANALYSIS_DIR / entry["png"]
        for entry in figure_comparison["figures"]
    ] + [
        ANALYSIS_DIR / entry["svg"]
        for entry in figure_comparison["figures"]
    ] + [
        ANALYSIS_DIR / entry["data"]
        for entry in figure_comparison["figures"]
    ]
    return {
        "status": "source_backed_diagnostic_complete",
        "validation_lane": "rezaee_2026_application_diagnostic",
        "row_count": int(replay["row_count"]),
        "source_files": source_files,
        "source_text_mismatch": {
            "available_si_equilibrium_rows": int(replay["row_count"]),
            "paper_text_equilibrium_data_points": 36,
            "conclusion": (
                "Use the 26 source-backed SI equilibrium rows as the benchmark basis until a 36-row "
                "source table is supplied."
            ),
        },
        "convention_scan": {
            "status": convention["status"],
            "best_variant": convention["best_variant"],
            "source_supported_variant": convention["source_supported_variant"],
            "acceptance": convention["acceptance"],
        },
        "option_scan": {
            "status": option_scan["status"],
            "best_source_supported_variant": option_scan.get("best_source_supported_variant"),
        },
        "paper_basis": {
            "status": paper_basis["status"],
            "row_count": paper_basis["row_count"],
            "best_mode_by_li_abs_error": paper_basis["best_mode_by_li_abs_error"],
        },
        "basis_inference": {
            "status": basis["status"],
            "row_count": basis["row_count"],
            "interpretation": basis["interpretation"],
        },
        "direct_closure": {
            "supported": closure_supported,
            "direct_case_id": direct["case_id"],
            "converged_rows": int(direct["converged_rows"]),
            "li_extraction_aard_pct": float(direct["Li_extraction_AARD_pct"]),
            "selectivity_aard_pct": float(direct["selectivity_AARD_pct"]),
            "paper_reference_aard_pct": section32["paper_reference_AARD_pct"]["after_table9_kij"],
        },
        "figure_comparisons": {
            "status": figure_comparison["status"],
            "figure_count": int(figure_comparison["figure_count"]),
            "figures": figure_comparison["figures"],
        },
        "residual_metrics": {
            "max_abs_charge_residual": float(replay["max_abs_charge_residual"]),
            "max_phase_charge_balance_norm": float(cross_phase["max_phase_charge_balance_norm"]),
            "max_element_balance_norm": float(cross_phase["max_element_balance_norm"]),
            "source_supported_combined_median_abs_ln_residual": float(
                convention["source_supported_variant"]["combined_median_abs_ln_residual"]
            ),
        },
        "retained_outputs": [_rel(path) for path in retained_outputs + tuple(figure_outputs)],
        "commands": command_results,
        "conclusion": (
            "Direct published-constant closure is not supported by the current source-backed inputs; "
            "the executable validation gate is the reproducible diagnostic evidence for that gap."
        ),
    }


def main() -> int:
    command_results = [_run_script(script_name) for script_name in SCRIPT_SEQUENCE]
    summary = _build_summary(command_results)
    SUMMARY_JSON.parent.mkdir(parents=True, exist_ok=True)
    SUMMARY_JSON.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
