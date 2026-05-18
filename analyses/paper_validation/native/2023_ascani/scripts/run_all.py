from __future__ import annotations

import json
import math
import sys
from pathlib import Path
from typing import Any

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[5]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts.dev.native_runtime_env import apply_to_current_process

apply_to_current_process()

import epcsaft

ANALYSIS_DIR = REPO_ROOT / "analyses" / "paper_validation" / "native" / "2023_ascani"
SOURCE_MD = (
    REPO_ROOT
    / "docs"
    / "papers"
    / "md"
    / "Ascani - 2023 - Simultaneous Predictions of Chemical and Phase Equilibria in Systems with an Esterif.md"
)
RESULTS_DIR = ANALYSIS_DIR / "results" / "reactive_phase_equilibrium"
SUMMARY_JSON = RESULTS_DIR / "summary.json"

STRICT_STATUSES = {
    "accepted_public_native_ipopt",
    "blocked_source_data",
    "blocked_solver",
    "blocked_capability",
    "failed_gate",
    "not_started",
}


def _rel(path: Path) -> str:
    return path.relative_to(REPO_ROOT).as_posix()


def _json_like(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): _json_like(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_like(item) for item in value]
    if isinstance(value, np.ndarray):
        return _json_like(value.tolist())
    if isinstance(value, (np.bool_, bool)):
        return bool(value)
    if isinstance(value, (np.integer, int)):
        return int(value)
    if isinstance(value, (np.floating, float)):
        number = float(value)
        return number if math.isfinite(number) else str(number)
    return value


def _phase_options_summary(options: epcsaft.EquilibriumOptions) -> dict[str, Any]:
    return {
        "max_iterations": int(options.max_iterations),
        "tolerance": float(options.tolerance),
        "min_composition": float(options.min_composition),
        "timeout_seconds": float(options.timeout_seconds),
    }


def _status_from_public_route(result: epcsaft.EquilibriumResult) -> str:
    diagnostics = result.diagnostics
    if (
        result.backend == "native_equilibrium_nlp"
        and diagnostics.get("solver_backend") == "ipopt"
        and diagnostics.get("derivative_backend") == "cppad_implicit"
        and diagnostics.get("density_backend") == "liquid_pressure_root"
    ):
        return "accepted_public_native_ipopt"
    return "failed_gate"


def _hypothetical_mixture() -> epcsaft.ePCSAFTMixture:
    return epcsaft.ePCSAFTMixture.from_params(
        {
            "MW": np.asarray([0.050, 0.050, 0.100]),
            "m": np.asarray([2.4, 1.08, 2.8]),
            "s": np.asarray([3.2, 3.0, 3.8]),
            "e": np.asarray([200.0, 400.0, 280.0]),
            "e_assoc": np.asarray([2500.0, 2500.0, 0.0]),
            "vol_a": np.asarray([0.05, 0.05, 0.0]),
            "assoc_scheme": ["2B", "2B", None],
            "k_ij": np.asarray(
                [
                    [0.0, -0.045, 0.045],
                    [-0.045, 0.0, -0.025],
                    [0.045, -0.025, 0.0],
                ],
                dtype=float,
            ),
            "z": np.zeros(3),
            "dielc": np.asarray([20.0, 20.0, 20.0]),
        },
        species=["A", "B", "C"],
    )


def _hypothetical_reaction(standard_state: str = "mole_fraction_activity") -> epcsaft.ReactionDefinition:
    return epcsaft.ReactionDefinition.from_literature_constant(
        {"A": -1.0, "B": -1.0, "C": 1.0},
        log_equilibrium_constant=math.log(2.25),
        name="ascani_2023_hypothetical_A_plus_B_to_C",
        standard_state=standard_state,
        source="Ascani 2023 Table 3 and Figure 3",
    )


def _run_hypothetical_homogeneous_ce() -> dict[str, Any]:
    mix = _hypothetical_mixture()
    feed = np.asarray([0.35, 0.35, 0.30], dtype=float)
    balances = {"A_total": {"A": 1.0, "C": 1.0}, "B_total": {"B": 1.0, "C": 1.0}}
    totals = {"A_total": float(feed[0] + feed[2]), "B_total": float(feed[1] + feed[2])}
    try:
        result = epcsaft.solve_reactive_speciation(
            species=["A", "B", "C"],
            mixture_factory=lambda x, T, P: mix,
            T=300.0,
            P=1.0e5,
            balances=balances,
            totals=totals,
            reactions=[_hypothetical_reaction("ideal_mole_fraction")],
            initial_x=feed,
            options=epcsaft.ReactiveSpeciationOptions(
                solver_backend="ipopt",
                jacobian_backend="cppad",
                tolerance=1.0e-7,
                max_iterations=80,
                activity_output="always",
            ),
        )
    except Exception as exc:
        return {
            "status": "blocked_solver",
            "exception_type": type(exc).__name__,
            "exception_message": str(exc),
        }
    return {
        "status": "accepted_public_native_ipopt" if result.success else "failed_gate",
        "public_api": "epcsaft.solve_reactive_speciation",
        "solver_backend": result.diagnostics.get("selected_solver_backend"),
        "derivative_backend": result.diagnostics.get("implicit_sensitivity_backend"),
        "reaction_residual_norm": result.diagnostics.get("reaction_residual_norm"),
        "composition": result.x,
    }


def _run_hypothetical_reactive_lle() -> dict[str, Any]:
    mix = _hypothetical_mixture()
    feed = np.asarray([0.20, 0.70, 0.10], dtype=float)
    balances = {"A_total": {"A": 1.0, "C": 1.0}, "B_total": {"B": 1.0, "C": 1.0}}
    totals = {"A_total": float(feed[0] + feed[2]), "B_total": float(feed[1] + feed[2])}
    phase_options = epcsaft.EquilibriumOptions(
        max_iterations=350,
        tolerance=2.0e-6,
        min_composition=1.0e-12,
        timeout_seconds=20.0,
    )
    try:
        result = mix.equilibrium(
            kind="reactive_lle",
            T=300.0,
            P=1.0e5,
            z=feed,
            balances=balances,
            totals=totals,
            reactions=[_hypothetical_reaction()],
            phase_options=phase_options,
        )
    except Exception as exc:
        return {
            "status": "blocked_solver",
            "public_api": 'mix.equilibrium(kind="reactive_lle")',
            "exception_type": type(exc).__name__,
            "exception_message": str(exc),
            "feed": feed.tolist(),
            "solver_options": _phase_options_summary(phase_options),
        }
    return {
        "status": _status_from_public_route(result),
        "public_api": 'mix.equilibrium(kind="reactive_lle")',
        "feed": feed.tolist(),
        "solver_options": _phase_options_summary(phase_options),
        "phase_distance": result.diagnostics.get("phase_distance"),
        "reaction_stationarity_norm": result.diagnostics.get("reaction_stationarity_norm"),
        "phase_equilibrium_norm": result.diagnostics.get("phase_equilibrium_norm"),
        "solver_status": result.diagnostics.get("solver_status"),
        "application_status": result.diagnostics.get("application_status"),
        "phases": [phase.to_dict() for phase in result.phases],
    }


def _system1_mixture() -> epcsaft.ePCSAFTMixture:
    temperature = 318.15
    water_sigma = 2.7927 + (10.11 * math.exp(-0.01775 * temperature) - 1.417 * math.exp(-0.01146 * temperature))
    species = ["Water", "Acetic acid", "1-Pentanol", "Pentyl Acetate"]
    return epcsaft.ePCSAFTMixture.from_params(
        {
            "MW": np.asarray([0.01801528, 0.060052, 0.088148, 0.130185]),
            "m": np.asarray([1.2047, 1.3402, 3.6260, 4.7077]),
            "s": np.asarray([water_sigma, 3.8582, 3.4508, 3.4729]),
            "e": np.asarray([353.95, 311.59, 247.28, 234.57]),
            "e_assoc": np.asarray([2425.7, 3044.4, 2252.1, 0.0]),
            "vol_a": np.asarray([0.04509, 0.07555, 0.01033, 0.04509]),
            "assoc_scheme": ["2B", "2B", "2B", None],
            "k_ij": np.asarray(
                [
                    [0.0, -0.1247, 0.001604, -0.0228],
                    [-0.1247, 0.0, -0.1, -0.1],
                    [0.001604, -0.1, 0.0, -0.0095],
                    [-0.0228, -0.1, -0.0095, 0.0],
                ],
                dtype=float,
            ),
            "z": np.zeros(4),
            "dielc": np.asarray([78.3, 6.2, 13.9, 5.0]),
        },
        species=species,
    )


def _run_system1_reactive_lle() -> dict[str, Any]:
    mix = _system1_mixture()
    feed = np.asarray([0.25, 0.25, 0.25, 0.25], dtype=float)
    balances = {
        "acetyl": {"Acetic acid": 1.0, "Pentyl Acetate": 1.0},
        "pentyl": {"1-Pentanol": 1.0, "Pentyl Acetate": 1.0},
        "water_oxygen_pool": {"Water": 1.0, "Pentyl Acetate": 1.0},
    }
    totals = {
        "acetyl": float(feed[1] + feed[3]),
        "pentyl": float(feed[2] + feed[3]),
        "water_oxygen_pool": float(feed[0] + feed[3]),
    }
    reaction = epcsaft.ReactionDefinition.from_literature_constant(
        {"Acetic acid": -1.0, "1-Pentanol": -1.0, "Pentyl Acetate": 1.0, "Water": 1.0},
        log_equilibrium_constant=math.log(43.99),
        name="ascani_2023_system1_esterification",
        standard_state="mole_fraction_activity",
        source="Ascani 2023 Table 4",
    )
    try:
        phase_options = epcsaft.EquilibriumOptions(
            max_iterations=500,
            tolerance=3.0e-6,
            min_composition=1.0e-12,
            timeout_seconds=30.0,
        )
        result = mix.equilibrium(
            kind="reactive_lle",
            T=318.15,
            P=101300.0,
            z=feed,
            balances=balances,
            totals=totals,
            reactions=[reaction],
            phase_options=phase_options,
        )
    except Exception as exc:
        return {
            "status": "blocked_solver",
            "public_api": 'mix.equilibrium(kind="reactive_lle")',
            "exception_type": type(exc).__name__,
            "exception_message": str(exc),
            "feed": feed.tolist(),
            "solver_options": _phase_options_summary(phase_options),
            "paper_match_status": "blocked_source_data",
        }
    return {
        "status": _status_from_public_route(result),
        "public_api": 'mix.equilibrium(kind="reactive_lle")',
        "feed": feed.tolist(),
        "solver_options": _phase_options_summary(phase_options),
        "phase_distance": result.diagnostics.get("phase_distance"),
        "reaction_stationarity_norm": result.diagnostics.get("reaction_stationarity_norm"),
        "phase_equilibrium_norm": result.diagnostics.get("phase_equilibrium_norm"),
        "solver_status": result.diagnostics.get("solver_status"),
        "application_status": result.diagnostics.get("application_status"),
        "phases": [phase.to_dict() for phase in result.phases],
        "paper_match_status": "blocked_source_data",
    }


def _overall_status(*statuses: str) -> str:
    if all(status == "accepted_public_native_ipopt" for status in statuses):
        return "accepted_public_native_ipopt"
    if any(status == "failed_gate" for status in statuses):
        return "failed_gate"
    if any(status == "blocked_solver" for status in statuses):
        return "blocked_solver"
    if any(status == "blocked_capability" for status in statuses):
        return "blocked_capability"
    if any(status == "blocked_source_data" for status in statuses):
        return "blocked_source_data"
    return "not_started"


def main() -> int:
    source_text = SOURCE_MD.read_text(encoding="utf-8", errors="replace")
    required_markers = (
        "Table 1. PC-SAFT pure-component parameters",
        "Table 2. Binary interaction parameters",
        "Table 3. PC-SAFT pure-component parameters of hypothetical mixture",
        "Table 4. Obtained",
    )
    present_markers = [marker for marker in required_markers if marker in source_text]
    homogeneous = _run_hypothetical_homogeneous_ce()
    hypothetical_lle = _run_hypothetical_reactive_lle()
    system1 = _run_system1_reactive_lle()
    stage_d_status = _overall_status(str(homogeneous["status"]), str(hypothetical_lle["status"]))
    stage_e_status = str(system1["status"])
    status = _overall_status(stage_d_status, stage_e_status)
    summary = {
        "schema_version": 1,
        "stage": "D-E",
        "lane_id": "ascani_2023_reactive_phase_equilibrium",
        "status": status,
        "status_reason": (
            "source-backed public native Ipopt reactive LLE route attempt remains blocked by solver closure"
            if status == "blocked_solver"
            else "strict route gates evaluated"
        ),
        "source_records": [_rel(SOURCE_MD)],
        "source_markers_present": present_markers,
        "stage_d_reactive_foundation": {
            "status": stage_d_status,
            "homogeneous_ce_branch": homogeneous,
            "hypothetical_reactive_tie_line": hypothetical_lle,
            "source_basis": "Ascani 2023 Table 3, Figure 3, K_a=2.25",
        },
        "stage_e_ascani_2023": {
            "status": stage_e_status,
            "system_1_route_gate": system1,
            "source_basis": {
                "reaction": "acetic acid + 1-pentanol <=> pentyl acetate + water",
                "temperature_K": 318.15,
                "pressure_bar": 1.013,
                "K_a": 43.99,
                "parameter_tables": ["Table 1", "Table 2", "Table 4"],
            },
            "paper_match_status": "blocked_source_data",
        },
        "solver_contract": {
            "solver_backend": "ipopt",
            "derivative_backend": "cppad_implicit",
            "density_backend": "liquid_pressure_root",
            "hessian_approximation": "limited-memory",
            "phase_volumes_are_nlp_variables": False,
            "exact_gradient_required": True,
            "exact_jacobian_required": True,
        },
        "allowed_statuses": sorted(STRICT_STATUSES),
        "blockers": [
            item
            for item in (hypothetical_lle if hypothetical_lle["status"] != "accepted_public_native_ipopt" else None,
                         system1 if system1["status"] != "accepted_public_native_ipopt" else None)
            if item is not None
        ],
        "claim_boundary": (
            "No Ascani 2023 Figure 8/9 or tie-line paper-match claim is proven. "
            "The retained route attempts use source constants and public APIs only."
        ),
    }
    if summary["status"] not in STRICT_STATUSES:
        raise AssertionError(f"non-strict status generated: {summary['status']}")
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    SUMMARY_JSON.write_text(json.dumps(_json_like(summary), indent=2) + "\n", encoding="utf-8")
    print(json.dumps(_json_like(summary), indent=2))
    return 0 if summary["status"] == "accepted_public_native_ipopt" else 1


if __name__ == "__main__":
    raise SystemExit(main())
