"""Native regression benchmark harness for fixed-shape production contracts."""

from __future__ import annotations

import json
import statistics
import time
from collections import OrderedDict
from collections.abc import Callable, Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import epcsaft


@dataclass(frozen=True)
class PreparedNativeRegressionCase:
    case: str
    description: str
    records: tuple[dict[str, Any], ...]
    parameters: tuple[dict[str, Any], ...]


def _record(
    row_id: str,
    family: str,
    target: str,
    predicted: float,
    observed: float,
    *,
    scale: float = 1.0,
    success: bool = True,
    sensitivities: Mapping[str, float] | None = None,
) -> dict[str, Any]:
    return {
        "row_id": row_id,
        "family": family,
        "target": target,
        "name": f"{row_id}.{family}.{target}",
        "predicted": float(predicted),
        "observed": float(observed),
        "scale": float(scale),
        "sensitivities": {str(key): float(value) for key, value in (sensitivities or {}).items()},
        "success": bool(success),
        "recoverable_failure": False,
        "failure_message": "",
    }


def _parameter(name: str, kind: str, initial: float, lower: float, upper: float) -> dict[str, Any]:
    return {
        "name": name,
        "kind": kind,
        "initial": float(initial),
        "lower": float(lower),
        "upper": float(upper),
        "scale": max(abs(float(initial)), 1.0),
    }


def _native_neutral_tiny() -> PreparedNativeRegressionCase:
    return PreparedNativeRegressionCase(
        case="native_neutral_density_tiny",
        description="single-row neutral density residual contract",
        records=(
            _record("neutral_1", "density", "rho_liq", 998.3, 997.9, scale=0.01, sensitivities={"water.sigma": 12.0}),
            _record(
                "neutral_1",
                "fugacity",
                "ln_phi_water",
                -0.002,
                0.0,
                scale=10.0,
                sensitivities={"water.epsilon_k": -1.0e-4},
            ),
        ),
        parameters=(
            _parameter("water.sigma", "pure_component", 2.79, 2.5, 3.2),
            _parameter("water.epsilon_k", "pure_component", 353.95, 250.0, 450.0),
        ),
    )


def _native_binary_kij_tiny() -> PreparedNativeRegressionCase:
    return PreparedNativeRegressionCase(
        case="native_binary_kij_tiny",
        description="binary interaction pressure/activity residual contract",
        records=(
            _record(
                "binary_1",
                "pressure",
                "P",
                101400.0,
                101325.0,
                scale=1.0e-5,
                sensitivities={"water:Na+.k_ij": 200.0},
            ),
            _record(
                "binary_1",
                "activity",
                "gamma_water",
                1.04,
                1.02,
                scale=10.0,
                sensitivities={"water:Na+.k_ij": 0.6, "water:Cl-.k_ij": -0.2},
            ),
            _record(
                "binary_1",
                "activity",
                "gamma_salt",
                0.93,
                0.95,
                scale=10.0,
                sensitivities={"water:Cl-.k_ij": 0.5},
            ),
        ),
        parameters=(
            _parameter("water:Na+.k_ij", "binary_interaction", 0.0045, -0.5, 0.5),
            _parameter("water:Cl-.k_ij", "binary_interaction", -0.25, -0.8, 0.5),
        ),
    )


def _native_reactive_born_ssmds_tiny() -> PreparedNativeRegressionCase:
    return PreparedNativeRegressionCase(
        case="native_reactive_born_ssmds_tiny",
        description="reactive electrolyte pressure/speciation residual contract with Born SSM+DS d_born and f_solv parameters",
        records=(
            _record(
                "reactive_1",
                "pressure",
                "P_total",
                101900.0,
                101325.0,
                scale=1.0e-5,
                sensitivities={"water.f_solv": 40.0},
            ),
            _record(
                "reactive_1",
                "speciation",
                "x_Na+",
                0.0102,
                0.0100,
                scale=100.0,
                sensitivities={"Na+.d_born": -2.0e-4},
            ),
            _record(
                "reactive_1",
                "speciation",
                "x_Cl-",
                0.0098,
                0.0100,
                scale=100.0,
                sensitivities={"Na+.d_born": 1.0e-4},
            ),
            _record(
                "reactive_1",
                "activity",
                "gamma_Na+",
                0.83,
                0.85,
                scale=10.0,
                sensitivities={"Na+.d_born": 0.02, "water.f_solv": -0.1},
            ),
        ),
        parameters=(
            _parameter("Na+.d_born", "born_radius", 3.445, 2.0, 6.0),
            _parameter("water.f_solv", "solvation_factor", 1.5, 0.5, 3.0),
        ),
    )


def _native_mea_35_row_surrogate() -> PreparedNativeRegressionCase:
    records: list[dict[str, Any]] = []
    for idx in range(35):
        row_id = f"mea_surrogate_{idx + 1:02d}"
        loading_shift = idx / 34.0
        records.append(
            _record(
                row_id,
                "pressure",
                "P_CO2",
                2400.0 + 120.0 * loading_shift,
                2380.0 + 110.0 * loading_shift,
                scale=1.0e-4,
                sensitivities={"MEA.f_solv": 20.0 + 5.0 * loading_shift},
            )
        )
        records.append(
            _record(
                row_id,
                "speciation",
                "x_MEAH+",
                0.08 + 0.02 * loading_shift,
                0.079 + 0.021 * loading_shift,
                scale=100.0,
                sensitivities={"MEAH+.d_born": -1.5e-3, "carbamate.log_k": 2.0e-3},
            )
        )
        records.append(
            _record(
                row_id,
                "speciation",
                "x_MEACOO-",
                0.04 + 0.015 * loading_shift,
                0.041 + 0.014 * loading_shift,
                scale=100.0,
                sensitivities={"MEAH+.d_born": 1.0e-3, "carbamate.log_k": -1.0e-3},
            )
        )
    return PreparedNativeRegressionCase(
        case="native_mea_pressure_speciation_35_row_surrogate",
        description="public MEA-style 35-row pressure/speciation residual contract",
        records=tuple(records),
        parameters=(
            _parameter("MEAH+.d_born", "born_radius", 3.2, 2.0, 6.0),
            _parameter("MEA.f_solv", "solvation_factor", 1.5, 0.5, 3.0),
            _parameter("carbamate.log_k", "reaction_equilibrium_constant", -3.0, -12.0, 4.0),
        ),
    )


CASE_BUILDERS: Mapping[str, Callable[[], PreparedNativeRegressionCase]] = OrderedDict(
    (
        ("native_neutral_density_tiny", _native_neutral_tiny),
        ("native_binary_kij_tiny", _native_binary_kij_tiny),
        ("native_reactive_born_ssmds_tiny", _native_reactive_born_ssmds_tiny),
        ("native_mea_pressure_speciation_35_row_surrogate", _native_mea_35_row_surrogate),
    )
)
DEFAULT_CASES: tuple[str, ...] = tuple(CASE_BUILDERS)


def _run_case_once(prepared: PreparedNativeRegressionCase) -> dict[str, Any]:
    return epcsaft.solve_native_regression_residual_records(
        list(prepared.records),
        list(prepared.parameters),
        options={"derivative_backend": "analytic", "optimizer_backend": "auto"},
    )


def _benchmark_case(prepared: PreparedNativeRegressionCase, *, warmup: int, repeat: int) -> dict[str, Any]:
    if warmup < 0 or repeat <= 0:
        raise ValueError("warmup must be nonnegative and repeat must be positive.")
    for _ in range(warmup):
        _run_case_once(prepared)
    timings: list[int] = []
    last_result: dict[str, Any] | None = None
    for _ in range(repeat):
        start = time.perf_counter_ns()
        last_result = _run_case_once(prepared)
        timings.append(time.perf_counter_ns() - start)
    assert last_result is not None
    objective = last_result["objective_result"]
    families = sorted({str(record["family"]) for record in prepared.records})
    parameter_kinds = sorted({str(parameter["kind"]) for parameter in prepared.parameters})
    return {
        "case": prepared.case,
        "description": prepared.description,
        "warmup": warmup,
        "repeat": repeat,
        "median_ns": int(statistics.median(timings)),
        "mean_ns": int(statistics.mean(timings)),
        "min_ns": min(timings),
        "max_ns": max(timings),
        "status": last_result["status"],
        "success": bool(last_result["success"]),
        "optimizer_backend": last_result["optimizer_backend"],
        "derivative_backend": last_result["derivative_backend"],
        "row_count": len({str(record["row_id"]) for record in prepared.records}),
        "residual_count": len(prepared.records),
        "parameter_count": len(prepared.parameters),
        "target_families": families,
        "parameter_kinds": parameter_kinds,
        "success_count": int(objective["success_count"]),
        "failure_count": int(objective["failure_count"]),
        "fixed_shape_residuals": bool(objective["fixed_shape_residuals"]),
        "residual_norm": float(last_result["residual_norm"]),
        "cost": float(last_result["final_cost"]),
        "production_finite_difference_allowed": False,
    }


def run_native_regression_benchmarks(
    *,
    warmup: int = 1,
    repeat: int = 3,
    case: str | None = None,
    baseline_json: str | Path | None = None,
) -> dict[str, Any]:
    cases = [case] if case is not None else list(DEFAULT_CASES)
    unknown = sorted(set(cases) - set(CASE_BUILDERS))
    if unknown:
        raise ValueError(f"Unknown native regression benchmark case(s): {', '.join(unknown)}")
    baseline = _load_baseline(baseline_json)
    payload = {"benchmark": "native_regression", "warmup": warmup, "repeat": repeat, "cases": []}
    for name in cases:
        case_payload = _benchmark_case(CASE_BUILDERS[name](), warmup=warmup, repeat=repeat)
        if name in baseline:
            baseline_median = int(baseline[name])
            case_payload["baseline_median_ns"] = baseline_median
            case_payload["speedup_vs_baseline"] = baseline_median / max(case_payload["median_ns"], 1)
        payload["cases"].append(case_payload)
    return payload


def _load_baseline(path: str | Path | None) -> dict[str, int]:
    if path is None:
        return {}
    payload = json.loads(Path(path).read_text(encoding="utf-8"))
    return {str(case["case"]): int(case["median_ns"]) for case in payload.get("cases", [])}


def render_benchmark_table(payload: Mapping[str, Any]) -> str:
    lines = ["case median_ms rows residuals params status backend"]
    for case in payload.get("cases", []):
        lines.append(
            " ".join(
                [
                    str(case["case"]),
                    f"{int(case['median_ns']) / 1.0e6:.3f}",
                    str(case["row_count"]),
                    str(case["residual_count"]),
                    str(case["parameter_count"]),
                    str(case["status"]),
                    str(case["optimizer_backend"]),
                ]
            )
        )
    return "\n".join(lines)
