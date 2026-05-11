"""Benchmark harness for the native Ceres thermodynamic regression slice."""

from __future__ import annotations

import math
import statistics
import time
from typing import Any

import numpy as np

import epcsaft


def _salt_speciation_mixture() -> epcsaft.ePCSAFTMixture:
    params = {
        "m": np.asarray([1.2047, 1.0, 1.0, 1.0]),
        "s": np.asarray([2.7927, 3.0, 2.8232, 2.7560]),
        "e": np.asarray([353.95, 200.0, 230.0, 170.0]),
        "z": np.asarray([0.0, 0.0, 1.0, -1.0]),
        "dielc": np.asarray([78.09, 8.0, 8.0, 8.0]),
        "d_born": np.asarray([0.0, 0.0, 3.445, 4.1]),
        "MW": np.asarray([18.01528e-3, 58.44e-3, 22.989e-3, 35.45e-3]),
        "k_ij": np.zeros((4, 4)),
    }
    return epcsaft.ePCSAFTMixture.from_params(params, species=["H2O", "NaCl", "Na+", "Cl-"])


def _case_payload() -> tuple[epcsaft.ePCSAFTMixture, dict[str, Any]]:
    initial_x = [0.998, 0.001, 0.0005, 0.0005]
    log_k = math.log(initial_x[2]) + math.log(initial_x[3]) - math.log(initial_x[1])
    row = {
        "row_id": "salt_speciation",
        "row_mode": "reactive_speciation",
        "T": 298.15,
        "P": 1.0e5,
        "initial_x": initial_x,
        "balance_matrix": [1, 0, 0, 0, 0, 1, 1, 0, 0, 1, 0, 1],
        "balance_rows": 3,
        "total_vector": [0.998, 0.0015, 0.0015],
        "reaction_stoichiometry": [0, -1, 1, 1],
        "reaction_rows": 1,
        "log_equilibrium_constants": [log_k - 1.0],
        "reaction_standard_states": [1],
        "options": {"jacobian_backend": "auto", "max_iterations": 50, "tolerance": 1.0e-10},
        "targets": [{"family": "speciation", "target": "Na+", "index": 2, "observed": 0.00065, "scale": 1000.0}],
    }
    request = {
        "species": ["H2O", "NaCl", "Na+", "Cl-"],
        "rows": [row],
        "parameters": [
            {
                "name": "salt.logK",
                "kind": "reaction_equilibrium_constant",
                "initial": log_k - 1.0,
                "lower": log_k - 5.0,
                "upper": log_k + 5.0,
                "metadata": {"row_id": "salt_speciation", "reaction_index": "0"},
            }
        ],
        "options": {"max_iterations": 20, "derivative_backend": "implicit"},
    }
    return _salt_speciation_mixture(), request


def run_native_ceres_thermo_regression_benchmark(*, warmup: int = 1, repeat: int = 3) -> dict[str, Any]:
    if warmup < 0 or repeat <= 0:
        raise ValueError("warmup must be nonnegative and repeat must be positive.")
    mixture, request = _case_payload()
    for _ in range(warmup):
        epcsaft.fit_native_thermo_regression(mixture, request)

    timings: list[int] = []
    last_result: dict[str, Any] | None = None
    for _ in range(repeat):
        start = time.perf_counter_ns()
        last_result = epcsaft.fit_native_thermo_regression(mixture, request)
        timings.append(time.perf_counter_ns() - start)
    assert last_result is not None

    return {
        "benchmark": "native_ceres_thermo_regression",
        "case": "reactive_speciation_logk_implicit",
        "warmup": warmup,
        "repeat": repeat,
        "median_ns": int(statistics.median(timings)),
        "mean_ns": int(statistics.mean(timings)),
        "min_ns": min(timings),
        "max_ns": max(timings),
        "success": bool(last_result["success"]),
        "status": str(last_result["status"]),
        "optimizer_backend": str(last_result["optimizer_backend"]),
        "derivative_backend": str(last_result["derivative_backend"]),
        "native_hot_loop": str(last_result["optimizer_backend"]) == "ceres",
        "python_objective_used": False,
        "finite_difference_used": False,
        "initial_cost": float(last_result["initial_cost"]),
        "final_cost": float(last_result["final_cost"]),
        "objective_decreased": float(last_result["final_cost"]) < float(last_result["initial_cost"]),
        "message": str(last_result["message"]),
    }


def render_benchmark_table(payload: dict[str, Any]) -> str:
    return (
        "case median_ms status backend derivative native_hot_loop initial_cost final_cost\n"
        f"{payload['case']} {payload['median_ns'] / 1.0e6:.3f} {payload['status']} "
        f"{payload['optimizer_backend']} {payload['derivative_backend']} "
        f"{payload['native_hot_loop']} {payload['initial_cost']:.6g} {payload['final_cost']:.6g}"
    )
