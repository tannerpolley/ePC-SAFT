"""Reactive benchmark harness for package-owned tiny electrolyte regression timing."""

from __future__ import annotations

import json
import statistics
import subprocess
import time
from collections import Counter, OrderedDict
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np

import epcsaft

REPO_ROOT = Path(__file__).resolve().parents[3]

SPECIES = ("water", "Na+", "Cl-")
BALANCES: dict[str, dict[str, float]] = {"water": {"water": 1.0}, "sodium": {"Na+": 1.0}, "chloride": {"Cl-": 1.0}}
REACTIONS: tuple = ()
VAPOR_SPECIES = ("water",)


def _ionic_mix_params() -> dict[str, Any]:
    temperature = 298.15
    s_water = 2.7927 + 10.11 * np.exp(-0.01775 * temperature) - 1.417 * np.exp(-0.01146 * temperature)
    return {
        "MW": np.asarray([18.01528e-3, 22.98e-3, 35.45e-3]),
        "m": np.asarray([1.2047, 1.0, 1.0]),
        "s": np.asarray([s_water, 2.8232, 2.7560]),
        "e": np.asarray([353.95, 230.0, 170.0]),
        "e_assoc": np.asarray([2425.7, 0.0, 0.0]),
        "vol_a": np.asarray([0.04509, 0.0, 0.0]),
        "assoc_scheme": ["2B", None, None],
        "z": np.asarray([0.0, 1.0, -1.0]),
        "dielc": np.asarray([78.09, 8.0, 8.0]),
        "d_born": np.asarray([0.0, 3.445, 4.1]),
        "f_solv": np.asarray([1.5, 1.0, 1.0]),
        "k_ij": np.asarray([[0.0, 0.0045, -0.25], [0.0045, 0.0, 0.317], [-0.25, 0.317, 0.0]]),
        "l_ij": np.zeros((3, 3)),
        "k_hb": np.zeros((3, 3)),
    }


@dataclass(frozen=True)
class BenchmarkObservation:
    fingerprint: dict[str, Any]
    fallback_used: bool
    diagnostics: dict[str, Any]
    row_count: int
    parameter_count: int
    success_count: int
    failure_count: int
    residual_count: int
    cache_hits: int
    cache_misses: int
    speciation_solves: int
    bubble_solves: int
    density_solves: int
    activity_calls: int
    fugacity_calls: int
    counter_details: Mapping[str, Any] = field(default_factory=dict)


@dataclass(frozen=True)
class PreparedBenchmarkCase:
    case: str
    description: str
    runner: Callable[[], BenchmarkObservation]
    baseline_runner: Callable[[], BenchmarkObservation] | None = None


def _to_int(value: Any) -> int:
    try:
        return int(value)
    except (TypeError, ValueError):
        return 0


def _to_float(value: Any) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return float("nan")


def _runtime_cache_stats(mix: epcsaft.ePCSAFTMixture | None = None) -> dict[str, int]:
    if mix is None:
        return {
            "reference_state_cache_hits": 0,
            "reference_state_cache_misses": 0,
            "density_warm_start_hits": 0,
            "density_warm_start_fallbacks": 0,
        }
    try:
        payload = mix.runtime_cache_stats()
    except Exception:
        return {
            "reference_state_cache_hits": 0,
            "reference_state_cache_misses": 0,
            "density_warm_start_hits": 0,
            "density_warm_start_fallbacks": 0,
        }
    return {
        "reference_state_cache_hits": _to_int(payload.get("reference_state_cache_hits", 0)),
        "reference_state_cache_misses": _to_int(payload.get("reference_state_cache_misses", 0)),
        "density_warm_start_hits": _to_int(payload.get("density_warm_start_hits", 0)),
        "density_warm_start_fallbacks": _to_int(payload.get("density_warm_start_fallbacks", 0)),
    }


def _cache_delta(start: Mapping[str, int], end: Mapping[str, int], key: str) -> int:
    return int(end.get(key, 0) - start.get(key, 0))


def _fallback_used_from_diagnostics(diagnostics: Mapping[str, Any] | None) -> bool:
    diag = diagnostics or {}
    return bool(
        diag.get("jacobian_fallback_used", False)
        or diag.get("hessian_fallback_used", False)
        or diag.get("state_solver_fallback", False)
        or diag.get("phase_handoff", False)
    )


def _extract_nested_diagnostics(diagnostics: Mapping[str, Any] | None) -> tuple[int, int, list[str]]:
    speciation_diagnostics = _mapping_value(diagnostics, "speciation")
    bubble_diagnostics = _mapping_value(diagnostics, "bubble")
    spec_inner = _mapping_value(speciation_diagnostics, "diagnostics")
    bubble_inner = _mapping_value(bubble_diagnostics, "diagnostics")
    density = (
        _to_int(spec_inner.get("density_solve_count"))
        + _to_int(bubble_inner.get("density_solve_count"))
        + _to_int(speciation_diagnostics.get("density_solve_count"))
    )
    activity = _to_int(spec_inner.get("activity_evaluation_count")) + _to_int(
        speciation_diagnostics.get("activity_evaluation_count")
    )
    keys = []
    for block in (speciation_diagnostics, spec_inner, bubble_diagnostics, bubble_inner):
        keys.extend(str(key) for key in block)
    return density, activity, keys


def _safe_mapping(payload: Mapping[str, Any] | None) -> Mapping[str, Any]:
    return payload if isinstance(payload, Mapping) else {}


def _mapping_value(payload: Mapping[str, Any] | None, key: str) -> Mapping[str, Any]:
    value = _safe_mapping(payload).get(key)
    return value if isinstance(value, Mapping) else {}


def _mixture_factory_and_instance() -> (
    tuple[Callable[[Any, float, float | None], epcsaft.ePCSAFTMixture], epcsaft.ePCSAFTMixture]
):
    mix = epcsaft.ePCSAFTMixture.from_params(_ionic_mix_params(), species=SPECIES)

    def _factory(_x: Any, _t: float, _p: float | None = None) -> epcsaft.ePCSAFTMixture:
        return mix

    return _factory, mix


def _run_solve_reactive_speciation(points: Sequence[Mapping[str, Any]]) -> Callable[[], BenchmarkObservation]:
    factory, mix = _mixture_factory_and_instance()

    def _run() -> BenchmarkObservation:
        cache_before = _runtime_cache_stats(mix)
        success_count = 0
        failure_count = 0
        speciation_solves = 0
        bubble_solves = 0
        density_solves = 0
        activity_calls = 0
        fugacity_calls = 0
        fallback_used = False
        row_fingerprints: list[dict[str, Any]] = []
        diagnostics_keys: list[str] = []
        residual_count = 0

        for point in points:
            try:
                result = epcsaft.solve_reactive_speciation(
                    species=list(SPECIES),
                    mixture_factory=factory,
                    T=float(point["T"]),
                    P=float(point["P"]),
                    balances=BALANCES,
                    totals=dict(point["totals"]),
                    reactions=REACTIONS,
                    initial_x=point["initial_x"],
                    options=epcsaft.ReactiveSpeciationOptions(error_mode="result"),
                )
                speciation_solves += 1
                success_count += 1 if result.success else 0
                failure_count += 0 if result.success else 1
                fallback_used = fallback_used or _fallback_used_from_diagnostics(result.diagnostics)
                d_count, a_count, diag_keys = _extract_nested_diagnostics(result.diagnostics)
                density_solves += d_count
                activity_calls += a_count
                diagnostics_keys.extend(diag_keys)
                residual_count += 0
                row_fingerprints.append(
                    {
                        "success": bool(result.success),
                        "x_liq": {str(key): _to_float(value) for key, value in result.x.items()},
                        "charge_residual": _to_float(result.charge_residual),
                    }
                )
            except Exception as exc:
                failure_count += 1
                row_fingerprints.append({"success": False, "message": f"{type(exc).__name__}: {exc}"})

        cache_after = _runtime_cache_stats(mix)
        return BenchmarkObservation(
            fingerprint={
                "case": "reactive_speciation_batch_tiny",
                "row_count": len(points),
                "parameter_count": len(SPECIES),
                "row_fingerprints": row_fingerprints,
            },
            fallback_used=fallback_used,
            diagnostics={"diagnostics_keys": sorted(set(diagnostics_keys))},
            row_count=len(points),
            parameter_count=len(SPECIES),
            success_count=success_count,
            failure_count=failure_count,
            residual_count=residual_count,
            cache_hits=_cache_delta(cache_before, cache_after, "reference_state_cache_hits"),
            cache_misses=_cache_delta(cache_before, cache_after, "reference_state_cache_misses"),
            speciation_solves=speciation_solves,
            bubble_solves=bubble_solves,
            density_solves=density_solves,
            activity_calls=activity_calls,
            fugacity_calls=fugacity_calls,
            counter_details={
                "native_reference_state_cache_hits": _cache_delta(
                    cache_before, cache_after, "reference_state_cache_hits"
                ),
                "native_reference_state_cache_misses": _cache_delta(
                    cache_before, cache_after, "reference_state_cache_misses"
                ),
                "density_warm_start_hits": _cache_delta(cache_before, cache_after, "density_warm_start_hits"),
                "density_warm_start_fallbacks": _cache_delta(cache_before, cache_after, "density_warm_start_fallbacks"),
            },
        )

    return _run


def _run_solve_reactive_bubble(points: Sequence[Mapping[str, Any]]) -> Callable[[], BenchmarkObservation]:
    factory, mix = _mixture_factory_and_instance()

    def _run() -> BenchmarkObservation:
        cache_before = _runtime_cache_stats(mix)
        success_count = 0
        failure_count = 0
        speciation_solves = 0
        bubble_solves = 0
        density_solves = 0
        activity_calls = 0
        fugacity_calls = 0
        fallback_used = False
        row_fingerprints: list[dict[str, Any]] = []
        diagnostics_keys: list[str] = []
        residual_count = 0

        results = epcsaft.solve_reactive_electrolyte_bubble_sweep(
            species=list(SPECIES),
            mixture_factory=factory,
            points=points,
            balances=BALANCES,
            reactions=REACTIONS,
            vapor_species=VAPOR_SPECIES,
            continuation="none",
            options=epcsaft.ReactiveElectrolyteBubbleOptions(error_mode="result"),
        )
        for point, result in zip(points, results, strict=False):
            if result.success:
                success_count += 1
            else:
                failure_count += 1
            speciation_solves += 1
            bubble_solves += 1
            d_count, a_count, diag_keys = _extract_nested_diagnostics(_safe_mapping(result.diagnostics))
            density_solves += d_count
            activity_calls += a_count
            fugacity_calls += len(_safe_mapping(result.partial_pressures))
            fallback_used = fallback_used or _fallback_used_from_diagnostics(result.diagnostics.get("speciation"))
            diagnostics_keys.extend(diag_keys)
            residual_count += (
                _to_int(result.diagnostics.get("speciation", {}).get("residual_count"))
                if isinstance(result.diagnostics, Mapping)
                else 0
            )
            row_fingerprints.append(
                {
                    "success": bool(result.success),
                    "T": _to_float(point.get("T")),
                    "P_total": _to_float(result.P_total),
                    "y_vap": {str(key): _to_float(value) for key, value in result.y_vap.items()},
                }
            )
        cache_after = _runtime_cache_stats(mix)
        return BenchmarkObservation(
            fingerprint={
                "case": "reactive_bubble_batch_tiny",
                "row_count": len(points),
                "parameter_count": len(SPECIES),
                "vapor_species": list(VAPOR_SPECIES),
                "row_fingerprints": row_fingerprints,
            },
            fallback_used=fallback_used,
            diagnostics={"diagnostics_keys": sorted(set(diagnostics_keys))},
            row_count=len(points),
            parameter_count=len(SPECIES),
            success_count=success_count,
            failure_count=failure_count,
            residual_count=residual_count,
            cache_hits=_cache_delta(cache_before, cache_after, "reference_state_cache_hits"),
            cache_misses=_cache_delta(cache_before, cache_after, "reference_state_cache_misses"),
            speciation_solves=speciation_solves,
            bubble_solves=bubble_solves,
            density_solves=density_solves,
            activity_calls=activity_calls,
            fugacity_calls=fugacity_calls,
            counter_details={
                "native_reference_state_cache_hits": _cache_delta(
                    cache_before, cache_after, "reference_state_cache_hits"
                ),
                "native_reference_state_cache_misses": _cache_delta(
                    cache_before, cache_after, "reference_state_cache_misses"
                ),
                "density_warm_start_hits": _cache_delta(cache_before, cache_after, "density_warm_start_hits"),
                "density_warm_start_fallbacks": _cache_delta(cache_before, cache_after, "density_warm_start_fallbacks"),
            },
        )

    return _run


def _objective_rows() -> tuple[dict[str, Any], ...]:
    return (
        {
            "row_id": "tiny-1",
            "T": 298.15,
            "P_seed": 101325.0,
            "totals": {"water": 0.98, "sodium": 0.01, "chloride": 0.01},
            "initial_x": [0.98, 0.01, 0.01],
            "target_partial_pressures": {"water": 3166.40625},
            "target_x": {"water": 0.98},
        },
        {
            "row_id": "tiny-2",
            "T": 298.15,
            "P_seed": 95000.0,
            "totals": {"water": 0.97, "sodium": 0.015, "chloride": 0.015},
            "initial_x": [0.97, 0.015, 0.015],
            "target_partial_pressures": {"water": 3122.0},
            "target_x": {"water": 0.97},
        },
    )


def _run_objective_legacy(
    rows: Sequence[Mapping[str, Any]],
    *,
    pressure_weight: float,
    speciation_weight: float,
) -> BenchmarkObservation:
    factory, mix = _mixture_factory_and_instance()
    cache_before = _runtime_cache_stats(mix)
    result = epcsaft.evaluate_reactive_electrolyte_bubble_residuals(
        rows,
        species=SPECIES,
        mixture_factory=factory,
        balances=BALANCES,
        reactions=REACTIONS,
        vapor_species=VAPOR_SPECIES,
        reaction_names=("reactive_regression_reaction",),
        pressure_species=VAPOR_SPECIES,
        speciation_species=("water",),
        pressure_weight=pressure_weight,
        speciation_weight=speciation_weight,
        reaction_weight=1.0e-12,
        options=epcsaft.ReactiveElectrolyteBubbleOptions(error_mode="result"),
        continuation="none",
    )
    cache_after = _runtime_cache_stats(mix)
    cache_hits = _cache_delta(cache_before, cache_after, "reference_state_cache_hits")
    cache_misses = _cache_delta(cache_before, cache_after, "reference_state_cache_misses")

    residual_count = int(result.residuals.size)
    diagnostics_keys: list[str] = []
    density_solves = 0
    activity_calls = 0
    fugacity_calls = 0
    speciation_solves = 0
    bubble_solves = 0
    fallback_used = False
    for row in result.record_results:
        row_diag = _mapping_value(dict(row), "diagnostics")
        d_count, a_count, diag_keys = _extract_nested_diagnostics(row_diag)
        density_solves += d_count
        activity_calls += a_count
        fugacity_calls += len(_mapping_value(dict(row), "partial_pressures"))
        speciation_solves += 1
        bubble_solves += 1
        diagnostics_keys.extend(diag_keys)
        fallback_used = fallback_used or _fallback_used_from_diagnostics(row_diag)

    return BenchmarkObservation(
        fingerprint={
            "case": "reactive_regression_objective_tiny",
            "row_count": len(rows),
            "parameter_count": len(SPECIES),
            "residual_norm": _to_float(result.residual_norm),
        },
        fallback_used=fallback_used,
        diagnostics={"diagnostics_keys": sorted(set(diagnostics_keys)), "residual_payload": result.to_dict()},
        row_count=len(rows),
        parameter_count=len(SPECIES),
        success_count=int(result.success_count),
        failure_count=int(result.failure_count),
        residual_count=residual_count,
        cache_hits=cache_hits,
        cache_misses=cache_misses,
        speciation_solves=speciation_solves,
        bubble_solves=bubble_solves,
        density_solves=density_solves,
        activity_calls=activity_calls,
        fugacity_calls=fugacity_calls,
        counter_details={
            "native_reference_state_cache_hits": cache_hits,
            "native_reference_state_cache_misses": cache_misses,
            "density_warm_start_hits": _cache_delta(cache_before, cache_after, "density_warm_start_hits"),
            "density_warm_start_fallbacks": _cache_delta(cache_before, cache_after, "density_warm_start_fallbacks"),
        },
    )


def _compile_objective_context(
    rows: Sequence[Mapping[str, Any]],
    *,
    pressure_weight: float,
    speciation_weight: float,
) -> epcsaft.ReactiveElectrolyteRegressionContext:
    batch_rows = [
        epcsaft.ReactiveElectrolyteRow(
            row_id=str(row["row_id"]),
            T=float(row["T"]),
            P_seed=float(row["P_seed"]),
            totals=dict(row["totals"]),
            initial_x=list(row["initial_x"]),
            balances=BALANCES,
            reactions=REACTIONS,
            vapor_species=VAPOR_SPECIES,
            target_partial_pressures=dict(row["target_partial_pressures"]),
            target_speciation=dict(row["target_x"]),
        )
        for row in rows
    ]
    batch = epcsaft.ReactiveElectrolyteBatch(
        species=SPECIES,
        rows=batch_rows,
        balances=BALANCES,
        reactions=REACTIONS,
        vapor_species=VAPOR_SPECIES,
        base_parameters=_ionic_mix_params(),
        options=epcsaft.ReactiveElectrolyteBatchOptions(
            warm_start_rows=True,
            warm_start_objective=True,
            include_state_outputs=False,
            penalty_value=8.0,
        ),
        reactive_bubble_options=epcsaft.ReactiveElectrolyteBubbleOptions(error_mode="result"),
    )
    objective = epcsaft.build_reactive_regression_objective(
        batch,
        residual_weights={
            "partial_pressure": float(pressure_weight),
            "speciation": float(speciation_weight),
            "reaction": 1.0e-12,
        },
        failure_penalty=8.0,
    )
    return epcsaft.ReactiveElectrolyteRegressionContext.from_batch(
        species=batch.species,
        rows=batch.rows,
        balances=batch.balances,
        reactions=batch.reactions,
        options=batch.options,
        objective=objective,
        vapor_species=batch.vapor_species,
        base_parameters=batch.base_parameters,
        reactive_bubble_options=batch.reactive_bubble_options,
    )


def _run_objective_compiled(
    rows: Sequence[Mapping[str, Any]],
    *,
    pressure_weight: float,
    speciation_weight: float,
    parameter_map: Mapping[str, float] | None = None,
) -> Callable[[], BenchmarkObservation]:
    context = _compile_objective_context(
        rows,
        pressure_weight=pressure_weight,
        speciation_weight=speciation_weight,
    )
    parameter_map = dict(parameter_map or {})

    def _run() -> BenchmarkObservation:
        result = context.evaluate_objective(parameter_map)
        batch_result = result.batch_result
        target_counter = _safe_mapping(batch_result.diagnostics).get("target_family_counts", {})
        solve_counter = _safe_mapping(batch_result.diagnostics).get("solve_counts", {})
        context_cache_hits = int(batch_result.cache_stats.get("context_cache_hits", 0))
        context_cache_misses = int(batch_result.cache_stats.get("context_cache_misses", 0))
        objective_seed_hits = int(batch_result.cache_stats.get("objective_seed_hits", 0))
        objective_seed_misses = int(batch_result.cache_stats.get("objective_seed_misses", 0))
        return BenchmarkObservation(
            fingerprint={
                "case": "reactive_regression_objective_tiny",
                "row_count": len(rows),
                "parameter_count": len(parameter_map),
                "objective": _to_float(result.objective),
                "residual_norm": _to_float(np.linalg.norm(result.residuals)),
            },
            fallback_used=bool(batch_result.failure_count),
            diagnostics={
                "diagnostics_keys": sorted(batch_result.diagnostics.keys()),
                "cache_stats": dict(batch_result.cache_stats),
                "target_family_counts": dict(target_counter),
            },
            row_count=len(rows),
            parameter_count=max(len(parameter_map), 1),
            success_count=batch_result.success_count,
            failure_count=batch_result.failure_count,
            residual_count=int(result.residuals.size),
            cache_hits=context_cache_hits + objective_seed_hits,
            cache_misses=context_cache_misses + objective_seed_misses,
            speciation_solves=int(solve_counter.get("speciation_solves", 0)),
            bubble_solves=int(solve_counter.get("bubble_solves", 0)),
            density_solves=int(solve_counter.get("density_solves", 0)),
            activity_calls=int(solve_counter.get("activity_calls", 0)),
            fugacity_calls=int(target_counter.get("partial_pressure", 0)),
            counter_details={
                "context_cache_hits": context_cache_hits,
                "context_cache_misses": context_cache_misses,
                "objective_seed_hits": objective_seed_hits,
                "objective_seed_misses": objective_seed_misses,
                "native_reference_state_cache_hits": None,
                "native_reference_state_cache_misses": None,
                "density_warm_start_hits": None,
                "density_warm_start_fallbacks": None,
                "unavailable_counters": [
                    "native_reference_state_cache_hits",
                    "native_reference_state_cache_misses",
                    "density_warm_start_hits",
                    "density_warm_start_fallbacks",
                ],
            },
        )

    return _run


def _speciation_surrogate_rows(row_count: int = 35) -> tuple[dict[str, Any], ...]:
    rows: list[dict[str, Any]] = []
    for idx in range(row_count):
        salt = 0.004 + 0.00025 * (idx % 7)
        water = 1.0 - 2.0 * salt
        rows.append(
            {
                "row_id": f"mea-surrogate-{idx + 1:02d}",
                "T": 298.15 + 0.2 * (idx % 5),
                "P": 101325.0,
                "totals": {"water": water, "sodium": salt, "chloride": salt},
                "initial_x": [water, salt, salt],
                "target_x": {"water": water},
                "target_activity": {"water": 1.0},
            }
        )
    return tuple(rows)


def _mixed_pressure_speciation_surrogate_rows(row_count: int = 35) -> tuple[dict[str, Any], ...]:
    if row_count < 3:
        raise ValueError("mixed pressure/speciation surrogate requires at least three rows.")
    rows: list[dict[str, Any]] = []
    for row in _objective_rows():
        mutated = dict(row)
        mutated["mode"] = "bubble"
        mutated["source"] = "public_pressure_surrogate"
        mutated["split"] = "benchmark"
        rows.append(mutated)
    for idx, row in enumerate(_speciation_surrogate_rows(row_count - len(rows))):
        mutated = dict(row)
        mutated["row_id"] = f"mixed-speciation-{idx + 1:02d}"
        mutated["mode"] = "speciation"
        mutated["source"] = "public_speciation_surrogate"
        mutated["split"] = "benchmark"
        rows.append(mutated)
    return tuple(rows)


def _compile_speciation_surrogate_context(
    rows: Sequence[Mapping[str, Any]],
) -> epcsaft.ReactiveElectrolyteRegressionContext:
    batch_rows = [
        epcsaft.ReactiveElectrolyteRow(
            row_id=str(row["row_id"]),
            T=float(row["T"]),
            P=float(row["P"]),
            totals=dict(row["totals"]),
            initial_x=list(row["initial_x"]),
            balances=BALANCES,
            reactions=REACTIONS,
            target_speciation=dict(row["target_x"]),
            target_activity=dict(row["target_activity"]),
            source="public_surrogate",
            split="benchmark",
            metadata={"surrogate": "MEA-style trace-carbonate workflow shape"},
            mode="speciation",
        )
        for row in rows
    ]
    batch = epcsaft.ReactiveElectrolyteBatch(
        species=SPECIES,
        rows=batch_rows,
        balances=BALANCES,
        reactions=REACTIONS,
        base_parameters=_ionic_mix_params(),
        options=epcsaft.ReactiveElectrolyteBatchOptions(
            warm_start_rows=True,
            warm_start_objective=True,
            include_state_outputs=False,
            penalty_value=8.0,
        ),
    )
    objective = epcsaft.build_reactive_regression_objective(
        batch,
        residual_weights={"speciation": 1.0, "activity": 1.0, "reaction": 0.0},
        failure_penalty=8.0,
    )
    return epcsaft.ReactiveElectrolyteRegressionContext.from_batch(
        species=batch.species,
        rows=batch.rows,
        balances=batch.balances,
        reactions=batch.reactions,
        options=batch.options,
        objective=objective,
        base_parameters=batch.base_parameters,
    )


def _compile_mixed_pressure_speciation_surrogate_context(
    rows: Sequence[Mapping[str, Any]],
) -> epcsaft.ReactiveElectrolyteRegressionContext:
    batch_rows: list[epcsaft.ReactiveElectrolyteRow] = []
    for row in rows:
        mode = str(row.get("mode", "speciation"))
        if mode == "bubble":
            batch_rows.append(
                epcsaft.ReactiveElectrolyteRow(
                    row_id=str(row["row_id"]),
                    T=float(row["T"]),
                    P_seed=float(row["P_seed"]),
                    totals=dict(row["totals"]),
                    initial_x=list(row["initial_x"]),
                    balances=BALANCES,
                    reactions=REACTIONS,
                    vapor_species=VAPOR_SPECIES,
                    target_partial_pressures=dict(row["target_partial_pressures"]),
                    target_speciation=dict(row["target_x"]),
                    source=str(row.get("source", "public_pressure_surrogate")),
                    split=str(row.get("split", "benchmark")),
                    metadata={"surrogate": "pressure/speciation mixed public workflow shape"},
                    mode="bubble",
                )
            )
            continue
        batch_rows.append(
            epcsaft.ReactiveElectrolyteRow(
                row_id=str(row["row_id"]),
                T=float(row["T"]),
                P=float(row["P"]),
                totals=dict(row["totals"]),
                initial_x=list(row["initial_x"]),
                balances=BALANCES,
                reactions=REACTIONS,
                target_speciation=dict(row["target_x"]),
                target_activity=dict(row["target_activity"]),
                source=str(row.get("source", "public_speciation_surrogate")),
                split=str(row.get("split", "benchmark")),
                metadata={"surrogate": "pressure/speciation mixed public workflow shape"},
                mode="speciation",
            )
        )
    batch = epcsaft.ReactiveElectrolyteBatch(
        species=SPECIES,
        rows=batch_rows,
        balances=BALANCES,
        reactions=REACTIONS,
        vapor_species=VAPOR_SPECIES,
        base_parameters=_ionic_mix_params(),
        options=epcsaft.ReactiveElectrolyteBatchOptions(
            warm_start_rows=True,
            warm_start_objective=True,
            include_state_outputs=False,
            penalty_value=8.0,
        ),
        reactive_bubble_options=epcsaft.ReactiveElectrolyteBubbleOptions(error_mode="result"),
    )
    objective = epcsaft.build_reactive_regression_objective(
        batch,
        residual_weights={
            "partial_pressure": 1.0,
            "speciation": 1.0,
            "activity": 1.0,
            "reaction": 0.0,
        },
        failure_penalty=8.0,
    )
    return epcsaft.ReactiveElectrolyteRegressionContext.from_batch(
        species=batch.species,
        rows=batch.rows,
        balances=batch.balances,
        reactions=batch.reactions,
        options=batch.options,
        objective=objective,
        vapor_species=batch.vapor_species,
        base_parameters=batch.base_parameters,
        reactive_bubble_options=batch.reactive_bubble_options,
    )


def _run_speciation_surrogate_compiled(
    rows: Sequence[Mapping[str, Any]],
) -> Callable[[], BenchmarkObservation]:
    context = _compile_speciation_surrogate_context(rows)

    def _run() -> BenchmarkObservation:
        result = context.evaluate_objective({"Na+.sigma": 2.8232})
        batch_result = result.batch_result
        target_counter = _safe_mapping(batch_result.diagnostics).get("target_family_counts", {})
        solve_counter = _safe_mapping(batch_result.diagnostics).get("solve_counts", {})
        context_cache_hits = int(batch_result.cache_stats.get("context_cache_hits", 0))
        context_cache_misses = int(batch_result.cache_stats.get("context_cache_misses", 0))
        objective_seed_hits = int(batch_result.cache_stats.get("objective_seed_hits", 0))
        objective_seed_misses = int(batch_result.cache_stats.get("objective_seed_misses", 0))
        return BenchmarkObservation(
            fingerprint={
                "case": "mea_trace_carbonate_35_row_public_surrogate",
                "row_count": len(rows),
                "parameter_count": 1,
                "objective": _to_float(result.objective),
                "residual_norm": _to_float(np.linalg.norm(result.residuals)),
                "surrogate_note": "Public synthetic rows; same compiled reactive-regression speciation workflow shape.",
            },
            fallback_used=bool(batch_result.failure_count),
            diagnostics={
                "diagnostics_keys": sorted(batch_result.diagnostics.keys()),
                "cache_stats": dict(batch_result.cache_stats),
                "target_family_counts": dict(target_counter),
            },
            row_count=len(rows),
            parameter_count=1,
            success_count=batch_result.success_count,
            failure_count=batch_result.failure_count,
            residual_count=int(result.residuals.size),
            cache_hits=context_cache_hits + objective_seed_hits,
            cache_misses=context_cache_misses + objective_seed_misses,
            speciation_solves=int(solve_counter.get("speciation_solves", 0)),
            bubble_solves=int(solve_counter.get("bubble_solves", 0)),
            density_solves=int(solve_counter.get("density_solves", 0)),
            activity_calls=int(solve_counter.get("activity_calls", 0)),
            fugacity_calls=int(target_counter.get("partial_pressure", 0)),
            counter_details={
                "context_cache_hits": context_cache_hits,
                "context_cache_misses": context_cache_misses,
                "objective_seed_hits": objective_seed_hits,
                "objective_seed_misses": objective_seed_misses,
                "native_reference_state_cache_hits": None,
                "native_reference_state_cache_misses": None,
                "density_warm_start_hits": None,
                "density_warm_start_fallbacks": None,
                "unavailable_counters": [
                    "native_reference_state_cache_hits",
                    "native_reference_state_cache_misses",
                    "density_warm_start_hits",
                    "density_warm_start_fallbacks",
                ],
            },
        )

    return _run


def _run_mixed_pressure_speciation_surrogate_compiled(
    rows: Sequence[Mapping[str, Any]],
) -> Callable[[], BenchmarkObservation]:
    context = _compile_mixed_pressure_speciation_surrogate_context(rows)

    def _run() -> BenchmarkObservation:
        result = context.evaluate_objective({"Na+.sigma": 2.8232})
        batch_result = result.batch_result
        target_counter = _safe_mapping(batch_result.diagnostics).get("target_family_counts", {})
        solve_counter = _safe_mapping(batch_result.diagnostics).get("solve_counts", {})
        context_cache_hits = int(batch_result.cache_stats.get("context_cache_hits", 0))
        context_cache_misses = int(batch_result.cache_stats.get("context_cache_misses", 0))
        objective_seed_hits = int(batch_result.cache_stats.get("objective_seed_hits", 0))
        objective_seed_misses = int(batch_result.cache_stats.get("objective_seed_misses", 0))
        return BenchmarkObservation(
            fingerprint={
                "case": "reactive_regression_pressure_speciation_35_row_surrogate",
                "row_count": len(rows),
                "parameter_count": 1,
                "objective": _to_float(result.objective),
                "residual_norm": _to_float(np.linalg.norm(result.residuals)),
                "surrogate_note": "Public synthetic rows with both bubble-pressure and speciation residual families.",
                "target_family_counts": dict(target_counter),
            },
            fallback_used=bool(batch_result.failure_count),
            diagnostics={
                "diagnostics_keys": sorted(batch_result.diagnostics.keys()),
                "cache_stats": dict(batch_result.cache_stats),
                "target_family_counts": dict(target_counter),
            },
            row_count=len(rows),
            parameter_count=1,
            success_count=batch_result.success_count,
            failure_count=batch_result.failure_count,
            residual_count=int(result.residuals.size),
            cache_hits=context_cache_hits + objective_seed_hits,
            cache_misses=context_cache_misses + objective_seed_misses,
            speciation_solves=int(solve_counter.get("speciation_solves", 0)),
            bubble_solves=int(solve_counter.get("bubble_solves", 0)),
            density_solves=int(solve_counter.get("density_solves", 0)),
            activity_calls=int(solve_counter.get("activity_calls", 0)),
            fugacity_calls=int(target_counter.get("partial_pressure", 0)),
            counter_details={
                "context_cache_hits": context_cache_hits,
                "context_cache_misses": context_cache_misses,
                "objective_seed_hits": objective_seed_hits,
                "objective_seed_misses": objective_seed_misses,
                "native_reference_state_cache_hits": None,
                "native_reference_state_cache_misses": None,
                "density_warm_start_hits": None,
                "density_warm_start_fallbacks": None,
                "unavailable_counters": [
                    "native_reference_state_cache_hits",
                    "native_reference_state_cache_misses",
                    "density_warm_start_hits",
                    "density_warm_start_fallbacks",
                ],
            },
        )

    return _run


def _case_builder_reactive_speciation_tiny() -> PreparedBenchmarkCase:
    points = (
        {
            "T": 298.15,
            "P": 101325.0,
            "totals": {"water": 0.98, "sodium": 0.01, "chloride": 0.01},
            "initial_x": [0.98, 0.01, 0.01],
        },
        {
            "T": 298.15,
            "P": 95000.0,
            "totals": {"water": 0.97, "sodium": 0.015, "chloride": 0.015},
            "initial_x": [0.97, 0.015, 0.015],
        },
    )
    return PreparedBenchmarkCase(
        case="reactive_speciation_batch_tiny",
        description="Reactive speciation batch over two ionic-water points.",
        runner=_run_solve_reactive_speciation(points),
    )


def _case_builder_reactive_bubble_tiny() -> PreparedBenchmarkCase:
    points = (
        {
            "T": 298.15,
            "P_seed": 101325.0,
            "totals": {"water": 0.98, "sodium": 0.01, "chloride": 0.01},
            "initial_x": [0.98, 0.01, 0.01],
        },
        {
            "T": 298.15,
            "P_seed": 95000.0,
            "totals": {"water": 0.97, "sodium": 0.015, "chloride": 0.015},
            "initial_x": [0.97, 0.015, 0.015],
        },
    )
    return PreparedBenchmarkCase(
        case="reactive_bubble_batch_tiny",
        description="Reactive electrolyte bubble batch over two ionic-water points.",
        runner=_run_solve_reactive_bubble(points),
    )


def _case_builder_reactive_regression_objective_tiny() -> PreparedBenchmarkCase:
    rows = _objective_rows()
    return PreparedBenchmarkCase(
        case="reactive_regression_objective_tiny",
        description="Reactive electrolyte pressure/speciation objective over a tiny two-row batch.",
        runner=_run_objective_compiled(rows, pressure_weight=1.0, speciation_weight=1.0),
        baseline_runner=lambda: _run_objective_legacy(rows, pressure_weight=1.0, speciation_weight=1.0),
    )


def _case_builder_reactive_regression_parameter_perturbation() -> PreparedBenchmarkCase:
    base_rows = list(_objective_rows())
    perturbed_rows = []
    for row in base_rows:
        mutated = dict(row)
        target = dict(mutated["target_partial_pressures"])
        for species_key, value in target.items():
            target[species_key] = float(value) * 1.001
        mutated["target_partial_pressures"] = target
        mutated["row_id"] = f"{mutated['row_id']}-perturbed"
        perturbed_rows.append(mutated)
    rows = tuple(base_rows + perturbed_rows)

    return PreparedBenchmarkCase(
        case="reactive_regression_parameter_perturbation",
        description="Tiny reaction-objective perturbation with baseline/perturbed pressure targets.",
        runner=_run_objective_compiled(
            rows,
            pressure_weight=1.0,
            speciation_weight=1.0,
            parameter_map={"Na+.sigma": 2.8232 * 1.001},
        ),
        baseline_runner=lambda: _run_objective_legacy(rows, pressure_weight=1.0, speciation_weight=1.0),
    )


def _case_builder_mea_trace_carbonate_35_row_public_surrogate() -> PreparedBenchmarkCase:
    rows = _speciation_surrogate_rows(35)
    return PreparedBenchmarkCase(
        case="mea_trace_carbonate_35_row_public_surrogate",
        description="Public 35-row MEA-style reactive-regression surrogate using native speciation solves.",
        runner=_run_speciation_surrogate_compiled(rows),
    )


def _case_builder_reactive_regression_pressure_speciation_35_row_surrogate() -> PreparedBenchmarkCase:
    rows = _mixed_pressure_speciation_surrogate_rows(35)
    return PreparedBenchmarkCase(
        case="reactive_regression_pressure_speciation_35_row_surrogate",
        description="Public 35-row reactive-regression surrogate with pressure and speciation residual families.",
        runner=_run_mixed_pressure_speciation_surrogate_compiled(rows),
    )


CASE_BUILDERS: OrderedDict[str, Callable[[], PreparedBenchmarkCase]] = OrderedDict(
    (
        ("reactive_speciation_batch_tiny", _case_builder_reactive_speciation_tiny),
        ("reactive_bubble_batch_tiny", _case_builder_reactive_bubble_tiny),
        ("reactive_regression_objective_tiny", _case_builder_reactive_regression_objective_tiny),
        ("reactive_regression_parameter_perturbation", _case_builder_reactive_regression_parameter_perturbation),
        (
            "reactive_regression_pressure_speciation_35_row_surrogate",
            _case_builder_reactive_regression_pressure_speciation_35_row_surrogate,
        ),
        (
            "mea_trace_carbonate_35_row_public_surrogate",
            _case_builder_mea_trace_carbonate_35_row_public_surrogate,
        ),
    )
)
DEFAULT_CASES: tuple[str, ...] = (
    "reactive_speciation_batch_tiny",
    "reactive_bubble_batch_tiny",
    "reactive_regression_objective_tiny",
    "reactive_regression_parameter_perturbation",
)


def _git_commit() -> str | None:
    try:
        completed = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=REPO_ROOT,
            capture_output=True,
            text=True,
            check=True,
        )
    except (OSError, subprocess.SubprocessError):
        return None
    commit = completed.stdout.strip()
    return commit or None


def _benchmark_case(prepared: PreparedBenchmarkCase, *, warmup: int, repeat: int) -> dict[str, Any]:
    for _ in range(max(0, warmup)):
        prepared.runner()

    timings_ns: list[int] = []
    success_count = 0
    failure_count = 0
    residual_count = 0
    cache_hits = 0
    cache_misses = 0
    speciation_solves = 0
    bubble_solves = 0
    density_solves = 0
    activity_calls = 0
    fugacity_calls = 0
    measured_success_repeat_count = 0
    failure_messages: list[str] = []
    row_count = 0
    parameter_count = len(SPECIES)
    fallback_used = False
    fingerprint: dict[str, Any] | None = None
    diagnostics_keys: list[str] = []
    target_family_counts: Counter[str] = Counter()
    counter_details_totals: dict[str, Any] = {
        "context_cache_hits": 0,
        "context_cache_misses": 0,
        "objective_seed_hits": 0,
        "objective_seed_misses": 0,
        "native_reference_state_cache_hits": 0,
        "native_reference_state_cache_misses": 0,
        "density_warm_start_hits": 0,
        "density_warm_start_fallbacks": 0,
        "unavailable_counters": [],
    }

    for _ in range(max(1, repeat)):
        start_ns = time.perf_counter_ns()
        try:
            observation = prepared.runner()
        except Exception as exc:
            failure_count += 1
            failure_messages.append(f"{type(exc).__name__}: {exc}")
            continue
        elapsed_ns = time.perf_counter_ns() - start_ns
        timings_ns.append(elapsed_ns)
        measured_success_repeat_count += 1

        success_count += observation.success_count
        failure_count += observation.failure_count
        residual_count += observation.residual_count
        cache_hits += observation.cache_hits
        cache_misses += observation.cache_misses
        speciation_solves += observation.speciation_solves
        bubble_solves += observation.bubble_solves
        density_solves += observation.density_solves
        activity_calls += observation.activity_calls
        fugacity_calls += observation.fugacity_calls
        fallback_used = fallback_used or bool(observation.fallback_used)
        row_count = observation.row_count
        parameter_count = observation.parameter_count
        diagnostics_keys.extend(observation.diagnostics.get("diagnostics_keys", []))
        target_family_counts.update(
            {
                str(key): int(value)
                for key, value in _safe_mapping(observation.diagnostics.get("target_family_counts")).items()
            }
        )
        for key, value in observation.counter_details.items():
            if value is None:
                unavailable = counter_details_totals.setdefault("unavailable_counters", [])
                if key not in unavailable:
                    unavailable.append(key)
                counter_details_totals[key] = None
            elif isinstance(value, list):
                unavailable = counter_details_totals.setdefault("unavailable_counters", [])
                for item in value:
                    if item not in unavailable:
                        unavailable.append(item)
            elif counter_details_totals.get(key) is not None:
                counter_details_totals[key] = int(counter_details_totals.get(key, 0)) + int(value)
        if fingerprint is None:
            fingerprint = observation.fingerprint

    if not timings_ns:
        raise RuntimeError(f"Reactive benchmark case {prepared.case} failed for every measured repeat.")

    timings = np.asarray(timings_ns, dtype=np.int64)
    payload = {
        "case": prepared.case,
        "description": prepared.description,
        "warmup": int(warmup),
        "repeat": int(repeat),
        "row_count": int(row_count),
        "parameter_count": int(parameter_count),
        "success_count": int(success_count),
        "failure_count": int(failure_count),
        "measured_success_repeat_count": int(measured_success_repeat_count),
        "failure_messages": failure_messages,
        "residual_count": int(residual_count),
        "median_ns": int(np.median(timings)),
        "mean_ns": round(statistics.fmean(timings_ns)),
        "min_ns": int(np.min(timings)),
        "max_ns": int(np.max(timings)),
        "p10_ns": int(np.percentile(timings, 10.0)),
        "p90_ns": int(np.percentile(timings, 90.0)),
        "fallback_used": bool(fallback_used),
        "fingerprint": fingerprint or {},
        "diagnostics_keys": sorted(set(diagnostics_keys)),
        "target_family_counts": dict(target_family_counts),
        "cache_hits": int(cache_hits),
        "cache_misses": int(cache_misses),
        "speciation_solves": int(speciation_solves),
        "bubble_solves": int(bubble_solves),
        "density_solves": int(density_solves),
        "activity_calls": int(activity_calls),
        "fugacity_calls": int(fugacity_calls),
        "context_cache_hits": counter_details_totals["context_cache_hits"],
        "context_cache_misses": counter_details_totals["context_cache_misses"],
        "objective_seed_hits": counter_details_totals["objective_seed_hits"],
        "objective_seed_misses": counter_details_totals["objective_seed_misses"],
        "native_reference_state_cache_hits": counter_details_totals["native_reference_state_cache_hits"],
        "native_reference_state_cache_misses": counter_details_totals["native_reference_state_cache_misses"],
        "density_warm_start_hits": counter_details_totals["density_warm_start_hits"],
        "density_warm_start_fallbacks": counter_details_totals["density_warm_start_fallbacks"],
        "unavailable_counters": sorted(set(counter_details_totals["unavailable_counters"])),
    }
    if prepared.baseline_runner is not None:
        baseline_timings_ns: list[int] = []
        baseline_repeat = min(max(1, repeat), 3)
        baseline_warmup = min(max(0, warmup), 1)
        for _ in range(baseline_warmup):
            prepared.baseline_runner()
        for _ in range(baseline_repeat):
            start_ns = time.perf_counter_ns()
            prepared.baseline_runner()
            baseline_timings_ns.append(time.perf_counter_ns() - start_ns)
        baseline_timings = np.asarray(baseline_timings_ns, dtype=np.int64)
        baseline_median = int(np.median(baseline_timings))
        payload["baseline_median_ns"] = baseline_median
        payload["baseline_repeat"] = int(baseline_repeat)
        payload["baseline_warmup"] = int(baseline_warmup)
        payload["speedup_vs_baseline"] = float(baseline_median / max(int(payload["median_ns"]), 1))
    return payload


def _augment_with_baseline(rows: list[dict[str, Any]], baseline_payload: dict[str, Any] | None) -> list[dict[str, Any]]:
    if baseline_payload is None:
        return rows
    baseline_rows = {str(row["case"]): row for row in baseline_payload.get("cases", [])}
    augmented: list[dict[str, Any]] = []
    for row in rows:
        current = dict(row)
        baseline = baseline_rows.get(str(row["case"]))
        if baseline is not None and "median_ns" in baseline:
            baseline_median = int(baseline["median_ns"])
            current["baseline_median_ns"] = baseline_median
            current["speedup_vs_baseline"] = float(baseline_median / max(int(row["median_ns"]), 1))
        augmented.append(current)
    return augmented


def run_reactive_regression_benchmarks(
    *,
    warmup: int = 3,
    repeat: int = 10,
    case: str | None = None,
    baseline_json: Path | None = None,
) -> dict[str, Any]:
    if warmup < 0:
        raise ValueError("warmup must be non-negative")
    if repeat <= 0:
        raise ValueError("repeat must be positive")
    if case is not None and case not in CASE_BUILDERS:
        raise ValueError(f"Unknown case {case!r}. Expected one of: {', '.join(CASE_BUILDERS)}")

    baseline_payload: dict[str, Any] | None = None
    if baseline_json is not None:
        baseline_payload = json.loads(Path(baseline_json).read_text(encoding="utf-8"))

    selected = [case] if case is not None else list(DEFAULT_CASES)
    rows = [_benchmark_case(CASE_BUILDERS[name](), warmup=warmup, repeat=repeat) for name in selected]
    rows = _augment_with_baseline(rows, baseline_payload)
    build_info = epcsaft.runtime_build_info()
    return {
        "warmup": int(warmup),
        "repeat": int(repeat),
        "selected_cases": selected,
        "package_version": str(epcsaft.__version__),
        "git_commit": _git_commit() or str(epcsaft.__git_commit__),
        "build_info": build_info,
        "cases": rows,
    }


def render_benchmark_table(payload: dict[str, Any]) -> str:
    headers = [
        "case",
        "row_count",
        "median_ms",
        "mean_ms",
        "p10_ms",
        "p90_ms",
        "success",
        "fallback",
    ]
    if any("speedup_vs_baseline" in row for row in payload["cases"]):
        headers.append("speedup")
    widths = {header: len(header) for header in headers}
    rows: list[list[str]] = []

    for row in payload["cases"]:
        values = [
            str(row["case"]),
            str(row["row_count"]),
            f"{row['median_ns'] / 1.0e6:.3f}",
            f"{row['mean_ns'] / 1.0e6:.3f}",
            f"{row['p10_ns'] / 1.0e6:.3f}",
            f"{row['p90_ns'] / 1.0e6:.3f}",
            f"{row['success_count']}",
            "yes" if row["fallback_used"] else "no",
        ]
        if "speedup_vs_baseline" in row:
            values.append(f"{row['speedup_vs_baseline']:.2f}x")
        rows.append(values)
        for header, value in zip(headers, values):
            widths[header] = max(widths[header], len(value))

    def _format(values: list[str]) -> str:
        return "  ".join(value.ljust(widths[header]) for header, value in zip(headers, values))

    header_line = _format(headers)
    divider = "  ".join("-" * widths[header] for header in headers)
    body = "\n".join(_format(values) for values in rows)
    return "\n".join((header_line, divider, body))
