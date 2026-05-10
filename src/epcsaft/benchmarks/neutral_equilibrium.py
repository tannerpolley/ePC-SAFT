from __future__ import annotations

import json
import statistics
import subprocess
import time
from collections import OrderedDict
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

import epcsaft

REPO_ROOT = Path(__file__).resolve().parents[3]


@dataclass(frozen=True)
class BenchmarkObservation:
    fingerprint: dict[str, Any]
    fallback_used: bool
    diagnostics: dict[str, Any]


@dataclass(frozen=True)
class PreparedBenchmarkCase:
    case: str
    description: str
    runner: Callable[[], BenchmarkObservation]


def _hydrocarbon_mixture() -> epcsaft.ePCSAFTMixture:
    params = {
        "m": np.asarray([1.0, 1.6069, 2.0020]),
        "s": np.asarray([3.7039, 3.5206, 3.6184]),
        "e": np.asarray([150.03, 191.42, 208.11]),
        "k_ij": np.asarray(
            [
                [0.0, 3.0e-4, 1.15e-2],
                [3.0e-4, 0.0, 5.10e-3],
                [1.15e-2, 5.10e-3, 0.0],
            ]
        ),
    }
    return epcsaft.ePCSAFTMixture.from_params(params, species=["Methane", "Ethane", "Propane"])


def _methanol_cyclohexane_mixture() -> epcsaft.ePCSAFTMixture:
    params = {
        "MW": np.asarray([32.042e-3, 84.147e-3]),
        "m": np.asarray([1.5255, 2.5303]),
        "s": np.asarray([3.2300, 3.8499]),
        "e": np.asarray([188.90, 278.11]),
        "e_assoc": np.asarray([2899.5, 0.0]),
        "vol_a": np.asarray([0.035176, 0.0]),
        "assoc_scheme": ["2B", None],
        "k_ij": np.asarray([[0.0, 0.051], [0.051, 0.0]]),
        "z": np.asarray([0.0, 0.0]),
        "dielc": np.asarray([33.05, 2.02]),
    }
    return epcsaft.ePCSAFTMixture.from_params(params, species=["Methanol", "Cyclohexane"])


def _round_scalar(value: Any, digits: int = 10) -> float:
    return round(float(value), digits)


def _round_array(values: Any, digits: int = 10) -> list[float]:
    arr = np.asarray(values, dtype=float)
    return [round(float(item), digits) for item in arr.tolist()]


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


def _fallback_used_from_diagnostics(diagnostics: dict[str, Any]) -> bool:
    return bool(
        diagnostics.get("neutral_fallback_used", False)
        or diagnostics.get("density_fallback_used", False)
        or diagnostics.get("jacobian_fallback_used", False)
        or diagnostics.get("finite_difference_fallback_used", False)
        or diagnostics.get("hessian_fallback_used", False)
    )


def _equilibrium_observation(result: epcsaft.EquilibriumResult) -> BenchmarkObservation:
    diagnostics = dict(result.diagnostics)
    fingerprint = {
        "backend": str(result.backend),
        "problem_kind": str(result.problem_kind),
        "phase_labels": list(result.phase_labels),
        "split_detected": bool(result.split_detected),
        "stable": bool(result.stable),
        "phase_fractions": _round_array([phase.phase_fraction for phase in result.phases]),
        "densities": _round_array([phase.density for phase in result.phases]),
        "pressures": _round_array([phase.pressure for phase in result.phases]),
        "temperatures": _round_array([phase.temperature for phase in result.phases]),
        "compositions": [_round_array(phase.composition) for phase in result.phases],
        "ln_phi": [_round_array(phase.ln_fugacity_coefficient) for phase in result.phases],
        "fugacity_residual_norm": _round_scalar(diagnostics.get("fugacity_residual_norm", 0.0)),
    }
    return BenchmarkObservation(
        fingerprint=fingerprint,
        fallback_used=_fallback_used_from_diagnostics(diagnostics),
        diagnostics=diagnostics,
    )


def _state_observation(state: epcsaft.ePCSAFTState) -> BenchmarkObservation:
    diagnostics = dict(getattr(state, "diagnostics", lambda: {})() or {})
    fingerprint = {
        "density": _round_scalar(state.density()),
        "pressure": _round_scalar(state.pressure()),
        "compressibility_factor": _round_scalar(state.compressibility_factor()),
        "ln_phi": _round_array(state.fugacity_coefficient()),
    }
    return BenchmarkObservation(
        fingerprint=fingerprint,
        fallback_used=_fallback_used_from_diagnostics(diagnostics),
        diagnostics=diagnostics,
    )


def _prepare_neutral_state() -> PreparedBenchmarkCase:
    mix = _hydrocarbon_mixture()
    feed = np.asarray([0.1, 0.3, 0.6], dtype=float)

    def runner() -> BenchmarkObservation:
        state = mix.state(T=220.0, P=1.0e5, x=feed, phase="vap")
        return _state_observation(state)

    return PreparedBenchmarkCase(
        case="neutral_state",
        description="Neutral hydrocarbon state density and fugacity coefficients at T=220 K, P=1e5 Pa.",
        runner=runner,
    )


def _prepare_tp_flash() -> PreparedBenchmarkCase:
    mix = _hydrocarbon_mixture()
    feed = np.asarray([0.1, 0.3, 0.6], dtype=float)

    def runner() -> BenchmarkObservation:
        result = mix.flash_tp(T=220.0, P=1.0e5, z=feed)
        return _equilibrium_observation(result)

    return PreparedBenchmarkCase(
        case="tp_flash",
        description="Methane/ethane/propane TP flash at T=220 K, P=1e5 Pa, z=[0.1, 0.3, 0.6].",
        runner=runner,
    )


def _prepare_bubble_p() -> PreparedBenchmarkCase:
    mix = _hydrocarbon_mixture()
    flash = mix.flash_tp(T=220.0, P=1.0e5, z=np.asarray([0.1, 0.3, 0.6], dtype=float))
    liquid = flash.phases[0]

    def runner() -> BenchmarkObservation:
        result = mix.bubble_p(T=220.0, x=liquid.composition)
        return _equilibrium_observation(result)

    return PreparedBenchmarkCase(
        case="bubble_p",
        description="Neutral bubble pressure from the ternary hydrocarbon TP-flash liquid endpoint.",
        runner=runner,
    )


def _prepare_dew_p() -> PreparedBenchmarkCase:
    mix = _hydrocarbon_mixture()
    flash = mix.flash_tp(T=220.0, P=1.0e5, z=np.asarray([0.1, 0.3, 0.6], dtype=float))
    vapor = flash.phases[1]

    def runner() -> BenchmarkObservation:
        result = mix.dew_p(T=220.0, y=vapor.composition)
        return _equilibrium_observation(result)

    return PreparedBenchmarkCase(
        case="dew_p",
        description="Neutral dew pressure from the ternary hydrocarbon TP-flash vapor endpoint.",
        runner=runner,
    )


def _prepare_lle_seeded() -> PreparedBenchmarkCase:
    mix = _methanol_cyclohexane_mixture()
    methanol_poor = np.asarray([0.05, 0.95], dtype=float)
    methanol_rich = np.asarray([0.85, 0.15], dtype=float)
    feed = 0.5 * methanol_poor + 0.5 * methanol_rich
    initial_phases = {"liq1": methanol_poor, "liq2": methanol_rich, "phase_fraction": 0.5}
    options = epcsaft.EquilibriumOptions(
        max_iterations=240,
        tolerance=1.0e-10,
        damping=0.5,
        jacobian_backend="finite_difference",
    )

    def runner() -> BenchmarkObservation:
        result = mix.lle_tp(T=298.15, P=1.013e5, z=feed, initial_phases=initial_phases, options=options)
        return _equilibrium_observation(result)

    return PreparedBenchmarkCase(
        case="lle_seeded",
        description="Seeded methanol/cyclohexane neutral LLE flash at T=298.15 K, P=1.013e5 Pa.",
        runner=runner,
    )


CASE_BUILDERS: OrderedDict[str, Callable[[], PreparedBenchmarkCase]] = OrderedDict(
    (
        ("neutral_state", _prepare_neutral_state),
        ("tp_flash", _prepare_tp_flash),
        ("bubble_p", _prepare_bubble_p),
        ("dew_p", _prepare_dew_p),
        ("lle_seeded", _prepare_lle_seeded),
    )
)


def _benchmark_case(prepared: PreparedBenchmarkCase, *, warmup: int, repeat: int) -> dict[str, Any]:
    for _ in range(warmup):
        prepared.runner()

    timings_ns: list[int] = []
    failures = 0
    fallback_used = False
    fingerprint: dict[str, Any] | None = None
    diagnostics_keys: list[str] = []
    fingerprint_consistent = True
    failure_messages: list[str] = []

    for _ in range(repeat):
        start = time.perf_counter_ns()
        try:
            observation = prepared.runner()
        except Exception as exc:  # pragma: no cover - failure path exercised through JSON payloads if needed
            failures += 1
            failure_messages.append(str(exc))
            continue
        elapsed_ns = time.perf_counter_ns() - start
        timings_ns.append(elapsed_ns)
        fallback_used = fallback_used or bool(observation.fallback_used)
        diagnostics_keys = sorted(set(diagnostics_keys).union(observation.diagnostics.keys()))
        if fingerprint is None:
            fingerprint = observation.fingerprint
        else:
            fingerprint_consistent = fingerprint_consistent and (fingerprint == observation.fingerprint)

    if not timings_ns:
        raise RuntimeError(f"Benchmark case {prepared.case} failed for every measured repeat.")

    timings = np.asarray(timings_ns, dtype=np.int64)
    return {
        "case": prepared.case,
        "description": prepared.description,
        "warmup": int(warmup),
        "repeat": int(repeat),
        "failures": int(failures),
        "fallback_used": bool(fallback_used),
        "median_ns": int(np.median(timings)),
        "mean_ns": int(round(statistics.fmean(timings_ns))),
        "min_ns": int(np.min(timings)),
        "max_ns": int(np.max(timings)),
        "p10_ns": int(np.percentile(timings, 10.0)),
        "p90_ns": int(np.percentile(timings, 90.0)),
        "iqr_ns": int(np.percentile(timings, 75.0) - np.percentile(timings, 25.0)),
        "fingerprint_consistent": bool(fingerprint_consistent),
        "fingerprint": fingerprint or {},
        "diagnostics_keys": diagnostics_keys,
        "failure_messages": failure_messages[:10],
    }


def _augment_with_baseline(
    rows: list[dict[str, Any]],
    baseline_payload: dict[str, Any] | None,
) -> list[dict[str, Any]]:
    if baseline_payload is None:
        return rows
    baseline_rows = {str(row["case"]): row for row in baseline_payload.get("cases", [])}
    augmented: list[dict[str, Any]] = []
    for row in rows:
        current = dict(row)
        baseline = baseline_rows.get(str(row["case"]))
        if baseline is not None:
            baseline_median = int(baseline["median_ns"])
            current["baseline_median_ns"] = baseline_median
            current["speedup_vs_baseline"] = float(baseline_median / max(int(row["median_ns"]), 1))
        augmented.append(current)
    return augmented


def run_neutral_equilibrium_benchmarks(
    *,
    warmup: int = 20,
    repeat: int = 100,
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

    selected = [case] if case is not None else list(CASE_BUILDERS)
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
    headers = ["case", "median_ms", "mean_ms", "p10_ms", "p90_ms", "failures", "fallback"]
    if any("speedup_vs_baseline" in row for row in payload["cases"]):
        headers.append("speedup")
    widths = {header: len(header) for header in headers}
    rows: list[list[str]] = []
    for row in payload["cases"]:
        values = [
            str(row["case"]),
            f"{row['median_ns'] / 1.0e6:.3f}",
            f"{row['mean_ns'] / 1.0e6:.3f}",
            f"{row['p10_ns'] / 1.0e6:.3f}",
            f"{row['p90_ns'] / 1.0e6:.3f}",
            str(row["failures"]),
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
