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


CASE_BUILDERS: OrderedDict[str, Callable[[], PreparedBenchmarkCase]] = OrderedDict(
    (
        ("neutral_state", _prepare_neutral_state),
    )
)


def _benchmark_case(prepared: PreparedBenchmarkCase, *, warmup: int, repeat: int) -> dict[str, Any]:
    for _ in range(warmup):
        prepared.runner()

    timings_ns: list[int] = []
    failures = 0
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
        "median_ns": int(np.median(timings)),
        "mean_ns": round(statistics.fmean(timings_ns)),
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
    headers = ["case", "median_ms", "mean_ms", "p10_ms", "p90_ms", "failures"]
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
