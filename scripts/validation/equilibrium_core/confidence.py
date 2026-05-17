"""Opt-in electrolyte LLE confidence validation helpers.

These helpers intentionally call the public ``ePCSAFTMixture.equilibrium`` facade
so the confidence suite validates the same route a user exercises.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from epcsaft._types import SolutionError
from epcsaft.epcsaft import ePCSAFTMixture
from epcsaft.equilibrium import EquilibriumOptions
from scripts.validation.equilibrium_core.phase_state import phase_state

REPO_ROOT = Path(__file__).resolve().parents[3]
BENCHMARK_ROOT = REPO_ROOT / "data" / "reference" / "equilibrium_benchmarks" / "electrolyte_lle"
KHUIDAIDA_ANALYSIS = REPO_ROOT / "analyses" / "paper_validation" / "application" / "2026_khudaida"
DEFAULT_OUTPUT_ROOT = KHUIDAIDA_ANALYSIS / "results" / "final" / "reports"
PRESSURE_PA = 100000.0
FORMULA_SPECIES = ("H2O", "Ethanol", "Butanol", "NaCl")
EXPLICIT_SPECIES = ("H2O", "Ethanol", "Butanol", "Na+", "Cl-")
CHARGES = np.asarray([0.0, 0.0, 0.0, 1.0, -1.0], dtype=float)
CSV_FLOAT_FORMAT = ".17g"


@dataclass(frozen=True)
class BenchmarkCase:
    suite: str
    case_key: str
    tie_line: int
    figure: int
    temperature_K: float
    salt_wtfrac: float
    feed_formula: np.ndarray
    experimental_organic_formula: np.ndarray
    experimental_aqueous_formula: np.ndarray
    source: str


@dataclass(frozen=True)
class BenchmarkSuite:
    name: str
    species: tuple[str, ...]
    formula_species: tuple[str, ...]
    dataset: str
    thresholds: dict[str, float]
    cases: tuple[BenchmarkCase, ...]


@dataclass(frozen=True)
class NativeBenchmarkInputs:
    feed: np.ndarray


@dataclass(frozen=True)
class BenchmarkMetrics:
    grand_aad: float
    max_abs_error: float
    organic_aad: float
    aqueous_aad: float
    phase_distance: float


@dataclass(frozen=True)
class BenchmarkPrediction:
    case: BenchmarkCase
    outcome: str
    diagnostics: dict[str, Any]
    metrics: BenchmarkMetrics | None
    organic_formula: np.ndarray | None
    aqueous_formula: np.ndarray | None
    beta_org: float | None


@dataclass(frozen=True)
class ContinuationMetrics:
    series_key: str
    from_case: str
    to_case: str
    composition_jump_norm: float
    beta_jump: float
    phase_label_swapped: bool
    branch_collapse_flag: bool


@dataclass(frozen=True)
class OracleCheck:
    case_key: str
    native_outcome: str
    native_gibbs_delta: float
    oracle_gibbs_delta: float
    native_minus_oracle_gibbs_delta: float
    oracle_flags_native: bool


@dataclass(frozen=True)
class ConfidenceReport:
    output_dir: Path
    summary_path: Path
    benchmark_csv: Path
    continuation_csv: Path
    oracle_csv: Path
    stress_csv: Path
    residual_gate_plot: Path
    error_plot: Path
    all_tielines_plot: Path
    continuation_plot: Path


def load_benchmark_suite(name: str = "khudaida_2026") -> BenchmarkSuite:
    if name != "khudaida_2026":
        raise ValueError("Only the khudaida_2026 electrolyte LLE benchmark suite is available.")
    root = BENCHMARK_ROOT / name
    metadata = _read_json(root / "metadata.json")
    thresholds = {str(k): float(v) for k, v in _read_json(root / "thresholds.json").items()}
    experimental = _load_phase_rows(root / "experimental_tielines.csv")
    feeds = _load_feed_rows(root / "feed_compositions.csv")
    cases: list[BenchmarkCase] = []
    for case_key in sorted(experimental, key=_case_sort_key):
        phases = experimental[case_key]
        feed = feeds[case_key]
        cases.append(
            BenchmarkCase(
                suite=name,
                case_key=case_key,
                tie_line=int(feed["tie_line"]),
                figure=int(feed["figure"]),
                temperature_K=float(feed["temperature_K"]),
                salt_wtfrac=float(feed["salt_wtfrac"]),
                feed_formula=_normalize(np.asarray(feed["formula"], dtype=float)),
                experimental_organic_formula=_normalize(np.asarray(phases["organic"], dtype=float)),
                experimental_aqueous_formula=_normalize(np.asarray(phases["aqueous"], dtype=float)),
                source=str(feed["source"]),
            )
        )
    return BenchmarkSuite(
        name=name,
        species=tuple(metadata["species"]),
        formula_species=tuple(metadata["formula_species"]),
        dataset=str(metadata["parameter_dataset"]),
        thresholds=thresholds,
        cases=tuple(cases),
    )


def benchmark_case_to_native_inputs(case: BenchmarkCase) -> NativeBenchmarkInputs:
    feed = formula_to_explicit(case.feed_formula)
    return NativeBenchmarkInputs(feed=feed)


def run_smoke_cases(suite: BenchmarkSuite, case_keys: Sequence[str]) -> list[BenchmarkPrediction]:
    case_map = {case.case_key: case for case in suite.cases}
    return [predict_case(suite, case_map[key], max_iterations=5) for key in case_keys]


def predict_case(
    suite: BenchmarkSuite,
    case: BenchmarkCase,
    *,
    max_iterations: int = 180,
    tolerance: float = 1.0e-8,
) -> BenchmarkPrediction:
    native = benchmark_case_to_native_inputs(case)
    mixture = ePCSAFTMixture.from_dataset(suite.dataset, suite.species, native.feed, case.temperature_K)
    try:
        result = mixture.equilibrium(
            kind="electrolyte_lle",
            T=case.temperature_K,
            P=PRESSURE_PA,
            z=native.feed,
            options=EquilibriumOptions(
                max_iterations=max_iterations,
                tolerance=tolerance,
            ),
        )
    except SolutionError as exc:
        diagnostics = _diagnostics_from_exception(exc)
        return BenchmarkPrediction(case, "rejected", diagnostics, None, None, None, None)

    diagnostics = dict(result.diagnostics)
    phases = {phase.label: phase for phase in result.phases}
    if "aq" not in phases or "org" not in phases:
        diagnostics["rejection_reason"] = "native result did not contain aq/org phase labels"
        return BenchmarkPrediction(case, "rejected", diagnostics, None, None, None, None)
    aq_formula = explicit_to_formula(phases["aq"].composition)
    org_formula = explicit_to_formula(phases["org"].composition)
    metrics = _metrics(case, org_formula, aq_formula, diagnostics)
    return BenchmarkPrediction(
        case=case,
        status="accepted",
        diagnostics=diagnostics,
        metrics=metrics,
        organic_formula=org_formula,
        aqueous_formula=aq_formula,
        beta_org=float(phases["org"].phase_fraction),
    )


def run_confidence_suite(
    suite: str | BenchmarkSuite = "khudaida_2026",
    *,
    mode: str = "full",
    output_root: str | Path | None = None,
    write_gallery: bool = False,
) -> ConfidenceReport:
    benchmark_suite = load_benchmark_suite(suite) if isinstance(suite, str) else suite
    out_dir = Path(output_root) if output_root is not None else DEFAULT_OUTPUT_ROOT
    out_dir = out_dir / benchmark_suite.name
    out_dir.mkdir(parents=True, exist_ok=True)

    selected_cases = _cases_for_mode(benchmark_suite, mode)
    predictions = [predict_case(benchmark_suite, case, max_iterations=5) for case in selected_cases]
    continuation = continuation_metrics(predictions)
    if mode == "full":
        oracle = oracle_checks(benchmark_suite, predictions[: min(6, len(predictions))])
        stress = stress_cases(benchmark_suite)
    else:
        oracle = oracle_checks(benchmark_suite, predictions[:1])
        stress = []

    report = ConfidenceReport(
        output_dir=out_dir,
        summary_path=out_dir / "summary.json",
        benchmark_csv=out_dir / "benchmark_predictions.csv",
        continuation_csv=out_dir / "continuation_metrics.csv",
        oracle_csv=out_dir / "oracle_checks.csv",
        stress_csv=out_dir / "stress_cases.csv",
        residual_gate_plot=out_dir / "residual_gate_summary.png",
        error_plot=out_dir / "per_species_error_summary.png",
        all_tielines_plot=out_dir / "khudaida_2026_all_tielines.png",
        continuation_plot=out_dir / "continuation_smoothness.png",
    )

    _write_benchmark_csv(report.benchmark_csv, predictions)
    _write_continuation_csv(report.continuation_csv, continuation)
    _write_oracle_csv(report.oracle_csv, oracle)
    _write_stress_csv(report.stress_csv, stress)
    _write_summary(report.summary_path, benchmark_suite, mode, predictions, continuation, oracle, stress)
    _write_plots(report, benchmark_suite, predictions, continuation)
    if write_gallery:
        _write_gallery_copies(report)
    return report


def _cases_for_mode(suite: BenchmarkSuite, mode: str) -> tuple[BenchmarkCase, ...]:
    if mode == "full":
        return suite.cases
    if mode == "smoke":
        return suite.cases[: min(2, len(suite.cases))]
    raise ValueError("mode must be 'smoke' or 'full'.")


def continuation_metrics(predictions: Sequence[BenchmarkPrediction]) -> list[ContinuationMetrics]:
    grouped: dict[str, list[BenchmarkPrediction]] = {}
    for prediction in predictions:
        key = f"{prediction.case.salt_wtfrac:.2f}:{prediction.case.temperature_K:.2f}"
        grouped.setdefault(key, []).append(prediction)
    rows: list[ContinuationMetrics] = []
    for series_key, series in grouped.items():
        ordered = sorted(series, key=lambda item: item.case.tie_line)
        for left, right in zip(ordered, ordered[1:]):
            if left.outcome != "accepted" or right.outcome != "accepted":
                rows.append(
                    ContinuationMetrics(
                        series_key, left.case.case_key, right.case.case_key, math.nan, math.nan, False, True
                    )
                )
                continue
            assert left.organic_formula is not None and left.aqueous_formula is not None
            assert right.organic_formula is not None and right.aqueous_formula is not None
            jump = max(
                float(np.max(np.abs(right.organic_formula - left.organic_formula))),
                float(np.max(np.abs(right.aqueous_formula - left.aqueous_formula))),
            )
            beta_jump = abs(float((right.beta_org or 0.0) - (left.beta_org or 0.0)))
            phase_distance = min(
                float(left.diagnostics.get("phase_distance", math.inf)),
                float(right.diagnostics.get("phase_distance", math.inf)),
            )
            rows.append(
                ContinuationMetrics(
                    series_key=series_key,
                    from_case=left.case.case_key,
                    to_case=right.case.case_key,
                    composition_jump_norm=jump,
                    beta_jump=beta_jump,
                    phase_label_swapped=bool(right.diagnostics.get("phase_labels_swapped", False)),
                    branch_collapse_flag=phase_distance < 1.0e-4 or jump > 0.25 or beta_jump > 0.35,
                )
            )
    return rows


def oracle_checks(suite: BenchmarkSuite, predictions: Sequence[BenchmarkPrediction]) -> list[OracleCheck]:
    checks: list[OracleCheck] = []
    for prediction in predictions:
        case = prediction.case
        native_delta = float(prediction.diagnostics.get("gibbs_delta", math.nan))
        oracle_delta = fixed_phase_gibbs_delta(
            suite, case, case.experimental_organic_formula, case.experimental_aqueous_formula
        )
        checks.append(
            OracleCheck(
                case_key=case.case_key,
                native_outcome=prediction.outcome,
                native_gibbs_delta=native_delta,
                oracle_gibbs_delta=oracle_delta,
                native_minus_oracle_gibbs_delta=(
                    float(native_delta - oracle_delta) if np.isfinite(native_delta) else math.nan
                ),
                oracle_flags_native=bool(np.isfinite(native_delta) and native_delta - oracle_delta > 1.0e-6),
            )
        )
    return checks


def fixed_phase_gibbs_delta(
    suite: BenchmarkSuite,
    case: BenchmarkCase,
    organic_formula: np.ndarray,
    aqueous_formula: np.ndarray,
) -> float:
    native = benchmark_case_to_native_inputs(case)
    org = formula_to_explicit(organic_formula)
    aq = formula_to_explicit(aqueous_formula)
    beta = _best_phase_fraction(native.feed, org, aq)
    mixture = ePCSAFTMixture.from_dataset(suite.dataset, suite.species, native.feed, case.temperature_K)
    feed_state = phase_state(mixture, case.temperature_K, PRESSURE_PA, native.feed, "liq", "oracle_feed")[
        "state"
    ]
    org_state = phase_state(mixture, case.temperature_K, PRESSURE_PA, org, "liq", "oracle_org")["state"]
    aq_state = phase_state(mixture, case.temperature_K, PRESSURE_PA, aq, "liq", "oracle_aq")["state"]
    g_feed = _gibbs_proxy(native.feed, feed_state.fugacity_coefficient())
    g_org = _gibbs_proxy(org, org_state.fugacity_coefficient())
    g_aq = _gibbs_proxy(aq, aq_state.fugacity_coefficient())
    return float(beta * g_org + (1.0 - beta) * g_aq - g_feed)


def stress_cases(suite: BenchmarkSuite) -> list[BenchmarkPrediction]:
    base = _preferred_case(suite, "0.10:303.15:1")
    formulas = {
        "very_dilute_salt": _with_salt(base.feed_formula, 1.0e-6),
        "high_salt": _with_salt(base.feed_formula, 0.12),
        "trace_organic": _normalize(np.asarray([0.86, 0.001, 0.05, 0.089], dtype=float)),
        "near_plait_like": _normalize(0.5 * (base.experimental_organic_formula + base.experimental_aqueous_formula)),
        "nearly_identical_phases": _normalize(base.feed_formula + np.asarray([1.0e-5, -5.0e-6, -5.0e-6, 0.0])),
    }
    predictions: list[BenchmarkPrediction] = []
    for name, feed in formulas.items():
        case = BenchmarkCase(
            suite=base.suite,
            case_key=f"stress:{name}",
            tie_line=base.tie_line,
            figure=base.figure,
            temperature_K=base.temperature_K,
            salt_wtfrac=base.salt_wtfrac,
            feed_formula=feed,
            experimental_organic_formula=base.experimental_organic_formula,
            experimental_aqueous_formula=base.experimental_aqueous_formula,
            source="synthetic_stress",
        )
        predictions.append(predict_case(suite, case, max_iterations=5, tolerance=1.0e-8))
    return predictions


def formula_to_explicit(formula: Sequence[float]) -> np.ndarray:
    x = _normalize(np.asarray(formula, dtype=float))
    expanded = np.asarray([x[0], x[1], x[2], x[3], x[3]], dtype=float)
    return _normalize(expanded)


def explicit_to_formula(explicit: Sequence[float]) -> np.ndarray:
    x = _normalize(np.asarray(explicit, dtype=float))
    salt = 0.5 * (x[3] + x[4])
    return _normalize(np.asarray([x[0], x[1], x[2], salt], dtype=float))


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Run opt-in electrolyte LLE confidence validation.")
    parser.add_argument("--suite", default="khudaida_2026")
    parser.add_argument("--mode", choices=("smoke", "full"), default="full")
    parser.add_argument("--output-root", default=str(DEFAULT_OUTPUT_ROOT))
    parser.add_argument("--write-gallery", action="store_true")
    args = parser.parse_args(argv)
    report = run_confidence_suite(
        args.suite, mode=args.mode, output_root=args.output_root, write_gallery=args.write_gallery
    )
    print(f"Wrote electrolyte LLE confidence report: {report.output_dir}")
    return 0


def _read_json(path: Path) -> dict[str, Any]:
    with path.open(encoding="utf-8") as handle:
        return json.load(handle)


def _load_phase_rows(path: Path) -> dict[str, dict[str, np.ndarray]]:
    grouped: dict[str, dict[str, np.ndarray]] = {}
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            grouped.setdefault(row["case_key"], {})[row["phase"]] = np.asarray(
                [float(row["x_water"]), float(row["x_ethanol"]), float(row["x_isobutanol"]), float(row["x_nacl"])],
                dtype=float,
            )
    return grouped


def _load_feed_rows(path: Path) -> dict[str, dict[str, Any]]:
    rows: dict[str, dict[str, Any]] = {}
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            rows[row["case_key"]] = {
                "tie_line": int(row["tie_line"]),
                "figure": int(row["figure"]),
                "temperature_K": float(row["temperature_K"]),
                "salt_wtfrac": float(row["salt_wtfrac"]),
                "formula": np.asarray(
                    [
                        float(row["x_water_total"]),
                        float(row["x_ethanol_total"]),
                        float(row["x_isobutanol_total"]),
                        float(row["x_nacl_total"]),
                    ],
                    dtype=float,
                ),
                "source": row["source"],
            }
    return rows


def _case_sort_key(case_key: str) -> tuple[float, float, int]:
    salt, temperature, tie_line = case_key.split(":")
    return float(salt), float(temperature), int(tie_line)


def _normalize(values: np.ndarray) -> np.ndarray:
    total = float(np.sum(values))
    if not np.isfinite(total) or total <= 0.0:
        raise ValueError("composition total must be positive")
    clipped = np.clip(np.asarray(values, dtype=float), 0.0, None)
    return clipped / float(np.sum(clipped))


def _best_phase_fraction(feed: np.ndarray, organic: np.ndarray, aqueous: np.ndarray) -> float:
    direction = organic - aqueous
    denom = float(np.dot(direction, direction))
    if denom <= 0.0:
        return 0.5
    return float(np.clip(np.dot(feed - aqueous, direction) / denom, 1.0e-8, 1.0 - 1.0e-8))


def _metrics(
    case: BenchmarkCase, organic: np.ndarray, aqueous: np.ndarray, diagnostics: Mapping[str, Any]
) -> BenchmarkMetrics:
    organic_delta = np.abs(organic - case.experimental_organic_formula)
    aqueous_delta = np.abs(aqueous - case.experimental_aqueous_formula)
    return BenchmarkMetrics(
        grand_aad=float((np.sum(organic_delta) + np.sum(aqueous_delta)) / 8.0),
        max_abs_error=float(max(np.max(organic_delta), np.max(aqueous_delta))),
        organic_aad=float(np.mean(organic_delta)),
        aqueous_aad=float(np.mean(aqueous_delta)),
        phase_distance=float(diagnostics.get("phase_distance", np.max(np.abs(organic - aqueous)))),
    )


def _diagnostics_from_exception(exc: SolutionError) -> dict[str, Any]:
    diagnostics = getattr(exc, "diagnostics", None)
    if diagnostics is None and len(exc.args) > 1 and isinstance(exc.args[1], dict):
        diagnostics = exc.args[1]
    return dict(diagnostics or {"message": str(exc), "acceptance_gate": "predictive_solve_failed"})


def _gibbs_proxy(composition: np.ndarray, ln_phi: Sequence[float]) -> float:
    comp = np.clip(np.asarray(composition, dtype=float), 1.0e-300, 1.0)
    return float(np.sum(comp * (np.log(comp) + np.asarray(ln_phi, dtype=float))))


def _preferred_case(suite: BenchmarkSuite, case_key: str) -> BenchmarkCase:
    for case in suite.cases:
        if case.case_key == case_key:
            return case
    return suite.cases[0]


def _with_salt(feed_formula: np.ndarray, salt: float) -> np.ndarray:
    neutrals = _normalize(np.asarray(feed_formula[:3], dtype=float))
    return _normalize(np.asarray([*(neutrals * (1.0 - salt)), salt], dtype=float))


def _write_benchmark_csv(path: Path, predictions: Sequence[BenchmarkPrediction]) -> None:
    fieldnames = [
        "case_key",
        "outcome",
        "temperature_K",
        "salt_wtfrac",
        "tie_line",
        "grand_aad",
        "max_abs_error",
        "organic_aad",
        "aqueous_aad",
        "beta_org",
        "solver_residual_norm",
        "material_balance_error",
        "charge_balance_error",
        "gibbs_delta",
        "phase_distance",
        "tpd_trial_count",
        "tpd_polish_iterations",
        "acceptance_gate",
        "rejection_reason",
    ]
    rows = []
    for prediction in predictions:
        diag = prediction.diagnostics
        metrics = prediction.metrics
        rows.append(
            {
                "case_key": prediction.case.case_key,
                "outcome": prediction.outcome,
                "temperature_K": prediction.case.temperature_K,
                "salt_wtfrac": prediction.case.salt_wtfrac,
                "tie_line": prediction.case.tie_line,
                "grand_aad": metrics.grand_aad if metrics else math.nan,
                "max_abs_error": metrics.max_abs_error if metrics else math.nan,
                "organic_aad": metrics.organic_aad if metrics else math.nan,
                "aqueous_aad": metrics.aqueous_aad if metrics else math.nan,
                "beta_org": prediction.beta_org if prediction.beta_org is not None else math.nan,
                "solver_residual_norm": diag.get("solver_residual_norm", math.nan),
                "material_balance_error": diag.get("material_balance_error", math.nan),
                "charge_balance_error": diag.get("charge_balance_error", math.nan),
                "gibbs_delta": diag.get("gibbs_delta", math.nan),
                "phase_distance": diag.get("phase_distance", math.nan),
                "tpd_trial_count": diag.get("tpd_trial_count", math.nan),
                "tpd_polish_iterations": diag.get("tpd_polish_iterations", math.nan),
                "acceptance_gate": diag.get("acceptance_gate", ""),
                "rejection_reason": diag.get("rejection_reason", ""),
            }
        )
    _write_csv(path, fieldnames, rows)


def _write_continuation_csv(path: Path, rows: Sequence[ContinuationMetrics]) -> None:
    _write_csv(path, list(ContinuationMetrics.__dataclass_fields__), [row.__dict__ for row in rows])


def _write_oracle_csv(path: Path, rows: Sequence[OracleCheck]) -> None:
    _write_csv(path, list(OracleCheck.__dataclass_fields__), [row.__dict__ for row in rows])


def _write_stress_csv(path: Path, rows: Sequence[BenchmarkPrediction]) -> None:
    _write_csv(
        path,
        [
            "case_key",
            "outcome",
            "acceptance_gate",
            "rejection_reason",
            "solver_residual_norm",
            "gibbs_delta",
            "phase_distance",
        ],
        [
            {
                "case_key": row.case.case_key,
                "outcome": row.outcome,
                "acceptance_gate": row.diagnostics.get("acceptance_gate", ""),
                "rejection_reason": row.diagnostics.get("rejection_reason", ""),
                "solver_residual_norm": row.diagnostics.get("solver_residual_norm", math.nan),
                "gibbs_delta": row.diagnostics.get("gibbs_delta", math.nan),
                "phase_distance": row.diagnostics.get("phase_distance", math.nan),
            }
            for row in rows
        ],
    )


def _write_summary(
    path: Path,
    suite: BenchmarkSuite,
    mode: str,
    predictions: Sequence[BenchmarkPrediction],
    continuation: Sequence[ContinuationMetrics],
    oracle: Sequence[OracleCheck],
    stress: Sequence[BenchmarkPrediction],
) -> None:
    accepted = [item for item in predictions if item.outcome == "accepted"]
    summary = {
        "suite": suite.name,
        "mode": mode,
        "case_count": len(predictions),
        "accepted_count": len(accepted),
        "rejected_count": len(predictions) - len(accepted),
        "max_grand_aad": max((item.metrics.grand_aad for item in accepted if item.metrics), default=math.nan),
        "max_abs_error": max((item.metrics.max_abs_error for item in accepted if item.metrics), default=math.nan),
        "continuation_rows": len(continuation),
        "continuation_branch_flags": sum(1 for row in continuation if row.branch_collapse_flag),
        "oracle_rows": len(oracle),
        "oracle_flags": sum(1 for row in oracle if row.oracle_flags_native),
        "stress_rows": len(stress),
    }
    path.write_text(json.dumps(_json_ready(summary), indent=2, allow_nan=False), encoding="utf-8")


def _write_plots(
    report: ConfidenceReport,
    suite: BenchmarkSuite,
    predictions: Sequence[BenchmarkPrediction],
    continuation: Sequence[ContinuationMetrics],
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    accepted = [item for item in predictions if item.outcome == "accepted"]
    residual_labels = ["solver_residual_norm", "material_balance_error", "charge_balance_error", "phase_distance"]
    residual_values = [
        max(
            (
                float(item.diagnostics.get(label, math.nan))
                for item in predictions
                if np.isfinite(float(item.diagnostics.get(label, math.nan)))
            ),
            default=math.nan,
        )
        for label in residual_labels
    ]
    threshold_values = [suite.thresholds.get(label, math.nan) for label in residual_labels]
    _bar_comparison_plot(
        report.residual_gate_plot,
        "Electrolyte LLE residual gate summary",
        residual_labels,
        residual_values,
        threshold_values,
    )

    error_labels = ["organic AAD", "aqueous AAD", "grand AAD", "max abs error"]
    error_values = [
        max((item.metrics.organic_aad for item in accepted if item.metrics), default=math.nan),
        max((item.metrics.aqueous_aad for item in accepted if item.metrics), default=math.nan),
        max((item.metrics.grand_aad for item in accepted if item.metrics), default=math.nan),
        max((item.metrics.max_abs_error for item in accepted if item.metrics), default=math.nan),
    ]
    _bar_plot(
        report.error_plot, "Khudaida benchmark composition errors", error_labels, error_values, "Mole fraction error"
    )
    _all_tielines_plot(report.all_tielines_plot, predictions)

    cont_labels = [row.to_case for row in continuation]
    cont_values = [row.composition_jump_norm if np.isfinite(row.composition_jump_norm) else 0.0 for row in continuation]
    _line_plot(
        report.continuation_plot,
        "Khudaida continuation smoothness",
        cont_labels,
        cont_values,
        "Adjacent composition jump",
    )

    plt.close("all")


def _bar_comparison_plot(
    path: Path, title: str, labels: Sequence[str], actual: Sequence[float], expected: Sequence[float]
) -> None:
    import matplotlib.pyplot as plt

    x = np.arange(len(labels))
    fig, ax = plt.subplots(figsize=(8.0, 4.6))
    ax.bar(x - 0.18, actual, width=0.36, label="Actual max")
    ax.bar(x + 0.18, expected, width=0.36, label="Report threshold")
    ax.set_yscale("log")
    ax.set_xticks(x, labels, rotation=30, ha="right")
    ax.set_title(title)
    ax.legend()
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    fig.savefig(path.with_suffix(".svg"))
    plt.close(fig)
    _write_plot_data(path, [("actual", labels, actual), ("threshold", labels, expected)])


def _bar_plot(path: Path, title: str, labels: Sequence[str], values: Sequence[float], ylabel: str) -> None:
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(8.0, 4.6))
    ax.bar(np.arange(len(labels)), values)
    ax.set_xticks(np.arange(len(labels)), labels, rotation=30, ha="right")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    fig.savefig(path.with_suffix(".svg"))
    plt.close(fig)
    _write_plot_data(path, [("value", labels, values)])


def _line_plot(path: Path, title: str, labels: Sequence[str], values: Sequence[float], ylabel: str) -> None:
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(9.0, 4.6))
    ax.plot(np.arange(len(labels)), values, marker="o")
    step = max(1, len(labels) // 8)
    ax.set_xticks(
        np.arange(0, len(labels), step), [labels[i] for i in range(0, len(labels), step)], rotation=30, ha="right"
    )
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    fig.savefig(path.with_suffix(".svg"))
    plt.close(fig)
    _write_plot_data(path, [("value", labels, values)])


def _all_tielines_plot(path: Path, predictions: Sequence[BenchmarkPrediction]) -> None:
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(8.4, 6.0))
    rows: list[dict[str, Any]] = []
    outcome_colors = {"accepted": "#2ca02c", "rejected": "#d62728"}
    grouped: dict[int, list[BenchmarkPrediction]] = {}
    for prediction in predictions:
        grouped.setdefault(prediction.case.figure, []).append(prediction)

    for figure, group in sorted(grouped.items()):
        color = plt.cm.tab10((figure - 1) % 10)
        for prediction in sorted(group, key=lambda item: item.case.tie_line):
            case = prediction.case
            exp_org = case.experimental_organic_formula
            exp_aq = case.experimental_aqueous_formula
            label = (
                f"Fig. {figure} experimental" if case.tie_line == min(item.case.tie_line for item in group) else None
            )
            ax.plot(
                [exp_aq[1], exp_org[1]],
                [exp_aq[2], exp_org[2]],
                color=color,
                alpha=0.55,
                linewidth=1.2,
                marker="o",
                markersize=3.0,
                label=label,
            )
            rows.extend(
                [
                    _tieline_plot_row(case, "experimental_aqueous", prediction.outcome, exp_aq),
                    _tieline_plot_row(case, "experimental_organic", prediction.outcome, exp_org),
                    _tieline_plot_row(case, "feed", prediction.outcome, case.feed_formula),
                ]
            )
            ax.scatter(
                [case.feed_formula[1]],
                [case.feed_formula[2]],
                color=outcome_colors.get(prediction.outcome, "#7f7f7f"),
                marker="x",
                s=24,
                linewidths=1.2,
                alpha=0.8,
            )
            if (
                prediction.outcome == "accepted"
                and prediction.organic_formula is not None
                and prediction.aqueous_formula is not None
            ):
                ax.plot(
                    [prediction.aqueous_formula[1], prediction.organic_formula[1]],
                    [prediction.aqueous_formula[2], prediction.organic_formula[2]],
                    color="#111111",
                    linestyle="--",
                    linewidth=1.0,
                    alpha=0.75,
                )
                rows.extend(
                    [
                        _tieline_plot_row(case, "native_aqueous", prediction.outcome, prediction.aqueous_formula),
                        _tieline_plot_row(case, "native_organic", prediction.outcome, prediction.organic_formula),
                    ]
                )

    ax.scatter([], [], color=outcome_colors["accepted"], marker="x", label="native accepted feed")
    ax.scatter([], [], color=outcome_colors["rejected"], marker="x", label="native rejected feed")
    ax.set_xlabel(r"$x_{\mathrm{ethanol}}$")
    ax.set_ylabel(r"$x_{\mathrm{isobutanol}}$")
    ax.set_title("Khudaida 2026 electrolyte LLE tie-lines and native result outcome")
    ax.grid(True, alpha=0.25)
    ax.legend(loc="best", fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    fig.savefig(path.with_suffix(".svg"))
    plt.close(fig)
    _write_csv(
        path.parent / f"{path.stem}.csv",
        ["case_key", "figure", "tie_line", "point_type", "outcome", "x_water", "x_ethanol", "x_isobutanol", "x_nacl"],
        rows,
    )
    _write_mpl_style_contract(path, "Khudaida 2026 electrolyte LLE tie-lines and native result outcome")


def _tieline_plot_row(
    case: BenchmarkCase, point_type: str, outcome: str, composition: Sequence[float]
) -> dict[str, Any]:
    values = list(composition)
    return {
        "case_key": case.case_key,
        "figure": case.figure,
        "tie_line": case.tie_line,
        "point_type": point_type,
        "outcome": outcome,
        "x_water": values[0],
        "x_ethanol": values[1],
        "x_isobutanol": values[2],
        "x_nacl": values[3],
    }


def _write_plot_data(path: Path, series_rows: Sequence[tuple[str, Sequence[str], Sequence[float]]]) -> None:
    rows = []
    for series, labels, values in series_rows:
        for label, value in zip(labels, values):
            rows.append({"plot": path.name, "series": series, "label": label, "value": value})
    _write_csv(path.parent / f"{path.stem}.csv", ["plot", "series", "label", "value"], rows)
    _write_mpl_style_contract(path, "")


def _write_mpl_style_contract(path: Path, title: str) -> None:
    style_path = path.parent / f"{path.stem}.mpl.yaml"
    if style_path.exists():
        return
    style_path.write_text(
        "\n".join(
            [
                "# Matplotlib plot-set style contract.",
                "# Edit this sidecar, then rerun the owning confidence workflow.",
                "figure:",
                f'  file: "{path.name}"',
                f'  format: "{path.suffix.lstrip(".")}"',
                "axes:",
                "  - index: 0",
                f'    title: "{title}"',
                "",
            ]
        ),
        encoding="utf-8",
    )


def _write_gallery_copies(report: ConfidenceReport) -> None:
    # External gallery tools should discover analysis-owned results directly.
    return None


def _write_csv(path: Path, fieldnames: Sequence[str], rows: Iterable[Mapping[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    frame = pd.DataFrame(
        [{key: _csv_value(row.get(key, "")) for key in fieldnames} for row in rows], columns=list(fieldnames)
    )
    frame.to_csv(path, index=False)


def _csv_value(value: Any) -> Any:
    if isinstance(value, float):
        if not np.isfinite(value):
            return ""
        return format(value, CSV_FLOAT_FORMAT)
    return value


def _json_ready(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): _json_ready(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_ready(item) for item in value]
    if isinstance(value, np.ndarray):
        return _json_ready(value.tolist())
    if isinstance(value, np.generic):
        return _json_ready(value.item())
    if isinstance(value, float) and not np.isfinite(value):
        return None
    return value
