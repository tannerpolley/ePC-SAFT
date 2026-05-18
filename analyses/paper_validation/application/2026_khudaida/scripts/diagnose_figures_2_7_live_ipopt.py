from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from collections.abc import Mapping
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np

for _candidate in Path(__file__).resolve().parents:
    if (_candidate / "src" / "epcsaft").is_dir():
        if str(_candidate / "src") not in sys.path:
            sys.path.insert(0, str(_candidate / "src"))
        break
else:
    raise ModuleNotFoundError("Could not locate repo root containing src/epcsaft")

import epcsaft
from epcsaft import _core
from epcsaft.equilibrium import _native_timeout_seconds, neutral_two_phase_eos_tolerances


FIGURE_CASES: dict[int, tuple[float, float]] = {
    2: (293.15, 0.05),
    3: (303.15, 0.05),
    4: (313.15, 0.05),
    5: (293.15, 0.10),
    6: (303.15, 0.10),
    7: (313.15, 0.10),
}
ANALYSIS_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_OUTPUT_DIR = ANALYSIS_ROOT / "diagnostics" / "figures_2_7_live_ipopt"
P_REF = 1.0e5
PARAMETER_DATASET = "2026_Khudaida"
SPECIES = ["H2O", "Ethanol", "Butanol", "Na+", "Cl-"]
MW_FORMULA = np.asarray([18.01528e-3, 46.068e-3, 74.1216e-3, 58.43e-3], dtype=float)
RESIDUAL_LABELS = [
    "phase_equilibrium:H2O",
    "phase_equilibrium:Ethanol",
    "phase_equilibrium:Butanol",
    "phase_equilibrium:NaCl_ion_pair",
    "material_balance:H2O",
    "material_balance:Ethanol",
    "material_balance:Butanol",
    "material_balance:Na+",
    "material_balance:Cl-",
]


def _jsonable(value: Any) -> Any:
    if isinstance(value, np.ndarray):
        return _jsonable(value.tolist())
    if isinstance(value, np.generic):
        return _jsonable(value.item())
    if isinstance(value, float):
        return value if math.isfinite(value) else None
    if isinstance(value, Mapping):
        return {str(key): _jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(item) for item in value]
    return value


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _format_vector(value: Any) -> str:
    if value is None:
        return ""
    array = np.asarray(value, dtype=float).reshape(-1)
    if array.size == 0:
        return ""
    return json.dumps(_jsonable(array), separators=(",", ":"))


def _safe_float(value: Any) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return float("nan")


def figure_root(figure_number: int) -> Path:
    return ANALYSIS_ROOT / "figures" / f"figure_{figure_number}"


def _source_feed_path(figure_number: int) -> Path:
    return figure_root(figure_number) / "input" / "feed_compositions_digitized.csv"


def _model_cache_path(figure_number: int) -> Path:
    return figure_root(figure_number) / "output" / "data" / "model_tielines.csv"


def _raw_digitized_feed_rows(figure_number: int) -> list[dict[str, str]]:
    path = _source_feed_path(figure_number)
    with path.open("r", newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def _build_feed_formula_from_salt_free_molefractions(salt_free_xyz: np.ndarray, salt_wt: float) -> np.ndarray:
    salt_free_xyz = np.asarray(salt_free_xyz, dtype=float)
    salt_free_xyz = salt_free_xyz / np.sum(salt_free_xyz)
    neutral_average_mw = float(np.dot(salt_free_xyz, MW_FORMULA[:3]))
    salt_moles_per_neutral_mole = salt_wt * neutral_average_mw / ((1.0 - salt_wt) * MW_FORMULA[3])
    total_moles = 1.0 + salt_moles_per_neutral_mole
    return np.asarray(
        [
            salt_free_xyz[0] / total_moles,
            salt_free_xyz[1] / total_moles,
            salt_free_xyz[2] / total_moles,
            salt_moles_per_neutral_mole / total_moles,
        ],
        dtype=float,
    )


def formula_to_ion_basis(x_formula: np.ndarray) -> np.ndarray:
    x_formula = np.asarray(x_formula, dtype=float)
    denom = 1.0 + float(x_formula[3])
    return np.asarray(
        [
            x_formula[0] / denom,
            x_formula[1] / denom,
            x_formula[2] / denom,
            x_formula[3] / denom,
            x_formula[3] / denom,
        ],
        dtype=float,
    )


def _selected_feed_rows(
    figure_number: int,
    *,
    all_rows: bool,
    max_rows_per_figure: int,
    tie_lines: set[int] | None,
) -> list[tuple[dict[str, Any], dict[str, str]]]:
    temperature_k, salt_wt = FIGURE_CASES[figure_number]
    if not _source_feed_path(figure_number).exists():
        raise FileNotFoundError(f"Missing digitized feed CSV for Figure {figure_number}: {_source_feed_path(figure_number)}")
    raw_rows = _raw_digitized_feed_rows(figure_number)
    feed_rows = []
    for idx, raw in enumerate(raw_rows, start=1):
        ethanol = float(raw["x_ethanol_salt_free"])
        isobutanol = float(raw["x_isobutanol_salt_free"])
        water = float(raw["x_water_salt_free"])
        feed_formula = _build_feed_formula_from_salt_free_molefractions(
            np.asarray([water, ethanol, isobutanol], dtype=float),
            salt_wt,
        )
        feed_rows.append(
            {
                "tie_line": idx,
                "temperature_K": float(temperature_k),
                "salt_wtfrac": float(salt_wt),
                "feed_formula": feed_formula,
                "source": raw.get("source", "digitized_user_supplied") or "digitized_user_supplied",
            }
        )
    pairs = [(feed, raw) for feed, raw in zip(feed_rows, raw_rows, strict=True)]
    if tie_lines is not None:
        pairs = [(feed, raw) for feed, raw in pairs if int(feed["tie_line"]) in tie_lines]
    if not all_rows:
        pairs = pairs[: max(0, max_rows_per_figure)]
    return pairs


def _parse_float(value: str | float | int | None) -> float:
    if value in (None, ""):
        return float("nan")
    return float(value)


def _load_model_rows_from_csv(path: Path) -> list[dict[str, Any]]:
    grouped: dict[int, dict[str, Any]] = {}
    with path.open("r", newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        for raw in reader:
            tie_line = int(raw["tie_line"])
            entry = grouped.setdefault(
                tie_line,
                {
                    "tie_line": tie_line,
                    "temperature_K": _parse_float(raw.get("temperature_K")),
                    "salt_wtfrac": _parse_float(raw.get("salt_wtfrac")),
                    "organic_formula": np.full(4, np.nan),
                    "aqueous_formula": np.full(4, np.nan),
                    "beta_organic": np.nan,
                    "beta_aqueous": np.nan,
                    "residual_norm": np.nan,
                    "objective": np.nan,
                    "converged": False,
                },
            )
            phase_vec = np.asarray(
                [
                    _parse_float(raw.get("x_water")),
                    _parse_float(raw.get("x_ethanol")),
                    _parse_float(raw.get("x_isobutanol")),
                    _parse_float(raw.get("x_nacl")),
                ],
                dtype=float,
            )
            phase_name = str(raw.get("phase", "")).strip().lower()
            if phase_name == "organic":
                entry["organic_formula"] = phase_vec
                entry["beta_organic"] = _parse_float(raw.get("beta"))
            elif phase_name == "aqueous":
                entry["aqueous_formula"] = phase_vec
                entry["beta_aqueous"] = _parse_float(raw.get("beta"))
            entry["residual_norm"] = _parse_float(raw.get("residual_norm"))
            entry["objective"] = _parse_float(raw.get("objective"))
            converged_raw = str(raw.get("converged", "")).strip().lower()
            entry["converged"] = converged_raw in {"true", "1", "yes"}
    return [grouped[key] for key in sorted(grouped)]


def _cached_model_row(figure_number: int, tie_line: int) -> dict[str, Any] | None:
    cache_path = _model_cache_path(figure_number)
    if not cache_path.exists():
        return None
    for row in _load_model_rows_from_csv(cache_path):
        if int(row["tie_line"]) == int(tie_line):
            return row
    return None


def _native_route_payload(
    *,
    mixture: epcsaft.ePCSAFTMixture,
    temperature_k: float,
    z_feed: np.ndarray,
    max_iterations: int,
    tolerance: float,
    timeout_seconds: float | None,
) -> dict[str, Any]:
    options = epcsaft.EquilibriumOptions(
        max_iterations=max_iterations,
        tolerance=tolerance,
        timeout_seconds=timeout_seconds,
    )
    material_tolerance, pressure_tolerance, chemical_potential_tolerance, phase_distance_tolerance = (
        neutral_two_phase_eos_tolerances(P_REF, options)
    )
    return _core._native_electrolyte_lle_eos_route_result(
        mixture._native,
        float(temperature_k),
        P_REF,
        z_feed.tolist(),
        options.max_iterations,
        options.tolerance,
        _native_timeout_seconds(options),
        material_tolerance,
        pressure_tolerance,
        min(options.tolerance, 1.0e-8),
        chemical_potential_tolerance,
        phase_distance_tolerance,
    )


def _residual_at_route_variables(
    *,
    mixture: epcsaft.ePCSAFTMixture,
    temperature_k: float,
    z_feed: np.ndarray,
    variables: Any,
    max_iterations: int,
    tolerance: float,
) -> dict[str, Any] | None:
    variables_array = np.asarray(variables, dtype=float).reshape(-1)
    if variables_array.size == 0 or not np.all(np.isfinite(variables_array)):
        return None
    request = {
        "T": float(temperature_k),
        "P": P_REF,
        "z": z_feed.tolist(),
        "species": SPECIES,
        "max_iterations": max_iterations,
        "tolerance": tolerance,
        "variables": variables_array.tolist(),
    }
    return _core._evaluate_electrolyte_lle_residual_native(mixture._native, request)


def _diagnose_feed(
    *,
    figure_number: int,
    feed_row: dict[str, Any],
    raw_feed_row: dict[str, str],
    max_iterations: int,
    tolerance: float,
    timeout_seconds: float | None,
) -> tuple[dict[str, Any], dict[str, Any]]:
    temperature_k, salt_wt = FIGURE_CASES[figure_number]
    feed_formula = np.asarray(feed_row["feed_formula"], dtype=float)
    z_feed = formula_to_ion_basis(feed_formula)
    mixture = epcsaft.ePCSAFTMixture.from_dataset(PARAMETER_DATASET, SPECIES, z_feed, temperature_k)
    route = _native_route_payload(
        mixture=mixture,
        temperature_k=temperature_k,
        z_feed=z_feed,
        max_iterations=max_iterations,
        tolerance=tolerance,
        timeout_seconds=timeout_seconds,
    )
    residual_eval = _residual_at_route_variables(
        mixture=mixture,
        temperature_k=temperature_k,
        z_feed=z_feed,
        variables=route.get("variables", ()),
        max_iterations=max_iterations,
        tolerance=tolerance,
    )
    cached = _cached_model_row(figure_number, int(feed_row["tie_line"]))
    residual = [] if residual_eval is None else list(residual_eval.get("residual", ()))
    residual_array = np.asarray(residual, dtype=float).reshape(-1)
    if residual_array.size:
        residual_max_index = int(np.argmax(np.abs(residual_array)))
        residual_max_value = float(residual_array[residual_max_index])
        residual_inf = float(np.max(np.abs(residual_array)))
    else:
        residual_max_index = -1
        residual_max_value = float("nan")
        residual_inf = float("nan")
    residual_diagnostics = {} if residual_eval is None else dict(residual_eval.get("diagnostics", {}))
    cached_org = None if cached is None else cached.get("organic_formula")
    cached_aq = None if cached is None else cached.get("aqueous_formula")
    summary = {
        "figure": figure_number,
        "tie_line": int(feed_row["tie_line"]),
        "temperature_K": temperature_k,
        "salt_wtfrac": salt_wt,
        "feed_source_path": str(_source_feed_path(figure_number).relative_to(ANALYSIS_ROOT)),
        "feed_source": feed_row.get("source", ""),
        "raw_x_water_salt_free": raw_feed_row.get("x_water_salt_free", ""),
        "raw_x_ethanol_salt_free": raw_feed_row.get("x_ethanol_salt_free", ""),
        "raw_x_isobutanol_salt_free": raw_feed_row.get("x_isobutanol_salt_free", ""),
        "feed_formula": _format_vector(feed_formula),
        "feed_ion_basis": _format_vector(z_feed),
        "route_accepted": bool(route.get("accepted", False)),
        "route_status": route.get("status", ""),
        "solver_status": route.get("solver_status", ""),
        "application_status": route.get("application_status", ""),
        "route_objective": _safe_float(route.get("objective")),
        "route_variable_count": int(np.asarray(route.get("variables", ()), dtype=float).reshape(-1).size),
        "route_variables": _format_vector(route.get("variables", ())),
        "route_phase_amount_rows": int(np.asarray(route.get("phase_amounts", ()), dtype=float).shape[0])
        if len(route.get("phase_amounts", ())) else 0,
        "residual_inf_norm_at_route_variables": residual_inf,
        "residual_max_index": residual_max_index,
        "residual_max_label": RESIDUAL_LABELS[residual_max_index]
        if 0 <= residual_max_index < len(RESIDUAL_LABELS)
        else "",
        "residual_max_value": residual_max_value,
        "phase_fraction_org_at_route_variables": _safe_float(
            None if residual_eval is None else residual_eval.get("phase_fraction_org")
        ),
        "phase_distance_at_route_variables": _safe_float(
            None if residual_eval is None else residual_eval.get("phase_distance")
        ),
        "gibbs_delta_at_route_variables": _safe_float(None if residual_eval is None else residual_eval.get("gibbs_delta")),
        "material_balance_error_at_route_variables": _safe_float(
            None if residual_eval is None else residual_eval.get("material_balance_error")
        ),
        "charge_balance_error_at_route_variables": _safe_float(
            None if residual_eval is None else residual_eval.get("charge_balance_error")
        ),
        "aq_formula_at_route_variables": _format_vector(residual_diagnostics.get("aq_formula")),
        "org_formula_at_route_variables": _format_vector(residual_diagnostics.get("org_formula")),
        "model_cache_path_read_only": str(_model_cache_path(figure_number).relative_to(ANALYSIS_ROOT)),
        "cached_model_converged": "" if cached is None else bool(cached.get("converged", False)),
        "cached_model_residual_norm": "" if cached is None else _safe_float(cached.get("residual_norm")),
        "cached_model_objective": "" if cached is None else _safe_float(cached.get("objective")),
        "cached_organic_formula": _format_vector(cached_org),
        "cached_aqueous_formula": _format_vector(cached_aq),
    }
    detail = {
        "figure": figure_number,
        "tie_line": int(feed_row["tie_line"]),
        "feed_source_path": str(_source_feed_path(figure_number)),
        "raw_digitized_feed_row": raw_feed_row,
        "feed_formula": feed_formula,
        "feed_ion_basis": z_feed,
        "route": route,
        "residual_at_route_variables": residual_eval,
        "cached_model_row": cached,
    }
    return summary, detail


def _write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(rows[0])
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def _write_json(path: Path, payload: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(_jsonable(payload), indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _write_jsonl(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        for row in rows:
            handle.write(json.dumps(_jsonable(row), sort_keys=True) + "\n")


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Diagnose live native-Ipopt recomputation for Khudaida Figures 2-7 using the "
            "figure-owned digitized feed CSVs without modifying existing figure caches."
        )
    )
    parser.add_argument("--figures", nargs="+", type=int, default=sorted(FIGURE_CASES))
    parser.add_argument("--tie-line", action="append", type=int, dest="tie_lines")
    parser.add_argument("--all-rows", action="store_true")
    parser.add_argument("--max-rows-per-figure", type=int, default=1)
    parser.add_argument("--max-iterations", type=int, default=180)
    parser.add_argument("--tolerance", type=float, default=1.0e-8)
    parser.add_argument("--timeout-seconds", type=float, default=60.0)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(sys.argv[1:] if argv is None else argv)
    requested_figures = [int(item) for item in args.figures]
    unknown = sorted(set(requested_figures) - set(FIGURE_CASES))
    if unknown:
        raise ValueError(f"Only Figures 2-7 are supported; got {unknown}.")
    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    tie_lines = None if args.tie_lines is None else set(args.tie_lines)
    timeout_seconds = None if args.timeout_seconds <= 0.0 else float(args.timeout_seconds)
    metadata_path = output_dir / "metadata.json"
    active_case_path = output_dir / "active_case.json"
    for stale_path in (output_dir / "summary.csv", output_dir / "route_payloads.jsonl"):
        stale_path.unlink(missing_ok=True)

    summary_rows: list[dict[str, Any]] = []
    detail_rows: list[dict[str, Any]] = []
    metadata = {
        "purpose": "diagnostic-only live native Ipopt recompute; does not update figure caches",
        "figures": requested_figures,
        "tie_lines": sorted(tie_lines) if tie_lines is not None else None,
        "all_rows": bool(args.all_rows),
        "max_rows_per_figure": int(args.max_rows_per_figure),
        "max_iterations": int(args.max_iterations),
        "tolerance": float(args.tolerance),
        "timeout_seconds": timeout_seconds,
        "output_files": ["summary.csv", "route_payloads.jsonl", "active_case.json"],
        "digitized_feed_sources": {
            str(figure): str(_source_feed_path(figure).relative_to(ANALYSIS_ROOT))
            for figure in requested_figures
        },
        "read_only_model_caches": {
            str(figure): str(_model_cache_path(figure).relative_to(ANALYSIS_ROOT))
            for figure in requested_figures
        },
        "status": "running",
        "started_at_utc": _utc_now(),
        "completed_rows": 0,
    }
    _write_json(metadata_path, metadata)
    for figure_number in requested_figures:
        for feed_row, raw_feed_row in _selected_feed_rows(
            figure_number,
            all_rows=bool(args.all_rows),
            max_rows_per_figure=int(args.max_rows_per_figure),
            tie_lines=tie_lines,
        ):
            active_case = {
                "status": "running",
                "started_at_utc": _utc_now(),
                "figure": figure_number,
                "tie_line": int(feed_row["tie_line"]),
                "timeout_seconds": timeout_seconds,
                "feed_source_path": str(_source_feed_path(figure_number).relative_to(ANALYSIS_ROOT)),
                "raw_digitized_feed_row": raw_feed_row,
                "feed_formula": feed_row["feed_formula"],
                "feed_ion_basis": formula_to_ion_basis(np.asarray(feed_row["feed_formula"], dtype=float)),
            }
            _write_json(active_case_path, active_case)
            print(f"[diagnose] figure {figure_number} tie-line {feed_row['tie_line']}", flush=True)
            summary, detail = _diagnose_feed(
                figure_number=figure_number,
                feed_row=feed_row,
                raw_feed_row=raw_feed_row,
                max_iterations=int(args.max_iterations),
                tolerance=float(args.tolerance),
                timeout_seconds=timeout_seconds,
            )
            summary_rows.append(summary)
            detail_rows.append(detail)
            active_case["status"] = "completed"
            active_case["completed_at_utc"] = _utc_now()
            active_case["route_accepted"] = summary["route_accepted"]
            active_case["route_status"] = summary["route_status"]
            active_case["solver_status"] = summary["solver_status"]
            active_case["application_status"] = summary["application_status"]
            _write_json(active_case_path, active_case)
            _write_csv(output_dir / "summary.csv", summary_rows)
            _write_jsonl(output_dir / "route_payloads.jsonl", detail_rows)
            metadata["completed_rows"] = len(summary_rows)
            metadata["last_completed_at_utc"] = _utc_now()
            _write_json(metadata_path, metadata)

    metadata["status"] = "completed"
    metadata["completed_at_utc"] = _utc_now()
    _write_json(metadata_path, metadata)
    print(f"[done] wrote {len(summary_rows)} diagnostic row(s) to {output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
