r"""Opt-in runtime profiling for native pure-neutral regression.

This compares:

- the current public native least-squares workflow
- the internal native least-squares path

Run directly with:

    set ePCSAFT_RUN_PERF=1
    uv run python scripts\benchmarks\profile_regression_runtime.py
"""

from __future__ import annotations

import csv
import json
import os
import sys
import time
from pathlib import Path
from typing import Any

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from epcsaft import fit_pure_neutral
from epcsaft.regression import _fit_pure_neutral_least_squares_internal
from tests.helpers.regression_cases import _load_workbook_reference_rows
from tests.helpers.regression_cases import _neutral_fixed_parameters
from tests.helpers.regression_cases import _real_saturation_records

REPORT_DIR = REPO_ROOT / "build" / "runtime_profile"
REPORT_CSV = REPORT_DIR / "regression_runtime_profile.csv"
REPORT_MD = REPORT_DIR / "regression_runtime_profile.md"


def _should_run_perf() -> bool:
    return os.environ.get("ePCSAFT_RUN_PERF", "").strip().lower() in {"1", "true", "yes", "on"}


def _benchmark_kwargs(component: str) -> dict[str, Any]:
    refs = _load_workbook_reference_rows()
    ref = refs[component]
    return dict(
        records=_real_saturation_records(component),
        component=component,
        assoc_scheme="",
        fixed_parameters=_neutral_fixed_parameters(component),
        initial_guess={
            "m": ref["m"] * 1.08,
            "s": ref["s"] * 0.96,
            "e": ref["e"] * 1.05,
        },
        bounds={
            "m": (0.5, 3.5),
            "s": (2.0, 5.0),
            "e": (50.0, 400.0),
        },
    )


def _benchmark_current_case(component: str, backend: str) -> dict[str, Any]:
    kwargs = _benchmark_kwargs(component)
    if backend == "public_default":
        solve = fit_pure_neutral
    elif backend == "least_squares_native":
        solve = _fit_pure_neutral_least_squares_internal
    else:
        raise ValueError(f"Unsupported benchmark backend {backend!r}")

    t0 = time.perf_counter()
    result = solve(**kwargs)
    elapsed = time.perf_counter() - t0

    workflow = {"selected": backend}
    chosen_diag = {
        "starts_tried": 0,
        "initial_cost": float("nan"),
    }
    return {
        "case": component,
        "backend": backend,
        "returned_backend": str(result.backend),
        "workflow_selected": str(workflow.get("selected", "")),
        "wall_s": float(elapsed),
        "nfev": int(result.nfev),
        "success": bool(result.success),
        "status": int(result.status),
        "message": str(result.message),
        "m": float(result.fitted_values["m"]),
        "s": float(result.fitted_values["s"]),
        "e": float(result.fitted_values["e"]),
        "density_rms": float(result.metrics_by_term["density"]),
        "pure_vle_rms": float(result.metrics_by_term["pure_vle_fugacity_balance"]),
        "starts_tried": int(chosen_diag.get("starts_tried", 0)),
        "initial_cost": float(chosen_diag.get("initial_cost", float("nan"))),
    }


def _benchmark_current_suite(backend: str) -> dict[str, Any]:
    t0 = time.perf_counter()
    rows = [_benchmark_current_case(component, backend) for component in ("Methane", "Ethane", "Propane")]
    elapsed = time.perf_counter() - t0
    return {
        "case": "hydrocarbon_suite",
        "backend": backend,
        "returned_backend": backend,
        "workflow_selected": backend,
        "wall_s": float(elapsed),
        "nfev": int(sum(int(row["nfev"]) for row in rows)),
        "success": all(bool(row["success"]) for row in rows),
        "status": 0,
        "message": "suite",
        "m": float("nan"),
        "s": float("nan"),
        "e": float("nan"),
        "density_rms": float(max(float(row["density_rms"]) for row in rows)),
        "pure_vle_rms": float(max(float(row["pure_vle_rms"]) for row in rows)),
        "starts_tried": int(sum(int(row.get("starts_tried", 0)) for row in rows)),
        "initial_cost": float("nan"),
    }


def _format_float(value: Any) -> str:
    try:
        numeric = float(value)
    except (TypeError, ValueError):
        return str(value)
    if numeric != numeric:
        return "nan"
    return f"{numeric:.6g}"


def _write_reports(rows: list[dict[str, Any]]) -> None:
    REPORT_DIR.mkdir(parents=True, exist_ok=True)

    fieldnames = [
        "case",
        "backend",
        "returned_backend",
        "workflow_selected",
        "wall_s",
        "nfev",
        "success",
        "status",
        "message",
        "m",
        "s",
        "e",
        "density_rms",
        "pure_vle_rms",
        "starts_tried",
        "initial_cost",
    ]
    with REPORT_CSV.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in fieldnames})

    lines = [
        "# Regression Runtime Profile",
        "",
        "| Case | Backend | Returned | Workflow | Wall (s) | NFEV | Starts | Success | Density RMS | Pure VLE RMS | Initial Cost | m | s | e |",
        "| --- | --- | --- | --- | ---: | ---: | ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for row in rows:
        lines.append(
            "| "
            + " | ".join(
                [
                    str(row["case"]),
                    str(row["backend"]),
                    str(row["returned_backend"]),
                    str(row["workflow_selected"]),
                    _format_float(row["wall_s"]),
                    str(row["nfev"]),
                    str(row.get("starts_tried", 0)),
                    "yes" if row["success"] else "no",
                    _format_float(row["density_rms"]),
                    _format_float(row["pure_vle_rms"]),
                    _format_float(row.get("initial_cost", float("nan"))),
                    _format_float(row["m"]),
                    _format_float(row["s"]),
                    _format_float(row["e"]),
                ]
            )
            + " |"
        )
    REPORT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run_regression_runtime_profile() -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    rows.append(_benchmark_current_case("Methane", "public_default"))
    rows.append(_benchmark_current_case("Methane", "least_squares_native"))
    rows.append(_benchmark_current_suite("public_default"))
    rows.append(_benchmark_current_suite("least_squares_native"))

    _write_reports(rows)
    return rows


def main() -> int:
    if not _should_run_perf():
        print("Set ePCSAFT_RUN_PERF=1 to run regression runtime profiling.")
        return 0
    rows = run_regression_runtime_profile()
    print(json.dumps(rows, indent=2))
    print(f"Wrote {REPORT_CSV}")
    print(f"Wrote {REPORT_MD}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
