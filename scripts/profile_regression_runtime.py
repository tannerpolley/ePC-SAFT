r"""Opt-in runtime profiling for native pure-neutral regression.

This compares:

- the current public native IPOPT-owned path
- the current internal native least-squares path
- the old Python/SciPy path from ``build/old-regression-bench`` when available

Run directly with:

    set ePCSAFT_RUN_PERF=1
    C:\ProgramData\Miniconda3\envs\ePC-SAFT\python.exe scripts\profile_regression_runtime.py
"""

from __future__ import annotations

import csv
import json
import os
import subprocess
import sys
import time
from pathlib import Path
from typing import Any

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from epcsaft import fit_pure_neutral
from epcsaft.regression import _fit_pure_neutral_least_squares_internal
from tests.test_regression import _load_workbook_reference_rows
from tests.test_regression import _neutral_fixed_parameters
from tests.test_regression import _real_saturation_records


REPORT_DIR = REPO_ROOT / "build" / "runtime_profile"
REPORT_CSV = REPORT_DIR / "regression_runtime_profile.csv"
REPORT_MD = REPORT_DIR / "regression_runtime_profile.md"
OLD_WORKTREE = REPO_ROOT / "build" / "old-regression-bench"


def _should_run_perf() -> bool:
    return os.environ.get("ePCSAFT_RUN_PERF", "").strip().lower() in {"1", "true", "yes", "on"}


def _benchmark_current_case(component: str, backend: str) -> dict[str, Any]:
    refs = _load_workbook_reference_rows()
    ref = refs[component]
    records = _real_saturation_records(component)
    kwargs = dict(
        records=records,
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
    solve = fit_pure_neutral if backend == "ipopt_native" else _fit_pure_neutral_least_squares_internal
    t0 = time.perf_counter()
    result = solve(**kwargs)
    elapsed = time.perf_counter() - t0
    return {
        "case": component,
        "backend": backend,
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
    }


def _benchmark_current_suite(backend: str) -> dict[str, Any]:
    t0 = time.perf_counter()
    rows = [_benchmark_current_case(component, backend) for component in ("Methane", "Ethane", "Propane")]
    elapsed = time.perf_counter() - t0
    return {
        "case": "hydrocarbon_suite",
        "backend": backend,
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
    }


def _run_old_subprocess(code: str) -> dict[str, Any]:
    completed = subprocess.run(
        [sys.executable, "-c", code],
        cwd=str(OLD_WORKTREE),
        capture_output=True,
        text=True,
        check=True,
    )
    stdout = completed.stdout.strip().splitlines()
    if not stdout:
        raise RuntimeError("Old regression subprocess produced no JSON output.")
    return json.loads(stdout[-1])


def _benchmark_old_case(component: str) -> dict[str, Any]:
    old_root = str(OLD_WORKTREE)
    old_src = str(OLD_WORKTREE / "src")
    code = f"""
import json, sys, time
sys.path.insert(0, {old_root!r})
sys.path.insert(0, {old_src!r})
from tests.test_regression import _load_workbook_reference_rows, _neutral_fixed_parameters, _real_saturation_records
from epcsaft import fit_pure_neutral
ref = _load_workbook_reference_rows()[{component!r}]
records = _real_saturation_records({component!r})
t0 = time.perf_counter()
result = fit_pure_neutral(
    records,
    {component!r},
    assoc_scheme='',
    fixed_parameters=_neutral_fixed_parameters({component!r}),
    initial_guess={{'m': ref['m'] * 1.08, 's': ref['s'] * 0.96, 'e': ref['e'] * 1.05}},
    bounds={{'m': (0.5, 3.5), 's': (2.0, 5.0), 'e': (50.0, 400.0)}},
)
print(json.dumps({{
    'case': {component!r},
    'backend': 'python_scipy_legacy',
    'wall_s': time.perf_counter() - t0,
    'nfev': int(result.nfev),
    'success': bool(result.success),
    'status': int(result.status),
    'message': str(result.message),
    'm': float(result.fitted_values['m']),
    's': float(result.fitted_values['s']),
    'e': float(result.fitted_values['e']),
    'density_rms': float(result.metrics_by_term['density']),
    'pure_vle_rms': float(result.metrics_by_term['pure_vle_fugacity_balance']),
}}))
"""
    return _run_old_subprocess(code)


def _benchmark_old_suite() -> dict[str, Any]:
    old_root = str(OLD_WORKTREE)
    old_src = str(OLD_WORKTREE / "src")
    code = f"""
import json, sys, time
sys.path.insert(0, {old_root!r})
sys.path.insert(0, {old_src!r})
from tests.test_regression import _load_workbook_reference_rows, _neutral_fixed_parameters, _real_saturation_records
from epcsaft import fit_pure_neutral
refs = _load_workbook_reference_rows()
rows = []
t0 = time.perf_counter()
for component in ('Methane', 'Ethane', 'Propane'):
    ref = refs[component]
    records = _real_saturation_records(component)
    result = fit_pure_neutral(
        records,
        component,
        assoc_scheme='',
        fixed_parameters=_neutral_fixed_parameters(component),
        initial_guess={{'m': ref['m'] * 1.08, 's': ref['s'] * 0.96, 'e': ref['e'] * 1.05}},
        bounds={{'m': (0.5, 3.5), 's': (2.0, 5.0), 'e': (50.0, 400.0)}},
    )
    rows.append({{
        'density_rms': float(result.metrics_by_term['density']),
        'pure_vle_rms': float(result.metrics_by_term['pure_vle_fugacity_balance']),
        'nfev': int(result.nfev),
        'success': bool(result.success),
    }})
print(json.dumps({{
    'case': 'hydrocarbon_suite',
    'backend': 'python_scipy_legacy',
    'wall_s': time.perf_counter() - t0,
    'nfev': int(sum(row['nfev'] for row in rows)),
    'success': all(row['success'] for row in rows),
    'status': 0,
    'message': 'suite',
    'm': float('nan'),
    's': float('nan'),
    'e': float('nan'),
    'density_rms': float(max(row['density_rms'] for row in rows)),
    'pure_vle_rms': float(max(row['pure_vle_rms'] for row in rows)),
}}))
"""
    return _run_old_subprocess(code)


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
    ]
    with REPORT_CSV.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in fieldnames})

    lines = [
        "# Regression Runtime Profile",
        "",
        "| Case | Backend | Wall (s) | NFEV | Success | Density RMS | Pure VLE RMS | m | s | e |",
        "| --- | --- | ---: | ---: | --- | ---: | ---: | ---: | ---: | ---: |",
    ]
    for row in rows:
        lines.append(
            "| "
            + " | ".join(
                [
                    str(row["case"]),
                    str(row["backend"]),
                    _format_float(row["wall_s"]),
                    str(row["nfev"]),
                    "yes" if row["success"] else "no",
                    _format_float(row["density_rms"]),
                    _format_float(row["pure_vle_rms"]),
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
    rows.append(_benchmark_current_case("Methane", "ipopt_native"))
    rows.append(_benchmark_current_case("Methane", "least_squares_native"))
    rows.append(_benchmark_current_suite("ipopt_native"))
    rows.append(_benchmark_current_suite("least_squares_native"))

    if OLD_WORKTREE.exists():
        rows.append(_benchmark_old_case("Methane"))
        rows.append(_benchmark_old_suite())

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
