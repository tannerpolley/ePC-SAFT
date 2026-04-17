r"""Opt-in runtime profiling for MIAC fit generation.

This focuses on where time is spent inside the MIAC validation workflow:
payload generation, per-point state construction, activity-coefficient
evaluation, and plotting.

Run directly with:

    set ePCSAFT_RUN_PERF=1
    C:\ProgramData\Miniconda3\envs\ePC-SAFT\python.exe scripts\fits\profile_miac_runtime.py
"""

from __future__ import annotations

import csv
import os
import sys
import time
from pathlib import Path
from statistics import mean
from typing import Callable, Dict, List, Tuple

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts.fits import validate_miac_fits as vmf
from scripts._epcsaft_oop import as_mixture


REPORT_DIR = Path(__file__).resolve().parents[2] / "build" / "runtime_profile"
REPORT_CSV = REPORT_DIR / "miac_runtime_profile.csv"
REPORT_MD = REPORT_DIR / "miac_runtime_profile.md"


def _should_run_perf() -> bool:
    return os.environ.get("ePCSAFT_RUN_PERF", "").strip().lower() in {"1", "true", "yes", "on"}


def _requested_grid_points() -> int:
    token = os.environ.get("ePCSAFT_MIAC_PROFILE_GRID", "201").strip()
    try:
        value = int(token)
    except ValueError as exc:
        raise ValueError("ePCSAFT_MIAC_PROFILE_GRID must be an integer.") from exc
    return max(21, value)


def _requested_repeats() -> int:
    token = os.environ.get("ePCSAFT_MIAC_PROFILE_REPEATS", "5").strip()
    try:
        value = int(token)
    except ValueError as exc:
        raise ValueError("ePCSAFT_MIAC_PROFILE_REPEATS must be an integer.") from exc
    return max(1, value)


def _to_ms(seconds: float) -> float:
    return 1000.0 * float(seconds)


def _bench(name: str, fn: Callable[[], object], repeats: int, warmup: int = 1, meta: Dict[str, object] | None = None) -> dict:
    for _ in range(max(0, int(warmup))):
        fn()

    t0 = time.perf_counter()
    fn()
    first_ms = _to_ms(time.perf_counter() - t0)

    times_ms: List[float] = []
    for _ in range(max(1, int(repeats))):
        t1 = time.perf_counter()
        fn()
        times_ms.append(_to_ms(time.perf_counter() - t1))

    arr = np.asarray(times_ms, dtype=float)
    row = {
        "name": str(name),
        "repeats": int(arr.size),
        "first_ms": float(first_ms),
        "mean_ms": float(np.mean(arr)),
        "median_ms": float(np.median(arr)),
        "p95_ms": float(np.percentile(arr, 95)) if arr.size > 1 else float(arr[0]),
        "min_ms": float(np.min(arr)),
        "max_ms": float(np.max(arr)),
        "first_over_median": float(first_ms / max(float(np.median(arr)), 1.0e-12)),
    }
    if meta:
        row.update(meta)
    return row


def _find_combo(solvent_system: str, salt: str, comp_signature=None) -> Dict[str, object]:
    for combo in vmf.discover_combos(solvent_scope=solvent_system, salt_scope=salt):
        if comp_signature is None:
            if not combo.get("comp_signature"):
                return combo
        elif combo.get("comp_signature") == comp_signature:
            return combo
    raise LookupError(f"Could not locate combo for {solvent_system}/{salt} with comp_signature={comp_signature!r}")


def _profile_combo(
    combo: Dict[str, object],
    dataset_name: str,
    grid_points: int,
    repeats: int,
) -> Tuple[List[dict], List[str]]:
    salt = str(combo["salt"])
    solvent_system = str(combo["solvent_system"])
    label = f"{solvent_system}/{salt}"
    comp = dict(combo.get("comp", {}))
    molal_grid = np.linspace(0.0, 3.0, grid_points)
    species = vmf._species_for_combo(salt, solvent_system)
    user_options = dict(vmf.DATASET_VARIANTS[dataset_name].get("user_options", {}))
    params = vmf.build_params_for_variant(dataset_name, combo, user_options=user_options)
    mixture = as_mixture(params, species=species)
    pair_key = vmf._pair_key(salt)
    x_probe = vmf._molality_to_molefraction_combo(0.5, salt, solvent_system, comp)
    state_probe = mixture.state(T=vmf.T_REF, x=x_probe, P=vmf.P_REF, phase="liq")

    rows = [
        _bench(
            f"{label}.prepare_combo_payload[{grid_points}]",
            lambda: vmf.prepare_combo_payload(combo, grid_points=grid_points),
            repeats=max(1, min(3, repeats)),
            warmup=0,
            meta={"case": label, "dataset": dataset_name, "category": "payload"},
        )
    ]
    payload = vmf.prepare_combo_payload(combo, grid_points=grid_points)
    rows.extend(
        [
            _bench(
                f"{label}.plot_combo.miac_m",
                lambda: vmf.plot_combo(combo, save=False, close=True, quantity="miac_m", payload=payload),
                repeats=max(1, repeats),
                warmup=1,
                meta={"case": label, "dataset": dataset_name, "category": "plot"},
            ),
            _bench(
                f"{label}.plot_combo.miac",
                lambda: vmf.plot_combo(combo, save=False, close=True, quantity="miac", payload=payload),
                repeats=max(1, repeats),
                warmup=1,
                meta={"case": label, "dataset": dataset_name, "category": "plot"},
            ),
            _bench(
                f"{label}.point.mole_fraction_conversion",
                lambda: vmf._molality_to_molefraction_combo(0.5, salt, solvent_system, comp),
                repeats=max(5, repeats * 4),
                warmup=2,
                meta={"case": label, "dataset": dataset_name, "category": "point"},
            ),
            _bench(
                f"{label}.point.state_from_pressure",
                lambda: mixture.state(T=vmf.T_REF, x=x_probe, P=vmf.P_REF, phase="liq"),
                repeats=max(3, repeats),
                warmup=1,
                meta={"case": label, "dataset": dataset_name, "category": "point"},
            ),
            _bench(
                f"{label}.point.activity_component",
                lambda: state_probe.activity_coefficient(species=species),
                repeats=max(3, repeats),
                warmup=1,
                meta={"case": label, "dataset": dataset_name, "category": "point"},
            ),
            _bench(
                f"{label}.point.activity_mean_ionic_molality",
                lambda: state_probe.activity_coefficient(species=species, mean_ionic_form=True, basis="molality")[pair_key],
                repeats=max(3, repeats),
                warmup=1,
                meta={"case": label, "dataset": dataset_name, "category": "point"},
            ),
            _bench(
                f"{label}.point.fugacity_coefficient",
                lambda: state_probe.fugacity_coefficient(),
                repeats=max(5, repeats * 3),
                warmup=2,
                meta={"case": label, "dataset": dataset_name, "category": "point"},
            ),
            _bench(
                f"{label}.point.osmotic_coefficient",
                lambda: state_probe.osmotic_coefficient(),
                repeats=max(3, repeats),
                warmup=1,
                meta={"case": label, "dataset": dataset_name, "category": "point"},
            ),
        ]
    )

    def _row_named(suffix: str) -> dict:
        for row in rows:
            if row["name"].endswith(suffix):
                return row
        raise KeyError(suffix)

    payload_ms = _row_named(f"prepare_combo_payload[{grid_points}]")["mean_ms"]
    plot_ms = _row_named("plot_combo.miac_m")["mean_ms"]
    state_ms = _row_named("point.state_from_pressure")["mean_ms"]
    gamma_ms = _row_named("point.activity_mean_ionic_molality")["mean_ms"]
    x_ms = _row_named("point.mole_fraction_conversion")["mean_ms"]
    fug_ms = _row_named("point.fugacity_coefficient")["mean_ms"]
    osm_ms = _row_named("point.osmotic_coefficient")["mean_ms"]

    approx_curve_ms = float(grid_points) * (state_ms + gamma_ms)
    notes = [
        f"{label}: prepare_combo_payload[{grid_points}] averages {payload_ms:.1f} ms; the rendered plot itself is only about {plot_ms:.1f} ms, so plotting is not the bottleneck.",
        f"{label}: per-point work is dominated by state_from_pressure ({state_ms:.2f} ms) and activity_mean_ionic_molality ({gamma_ms:.2f} ms); mole-fraction conversion is only {x_ms:.3f} ms.",
        f"{label}: activity_mean_ionic_molality ({gamma_ms:.2f} ms) is roughly {gamma_ms / max(fug_ms, 1.0e-12):.0f}x the plain fugacity_coefficient call ({fug_ms:.3f} ms), which indicates the internal infinite-dilution reference-state work still dominates the activity path.",
        f"{label}: one model curve of {grid_points} points costs about {approx_curve_ms/1000.0:.2f} s from just the repeated pressure-state plus activity calls before Python bookkeeping or plotting.",
    ]
    if osm_ms > gamma_ms * 0.5:
        notes.append(f"{label}: osmotic_coefficient remains expensive ({osm_ms:.2f} ms) relative to activity, so any MIAC path that re-enters auxiliary electrolyte properties will be slow.")
    return rows, notes


def _write_reports(rows: List[dict], notes: List[str], grid_points: int, repeats: int) -> None:
    REPORT_DIR.mkdir(parents=True, exist_ok=True)

    fields = [
        "name",
        "case",
        "dataset",
        "category",
        "repeats",
        "first_ms",
        "mean_ms",
        "median_ms",
        "p95_ms",
        "min_ms",
        "max_ms",
        "first_over_median",
    ]
    with REPORT_CSV.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})

    lines = [
        "# MIAC runtime profile",
        "",
        f"- Grid points: {grid_points}",
        f"- Repeats: {repeats}",
        f"- Rows profiled: {len(rows)}",
        f"- CSV: `{REPORT_CSV}`",
        "",
        "## Slowest profiled calls",
        "",
        "| Call | Mean ms | Median ms | P95 ms | Category |",
        "|---|---:|---:|---:|---|",
    ]
    for row in sorted(rows, key=lambda item: item["mean_ms"], reverse=True)[:20]:
        lines.append(
            "| {name} | {mean_ms:.3f} | {median_ms:.3f} | {p95_ms:.3f} | {category} |".format(**row)
        )

    lines.extend(["", "## Notes", ""])
    for note in notes:
        lines.append(f"- {note}")

    REPORT_MD.write_text("\n".join(lines), encoding="utf-8")


def run_miac_runtime_profile() -> List[dict]:
    grid_points = _requested_grid_points()
    repeats = _requested_repeats()

    cases = [
        (_find_combo("water", "NaCl"), "2025_Figiel"),
        (_find_combo("methanol", "KCl"), "2025_Figiel"),
        (_find_combo("water-ethanol", "NaBr", comp_signature=(("water", 0.0), ("ethanol", 1.0))), "2025_Figiel"),
    ]

    rows: List[dict] = []
    notes: List[str] = []
    for combo, dataset_name in cases:
        case_rows, case_notes = _profile_combo(combo, dataset_name, grid_points, repeats)
        rows.extend(case_rows)
        notes.extend(case_notes)

    _write_reports(rows, notes, grid_points, repeats)
    print(f"MIAC runtime profile written:\n- {REPORT_CSV}\n- {REPORT_MD}")
    return rows


if __name__ == "__main__":
    if not _should_run_perf():
        raise SystemExit("Set ePCSAFT_RUN_PERF=1 to run MIAC runtime profiling.")
    run_miac_runtime_profile()
