"""Performance profiling test for the object-oriented PC-SAFT API.

This module is intentionally opt-in. Run with:

    set PCSAFT_RUN_PERF=1
    python -m pytest tests/test_runtime_profile.py -s
"""

from __future__ import annotations

import csv
import inspect
import os
import time
from copy import deepcopy
from pathlib import Path

import numpy as np
import pytest

from pcsaft import PCSAFTMixture


REPORT_DIR = Path(__file__).resolve().parents[1] / "build" / "runtime_profile"
REPORT_CSV = REPORT_DIR / "runtime_profile.csv"
REPORT_MD = REPORT_DIR / "runtime_profile.md"


def _should_run_perf() -> bool:
    return os.environ.get("PCSAFT_RUN_PERF", "").strip().lower() in {"1", "true", "yes", "on"}


def _to_ms(seconds: float) -> float:
    return 1000.0 * float(seconds)


def _bench(name: str, fn, repeats: int, warmup: int = 1) -> dict:
    for _ in range(max(0, int(warmup))):
        fn()

    t0 = time.perf_counter()
    fn()
    first_ms = _to_ms(time.perf_counter() - t0)

    times_ms = []
    for _ in range(max(1, int(repeats))):
        t1 = time.perf_counter()
        fn()
        times_ms.append(_to_ms(time.perf_counter() - t1))

    arr = np.asarray(times_ms, dtype=float)
    median_ms = float(np.median(arr))
    p95_ms = float(np.percentile(arr, 95)) if arr.size > 1 else float(arr[0])
    return {
        "name": str(name),
        "repeats": int(arr.size),
        "first_ms": float(first_ms),
        "mean_ms": float(np.mean(arr)),
        "median_ms": median_ms,
        "p95_ms": p95_ms,
        "min_ms": float(np.min(arr)),
        "max_ms": float(np.max(arr)),
        "first_over_median": float(first_ms / max(median_ms, 1.0e-12)),
    }


def _write_reports(rows: list[dict], bottleneck_notes: list[str]) -> None:
    REPORT_DIR.mkdir(parents=True, exist_ok=True)

    fields = [
        "name",
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
            writer.writerow(row)

    lines = [
        "# PC-SAFT runtime profile",
        "",
        f"- Rows profiled: {len(rows)}",
        f"- CSV: `{REPORT_CSV}`",
        "",
        "## Slowest by mean runtime (ms)",
        "",
        "| Method | First | Mean | Median | P95 | first/median |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for row in sorted(rows, key=lambda item: item["mean_ms"], reverse=True)[:15]:
        lines.append(
            "| {name} | {first_ms:.3f} | {mean_ms:.3f} | {median_ms:.3f} | {p95_ms:.3f} | {first_over_median:.2f} |".format(
                **row
            )
        )

    lines.extend(["", "## Bottleneck notes", ""])
    if bottleneck_notes:
        for note in bottleneck_notes:
            lines.append(f"- {note}")
    else:
        lines.append("- No obvious bottleneck heuristic triggered.")

    REPORT_MD.write_text("\n".join(lines), encoding="utf-8")


def _public_callables(cls) -> set[str]:
    names = set()
    for name, member in inspect.getmembers(cls, predicate=callable):
        if name.startswith("_"):
            continue
        names.add(name)
    return names


def _build_states() -> dict:
    neutral_species = ["Toluene"]
    neutral_params = {
        "m": np.asarray([2.8149]),
        "s": np.asarray([3.7169]),
        "e": np.asarray([285.69]),
    }
    neutral_mix = PCSAFTMixture.from_params(deepcopy(neutral_params), species=neutral_species)
    neutral_state_tp = neutral_mix.state(T=325.0, x=np.asarray([1.0]), P=101325.0, phase="liq")
    neutral_state_trho = neutral_mix.state(T=325.0, x=np.asarray([1.0]), rho=neutral_state_tp.density(), phase="liq")

    ionic_species = ["H2O-2B-Li", "Na+", "Cl-"]
    ionic_params = {
        "MW": np.asarray([18.01528e-3, 22.98e-3, 35.45e-3]),
        "m": np.asarray([1.2047, 1.0, 1.0]),
        "s": np.asarray([3.0, 2.8232, 2.7560]),
        "e": np.asarray([353.95, 230.0, 170.0]),
        "e_assoc": np.asarray([2425.7, 0.0, 0.0]),
        "vol_a": np.asarray([0.04509, 0.0, 0.0]),
        "assoc_scheme": ["2B", None, None],
        "dipm": np.asarray([0.0, 0.0, 0.0]),
        "dip_num": np.asarray([1.0, 1.0, 1.0]),
        "z": np.asarray([0.0, 1.0, -1.0]),
        "dielc": np.asarray([78.09, 8.0, 8.0]),
        "d_born": np.asarray([0.0, 3.445, 4.1]),
        "f_solv": np.asarray([1.5, 1.0, 1.0]),
        "k_ij": np.asarray(
            [
                [0.0, 0.0045, -0.25],
                [0.0045, 0.0, 0.317],
                [-0.25, 0.317, 0.0],
            ]
        ),
        "l_ij": np.zeros((3, 3)),
        "k_hb": np.zeros((3, 3)),
        "elec_model": {
            "rel_perm": {"rule": 1, "differential_mode": 0},
            "DH_model": {
                "d_ion_mode": 1,
                "bjeruum_treatment": False,
                "mu_DH_model": {
                    "differential_mode": 0,
                    "comp_dep_rel_perm": True,
                    "include_sum_term": True,
                },
            },
            "include_born_model": True,
            "born_model": {
                "d_Born_mode": 0,
                "solvation_shell_model": False,
                "dielectric_saturation": False,
                "bulk_mode": 0,
                "mu_born_model": {
                    "differential_mode": 0,
                    "comp_dep_rel_perm": True,
                    "include_sum_term": True,
                    "comp_dep_delta_d": False,
                },
            },
        },
        "debug": False,
    }
    z_feed = np.asarray([1.0 / 0.01801528, 1.0e-4, 1.0e-4], dtype=float)
    z_feed = z_feed / np.sum(z_feed)
    ionic_mix = PCSAFTMixture.from_params(deepcopy(ionic_params), species=ionic_species)
    ionic_state_tp = ionic_mix.state(T=298.15, x=z_feed, P=1.0e5, phase="liq")
    ionic_state_trho = ionic_mix.state(T=298.15, x=z_feed, rho=ionic_state_tp.density(), phase="liq")

    return {
        "neutral_species": neutral_species,
        "neutral_params": neutral_params,
        "neutral_mix": neutral_mix,
        "neutral_state_tp": neutral_state_tp,
        "neutral_state_trho": neutral_state_trho,
        "ionic_species": ionic_species,
        "ionic_params": ionic_params,
        "ionic_mix": ionic_mix,
        "ionic_state_tp": ionic_state_tp,
        "ionic_state_trho": ionic_state_trho,
        "z_feed": z_feed,
    }


def test_runtime_profile_oop_methods():
    if not _should_run_perf():
        pytest.skip("Set PCSAFT_RUN_PERF=1 to run runtime profiling.")

    ctx = _build_states()
    neutral_state_tp = ctx["neutral_state_tp"]
    neutral_state_trho = ctx["neutral_state_trho"]
    ionic_state_tp = ctx["ionic_state_tp"]
    ionic_state_trho = ctx["ionic_state_trho"]
    ionic_species = ctx["ionic_species"]
    z_feed = ctx["z_feed"]

    aly_lee_params = np.asarray([1.2, 0.8, -0.01, 2.0e-5, -3.0e-8], dtype=float)
    benches = [
        ("ctor.from_params", lambda: PCSAFTMixture.from_params(deepcopy(ctx["neutral_params"]), species=ctx["neutral_species"]), 3, 0),
        ("ctor.from_dataset.water", lambda: PCSAFTMixture.from_dataset("2012_Held", ["Water"], np.asarray([1.0]), 298.15), 3, 0),
        ("state.from_P", lambda: ctx["neutral_mix"].state(T=325.0, x=np.asarray([1.0]), P=101325.0, phase="liq"), 8, 1),
        ("state.from_rho", lambda: ctx["neutral_mix"].state(T=325.0, x=np.asarray([1.0]), rho=neutral_state_tp.density(), phase="liq"), 8, 1),
        ("state.pressure", lambda: neutral_state_trho.pressure(), 25, 2),
        ("state.density", lambda: neutral_state_tp.density(), 25, 2),
        ("state.Z", lambda: neutral_state_trho.Z(), 20, 2),
        ("state.a_res", lambda: neutral_state_trho.a_res(), 20, 2),
        ("state.dadt", lambda: neutral_state_trho.dadt(), 12, 2),
        ("state.h_res", lambda: neutral_state_trho.h_res(), 20, 2),
        ("state.s_res", lambda: neutral_state_trho.s_res(), 20, 2),
        ("state.g_res", lambda: neutral_state_trho.g_res(), 20, 2),
        ("state.mu_res", lambda: neutral_state_trho.mu_res(), 20, 2),
        ("state.gamma", lambda: neutral_state_trho.gamma(), 20, 2),
        ("state.lnfugcoef", lambda: neutral_state_trho.lnfugcoef(), 20, 2),
        ("state.lnfugcoef_terms", lambda: neutral_state_trho.lnfugcoef_terms(), 10, 1),
        ("state.breakdown.neutral", lambda: neutral_state_trho.breakdown(species=ctx["neutral_species"]), 6, 1),
        ("state.cp", lambda: neutral_state_tp.cp(aly_lee_params), 6, 1),
        ("state.flashTQ", lambda: neutral_state_tp.flashTQ(q=0.0), 3, 0),
        ("state.flashPQ", lambda: neutral_state_tp.flashPQ(p=101325.0, q=0.0), 3, 0),
        ("state.Hvap", lambda: neutral_state_tp.Hvap(), 3, 0),
        ("ionic.dielectric_eval", lambda: ionic_state_tp.dielectric_eval(), 10, 1),
        ("ionic.osmoticC", lambda: ionic_state_tp.osmoticC(), 10, 1),
        ("ionic.actcoeff", lambda: ionic_state_tp.actcoeff(species=ionic_species), 10, 1),
        ("ionic.miac_m", lambda: ionic_state_trho.miac_m(species=ionic_species), 8, 1),
        ("ionic.miac", lambda: ionic_state_trho.miac(species=ionic_species), 8, 1),
        ("ionic.gsolv", lambda: ionic_state_trho.gsolv(species=ionic_species), 8, 1),
        ("ionic.breakdown", lambda: ionic_state_tp.breakdown(species=ionic_species), 4, 0),
        (
            "ionic.miac_m.full_rebuild",
            lambda: PCSAFTMixture.from_params(deepcopy(ctx["ionic_params"]), species=ionic_species)
            .state(T=298.15, x=z_feed, P=1.0e5, phase="liq")
            .miac_m(species=ionic_species),
            3,
            0,
        ),
    ]

    profiled_aliases = {
        "PCSAFTMixture.from_params": "ctor.from_params",
        "PCSAFTMixture.from_dataset": "ctor.from_dataset.water",
        "PCSAFTMixture.state": "state.from_P",
        "PCSAFTState.pressure": "state.pressure",
        "PCSAFTState.density": "state.density",
        "PCSAFTState.Z": "state.Z",
        "PCSAFTState.a_res": "state.a_res",
        "PCSAFTState.dadt": "state.dadt",
        "PCSAFTState.h_res": "state.h_res",
        "PCSAFTState.s_res": "state.s_res",
        "PCSAFTState.g_res": "state.g_res",
        "PCSAFTState.mu_res": "state.mu_res",
        "PCSAFTState.gamma": "state.gamma",
        "PCSAFTState.lnfugcoef": "state.lnfugcoef",
        "PCSAFTState.lnfugcoef_terms": "state.lnfugcoef_terms",
        "PCSAFTState.breakdown": "state.breakdown.neutral",
        "PCSAFTState.cp": "state.cp",
        "PCSAFTState.flashTQ": "state.flashTQ",
        "PCSAFTState.flashPQ": "state.flashPQ",
        "PCSAFTState.Hvap": "state.Hvap",
        "PCSAFTState.dielectric_eval": "ionic.dielectric_eval",
        "PCSAFTState.osmoticC": "ionic.osmoticC",
        "PCSAFTState.actcoeff": "ionic.actcoeff",
        "PCSAFTState.miac_m": "ionic.miac_m",
        "PCSAFTState.miac": "ionic.miac",
        "PCSAFTState.gsolv": "ionic.gsolv",
    }
    bench_names = {name for name, _, _, _ in benches}
    missing_coverage = [api for api, alias in profiled_aliases.items() if alias not in bench_names]
    assert not missing_coverage, f"Missing benchmark coverage for API methods: {missing_coverage}"

    expected_state_methods = _public_callables(type(neutral_state_tp))
    expected_mixture_methods = _public_callables(type(ctx["neutral_mix"]))
    missing_state_methods = sorted(
        method for method in expected_state_methods if f"PCSAFTState.{method}" not in profiled_aliases
    )
    missing_mixture_methods = sorted(
        method for method in expected_mixture_methods if f"PCSAFTMixture.{method}" not in profiled_aliases
    )
    assert not missing_state_methods, (
        "Unprofiled public PCSAFTState methods found; extend tests/test_runtime_profile.py: "
        + ", ".join(missing_state_methods)
    )
    assert not missing_mixture_methods, (
        "Unprofiled public PCSAFTMixture methods found; extend tests/test_runtime_profile.py: "
        + ", ".join(missing_mixture_methods)
    )

    rows = [_bench(name, fn, repeats=repeats, warmup=warmup) for name, fn, repeats, warmup in benches]

    slow_rows = [row for row in rows if row["mean_ms"] >= 50.0]
    init_heavy = [row for row in rows if row["first_over_median"] >= 3.0 and row["median_ms"] >= 1.0]
    bottleneck_notes = []
    if slow_rows:
        bottleneck_notes.append(
            "High-mean methods (>=50 ms): "
            + ", ".join(row["name"] for row in sorted(slow_rows, key=lambda item: item["mean_ms"], reverse=True)[:8])
        )
    if init_heavy:
        bottleneck_notes.append(
            "High first-call overhead (first/median >= 3): "
            + ", ".join(row["name"] for row in sorted(init_heavy, key=lambda item: item["first_over_median"], reverse=True)[:8])
        )

    rebuild_row = next((row for row in rows if row["name"] == "ionic.miac_m.full_rebuild"), None)
    reuse_row = next((row for row in rows if row["name"] == "ionic.miac_m"), None)
    if rebuild_row is not None and reuse_row is not None and reuse_row["mean_ms"] > 0.0:
        ratio = rebuild_row["mean_ms"] / reuse_row["mean_ms"]
        bottleneck_notes.append(
            "Full rebuild vs reused-state MIAC ratio: {:.2f}x ({} / {}).".format(
                ratio, "ionic.miac_m.full_rebuild", "ionic.miac_m"
            )
        )

    _write_reports(rows, bottleneck_notes)

    print("\nRuntime profile report written:")
    print(f"- {REPORT_CSV}")
    print(f"- {REPORT_MD}")
    print("\nTop 10 by mean runtime (ms):")
    for row in sorted(rows, key=lambda item: item["mean_ms"], reverse=True)[:10]:
        print(
            "{:<30} mean={:>9.3f}  median={:>9.3f}  first={:>9.3f}  p95={:>9.3f}".format(
                row["name"], row["mean_ms"], row["median_ms"], row["first_ms"], row["p95_ms"]
            )
        )

    assert len(rows) >= 20
    assert REPORT_CSV.exists()
    assert REPORT_MD.exists()
