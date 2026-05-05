"""Performance profiling test for the object-oriented ePC-SAFT API.

This module is intentionally opt-in. Run with:

    set ePCSAFT_RUN_PERF=1
    python -m pytest tests/profile/test_runtime_profile.py -s
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

from epcsaft import ePCSAFTMixture

REPORT_DIR = Path(__file__).resolve().parents[2] / "build" / "runtime_profile"
REPORT_CSV = REPORT_DIR / "runtime_profile.csv"
REPORT_MD = REPORT_DIR / "runtime_profile.md"


def _should_run_perf() -> bool:
    return os.environ.get("ePCSAFT_RUN_PERF", "").strip().lower() in {"1", "true", "yes", "on"}


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
        "# ePC-SAFT runtime profile",
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


def _full_rebuild_reuse_note(rebuild_row: dict, reuse_row: dict, ratio: float) -> str:
    return (
        "Full rebuild vs reused-state MIAC ratio: {:.2f}x ({} / {}). "
        "Optimization hint: reuse ePCSAFTMixture and ePCSAFTState objects in hot loops instead of rebuilding them."
    ).format(ratio, rebuild_row["name"], reuse_row["name"])


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
    neutral_mix = ePCSAFTMixture.from_params(deepcopy(neutral_params), species=neutral_species)
    # Pressure-based states now pay the T, P -> rho closure cost during construction.
    neutral_state_tp = neutral_mix.state(T=325.0, x=np.asarray([1.0]), P=101325.0, phase="liq")
    neutral_state_trho = neutral_mix.state(T=325.0, x=np.asarray([1.0]), rho=neutral_state_tp.density(), phase="liq")

    equilibrium_species = ["Methane", "Ethane", "Propane"]
    equilibrium_params = {
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
    equilibrium_mix = ePCSAFTMixture.from_params(deepcopy(equilibrium_params), species=equilibrium_species)
    equilibrium_feed = np.asarray([0.1, 0.3, 0.6], dtype=float)

    ionic_species = ["H2O-2B-Li", "Na+", "Cl-"]
    ionic_params = {
        "MW": np.asarray([18.01528e-3, 22.98e-3, 35.45e-3]),
        "m": np.asarray([1.2047, 1.0, 1.0]),
        "s": np.asarray([3.0, 2.8232, 2.7560]),
        "e": np.asarray([353.95, 230.0, 170.0]),
        "e_assoc": np.asarray([2425.7, 0.0, 0.0]),
        "vol_a": np.asarray([0.04509, 0.0, 0.0]),
        "assoc_scheme": ["2B", None, None],
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
    }
    z_feed = np.asarray([1.0 / 0.01801528, 1.0e-4, 1.0e-4], dtype=float)
    z_feed = z_feed / np.sum(z_feed)
    ionic_mix = ePCSAFTMixture.from_params(deepcopy(ionic_params), species=ionic_species)
    ionic_state_tp = ionic_mix.state(T=298.15, x=z_feed, P=1.0e5, phase="liq")
    ionic_state_trho = ionic_mix.state(T=298.15, x=z_feed, rho=ionic_state_tp.density(), phase="liq")

    return {
        "neutral_species": neutral_species,
        "neutral_params": neutral_params,
        "neutral_mix": neutral_mix,
        "neutral_state_tp": neutral_state_tp,
        "neutral_state_trho": neutral_state_trho,
        "equilibrium_mix": equilibrium_mix,
        "equilibrium_feed": equilibrium_feed,
        "ionic_species": ionic_species,
        "ionic_params": ionic_params,
        "ionic_mix": ionic_mix,
        "ionic_state_tp": ionic_state_tp,
        "ionic_state_trho": ionic_state_trho,
        "z_feed": z_feed,
    }


def test_runtime_profile_oop_methods():
    if not _should_run_perf():
        pytest.skip("Set ePCSAFT_RUN_PERF=1 to run runtime profiling.")

    ctx = _build_states()
    neutral_state_tp = ctx["neutral_state_tp"]
    neutral_state_trho = ctx["neutral_state_trho"]
    ionic_state_tp = ctx["ionic_state_tp"]
    ionic_state_trho = ctx["ionic_state_trho"]
    ionic_species = ctx["ionic_species"]
    z_feed = ctx["z_feed"]

    benches = [
        (
            "ctor.from_params",
            lambda: ePCSAFTMixture.from_params(deepcopy(ctx["neutral_params"]), species=ctx["neutral_species"]),
            3,
            0,
        ),
        (
            "ctor.from_dataset.water",
            lambda: ePCSAFTMixture.from_dataset("2012_Held", ["Water"], np.asarray([1.0]), 298.15),
            3,
            0,
        ),
        ("mixture.clear_runtime_caches", lambda: ctx["neutral_mix"].clear_runtime_caches(), 10, 1),
        ("mixture.reset_runtime_cache_stats", lambda: ctx["neutral_mix"].reset_runtime_cache_stats(), 10, 1),
        ("mixture.runtime_cache_stats", lambda: ctx["neutral_mix"].runtime_cache_stats(), 10, 1),
        (
            "mixture.equilibrium.tp_flash",
            lambda: ctx["equilibrium_mix"].equilibrium(
                kind="tp_flash",
                T=220.0,
                P=1.0e5,
                z=ctx["equilibrium_feed"],
            ),
            3,
            0,
        ),
        ("state.from_P", lambda: ctx["neutral_mix"].state(T=325.0, x=np.asarray([1.0]), P=101325.0, phase="liq"), 8, 1),
        (
            "state.from_rho",
            lambda: ctx["neutral_mix"].state(T=325.0, x=np.asarray([1.0]), rho=neutral_state_tp.density(), phase="liq"),
            8,
            1,
        ),
        ("state.pressure", lambda: neutral_state_trho.pressure(), 25, 2),
        ("state.density", lambda: neutral_state_tp.density(), 25, 2),
        ("state.molar_density", lambda: neutral_state_trho.molar_density(), 25, 2),
        ("ionic.mass_density", lambda: ionic_state_trho.mass_density(), 25, 2),
        ("state.method_aliases", lambda: neutral_state_trho.method_aliases(), 25, 2),
        ("state.compressibility_factor", lambda: neutral_state_trho.compressibility_factor(), 20, 2),
        (
            "state.compressibility_factor.terms",
            lambda: neutral_state_trho.compressibility_factor(return_contribution_terms=True),
            10,
            1,
        ),
        ("state.residual_helmholtz", lambda: neutral_state_trho.residual_helmholtz(), 20, 2),
        (
            "state.residual_helmholtz.terms",
            lambda: neutral_state_trho.residual_helmholtz(return_contribution_terms=True),
            10,
            1,
        ),
        (
            "state.temperature_derivative_residual_helmholtz",
            lambda: neutral_state_trho.temperature_derivative_residual_helmholtz(),
            12,
            2,
        ),
        (
            "state.temperature_derivative_residual_helmholtz.terms",
            lambda: neutral_state_trho.temperature_derivative_residual_helmholtz(return_contribution_terms=True),
            8,
            1,
        ),
        ("state.residual_enthalpy", lambda: neutral_state_trho.residual_enthalpy(), 20, 2),
        ("state.residual_entropy", lambda: neutral_state_trho.residual_entropy(), 20, 2),
        ("state.residual_gibbs", lambda: neutral_state_trho.residual_gibbs(), 20, 2),
        (
            "state.composition_derivative_residual_helmholtz",
            lambda: neutral_state_trho.composition_derivative_residual_helmholtz(),
            10,
            1,
        ),
        ("state.residual_chemical_potential", lambda: neutral_state_trho.residual_chemical_potential(), 20, 2),
        (
            "state.residual_chemical_potential.terms",
            lambda: neutral_state_trho.residual_chemical_potential(return_contribution_terms=True),
            10,
            1,
        ),
        ("state.fugacity_coefficient", lambda: neutral_state_trho.fugacity_coefficient(), 20, 2),
        (
            "state.fugacity_coefficient.coefficient",
            lambda: neutral_state_trho.fugacity_coefficient(natural_log=False),
            20,
            2,
        ),
        (
            "state.fugacity_coefficient.terms",
            lambda: neutral_state_trho.fugacity_coefficient(return_contribution_terms=True),
            8,
            1,
        ),
        (
            "state.state_diagnostics.neutral",
            lambda: neutral_state_trho.state_diagnostics(species=ctx["neutral_species"]),
            6,
            1,
        ),
        ("ionic.relative_permittivity", lambda: ionic_state_tp.relative_permittivity(), 10, 1),
        ("ionic.osmotic_coefficient", lambda: ionic_state_tp.osmotic_coefficient(), 10, 1),
        (
            "ionic.activity_coefficient.component",
            lambda: ionic_state_tp.activity_coefficient(species=ionic_species),
            10,
            1,
        ),
        (
            "ionic.activity_coefficient.mean_ionic_molality",
            lambda: ionic_state_trho.activity_coefficient(
                species=ionic_species, mean_ionic_form=True, basis="molality"
            ),
            8,
            1,
        ),
        (
            "ionic.activity_coefficient.mean_ionic_mole",
            lambda: ionic_state_trho.activity_coefficient(species=ionic_species, mean_ionic_form=True, basis="mole"),
            8,
            1,
        ),
        ("ionic.solvation_free_energy", lambda: ionic_state_trho.solvation_free_energy(species=ionic_species), 8, 1),
        ("ionic.state_diagnostics", lambda: ionic_state_tp.state_diagnostics(species=ionic_species), 4, 0),
        (
            "ionic.activity_coefficient.mean_ionic_molality.full_rebuild",
            lambda: ePCSAFTMixture.from_params(deepcopy(ctx["ionic_params"]), species=ionic_species)
            .state(T=298.15, x=z_feed, P=1.0e5, phase="liq")
            .activity_coefficient(species=ionic_species, mean_ionic_form=True, basis="molality"),
            3,
            0,
        ),
    ]

    profiled_aliases = {
        "ePCSAFTMixture.from_params": "ctor.from_params",
        "ePCSAFTMixture.from_dataset": "ctor.from_dataset.water",
        "ePCSAFTMixture.clear_runtime_caches": "mixture.clear_runtime_caches",
        "ePCSAFTMixture.reset_runtime_cache_stats": "mixture.reset_runtime_cache_stats",
        "ePCSAFTMixture.runtime_cache_stats": "mixture.runtime_cache_stats",
        "ePCSAFTMixture.equilibrium": "mixture.equilibrium.tp_flash",
        "ePCSAFTMixture.state": "state.from_P",
        "ePCSAFTState.pressure": "state.pressure",
        "ePCSAFTState.density": "state.density",
        "ePCSAFTState.molar_density": "state.molar_density",
        "ePCSAFTState.mass_density": "ionic.mass_density",
        "ePCSAFTState.method_aliases": "state.method_aliases",
        "ePCSAFTState.compressibility_factor": "state.compressibility_factor",
        "ePCSAFTState.residual_helmholtz": "state.residual_helmholtz",
        "ePCSAFTState.temperature_derivative_residual_helmholtz": "state.temperature_derivative_residual_helmholtz",
        "ePCSAFTState.residual_enthalpy": "state.residual_enthalpy",
        "ePCSAFTState.residual_entropy": "state.residual_entropy",
        "ePCSAFTState.residual_gibbs": "state.residual_gibbs",
        "ePCSAFTState.composition_derivative_residual_helmholtz": "state.composition_derivative_residual_helmholtz",
        "ePCSAFTState.residual_chemical_potential": "state.residual_chemical_potential",
        "ePCSAFTState.fugacity_coefficient": "state.fugacity_coefficient",
        "ePCSAFTState.state_diagnostics": "state.state_diagnostics.neutral",
        "ePCSAFTState.relative_permittivity": "ionic.relative_permittivity",
        "ePCSAFTState.osmotic_coefficient": "ionic.osmotic_coefficient",
        "ePCSAFTState.activity_coefficient": "ionic.activity_coefficient.component",
        "ePCSAFTState.solvation_free_energy": "ionic.solvation_free_energy",
    }
    bench_names = {name for name, _, _, _ in benches}
    missing_coverage = [api for api, alias in profiled_aliases.items() if alias not in bench_names]
    assert not missing_coverage, f"Missing benchmark coverage for API methods: {missing_coverage}"
    required_variant_benches = {
        "state.compressibility_factor.terms",
        "state.residual_helmholtz.terms",
        "state.temperature_derivative_residual_helmholtz.terms",
        "state.residual_chemical_potential.terms",
        "state.fugacity_coefficient.terms",
    }
    missing_variant_benches = sorted(required_variant_benches - bench_names)
    assert (
        not missing_variant_benches
    ), "Missing contribution-term benchmark coverage; extend tests/profile/test_runtime_profile.py: " + ", ".join(
        missing_variant_benches
    )

    expected_state_methods = _public_callables(type(neutral_state_tp))
    expected_mixture_methods = _public_callables(type(ctx["neutral_mix"]))
    missing_state_methods = sorted(
        method for method in expected_state_methods if f"ePCSAFTState.{method}" not in profiled_aliases
    )
    missing_mixture_methods = sorted(
        method for method in expected_mixture_methods if f"ePCSAFTMixture.{method}" not in profiled_aliases
    )
    assert (
        not missing_state_methods
    ), "Unprofiled public ePCSAFTState methods found; extend tests/profile/test_runtime_profile.py: " + ", ".join(
        missing_state_methods
    )
    assert (
        not missing_mixture_methods
    ), "Unprofiled public ePCSAFTMixture methods found; extend tests/profile/test_runtime_profile.py: " + ", ".join(
        missing_mixture_methods
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
            + ", ".join(
                row["name"] for row in sorted(init_heavy, key=lambda item: item["first_over_median"], reverse=True)[:8]
            )
        )

    rebuild_row = next(
        (row for row in rows if row["name"] == "ionic.activity_coefficient.mean_ionic_molality.full_rebuild"), None
    )
    reuse_row = next((row for row in rows if row["name"] == "ionic.activity_coefficient.mean_ionic_molality"), None)
    if rebuild_row is not None and reuse_row is not None and reuse_row["mean_ms"] > 0.0:
        ratio = rebuild_row["mean_ms"] / reuse_row["mean_ms"]
        bottleneck_notes.append(_full_rebuild_reuse_note(rebuild_row, reuse_row, ratio))

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
    report_text = REPORT_MD.read_text(encoding="utf-8")
    assert "Full rebuild vs reused-state MIAC ratio" in report_text
    assert "reuse ePCSAFTMixture and ePCSAFTState objects in hot loops" in report_text
