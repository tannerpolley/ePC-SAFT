from __future__ import annotations

import csv
import json
import math
import sys
from pathlib import Path
from typing import Any

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[5]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts.dev.native_runtime_env import apply_to_current_process

ANALYSIS_DIR = REPO_ROOT / "analyses" / "paper_validation" / "native" / "2022_ascani"
SOURCE_CSV = REPO_ROOT / "data" / "reference" / "multiphase" / "ascani_case2_model_comparison.csv"
PROCESSED_DIR = ANALYSIS_DIR / "data" / "processed"
RESULTS_DIR = ANALYSIS_DIR / "results" / "electrolyte_lle"
NORMALIZED_SOURCE_CSV = PROCESSED_DIR / "source_expected_phase_compositions.csv"
FEED_CONVERSION_CSV = PROCESSED_DIR / "feed_conversion_table.csv"
SUMMARY_JSON = RESULTS_DIR / "summary.json"

R_GAS = 8.31446261815324
TEMPERATURE_K = 298.15
PRESSURE_PA = 1.0e5
PRESSURE_BAR = PRESSURE_PA / 1.0e5
TRACE_FLOOR = 1.0e-10

SPECIES = ["H2O", "Butanol", "Na+", "K+", "Cl-"]
PAPER_COMPONENTS = ["H2O", "Butanol", "NaCl", "KCl"]
PSEUDO_TERNARY_COMPONENTS = ["H2O", "Butanol", "total_salt"]
STAGE_STATUS_ACCEPTED = "accepted_public_native_ipopt"

MW_FORMULA_KG_PER_MOL = {
    "H2O": 18.01528e-3,
    "Butanol": 74.1216e-3,
    "NaCl": 58.44277e-3,
    "KCl": 74.5513e-3,
}

FEED_MASS_FRACTIONS = {
    "H2O": 0.8094,
    "Butanol": 0.1728,
    "NaCl": 0.0054,
    "KCl": 0.0124,
}

RESOLVED_TOLERANCES = {
    "material_balance_abs": 1.0e-8,
    "charge_balance_abs": 1.0e-8,
    "neutral_fugacity_abs": 1.0e-7,
    "salt_pair_fugacity_abs": 1.0e-7,
    "density_recompute_rel": 1.0e-8,
    "density_min_mol_m3": 1000.0,
    "phase_distance_min": 0.1,
    "minimum_phase_fraction_min": 1.0e-4,
    "ghat_delta_max": -1.0e-8,
    "tpdf_tolerance": 1.0e-8,
    "trace_floor": TRACE_FLOOR,
}

COMMAND_LIST = [
    "uv run python scripts/dev/doctor.py --require-ipopt",
    "uv run python run_pytest.py tests/native/equilibrium/test_route_builders.py tests/native/equilibrium/test_electrolyte_lle_residual_surface.py tests/native/equilibrium/test_electrolyte_lle_residual_jacobian.py tests/equilibrium/electrolyte -q",
    r"uv run python analyses\paper_validation\native\2022_ascani\scripts\run_all.py",
    "uv run python run_pytest.py tests/workflows/paper_validation/test_ascani_2022_lle_validation.py tests/workflows/benchmarks/test_benchmark_literature_suite.py -q",
    "uv run python scripts/benchmarks/benchmark_literature_suite.py --case ascani_2022_distributed_ion_lle --json build/validation/ascani_2022_lle.json",
]

PAPER_GIBBS = {
    "g_hat_feed_j_per_mol": -27361.317,
    "g_hat_eq_j_per_mol": -27479.860,
    "delta_g_hat_j_per_mol": -118.543,
}


def rel(path: Path) -> str:
    return path.relative_to(REPO_ROOT).as_posix()


def feed_from_case2_mass_fractions() -> np.ndarray:
    formula_moles = {
        component: FEED_MASS_FRACTIONS[component] / MW_FORMULA_KG_PER_MOL[component]
        for component in PAPER_COMPONENTS
    }
    explicit_moles = {
        "H2O": formula_moles["H2O"],
        "Butanol": formula_moles["Butanol"],
        "Na+": formula_moles["NaCl"],
        "K+": formula_moles["KCl"],
        "Cl-": formula_moles["NaCl"] + formula_moles["KCl"],
    }
    values = np.asarray([explicit_moles[label] for label in SPECIES], dtype=float)
    return values / float(np.sum(values))


FEED = feed_from_case2_mass_fractions().tolist()


def feed_conversion_rows() -> list[dict[str, Any]]:
    formula_moles = {
        component: FEED_MASS_FRACTIONS[component] / MW_FORMULA_KG_PER_MOL[component]
        for component in PAPER_COMPONENTS
    }
    formula_total = sum(formula_moles.values())
    explicit = {
        "H2O": formula_moles["H2O"],
        "Butanol": formula_moles["Butanol"],
        "Na+": formula_moles["NaCl"],
        "K+": formula_moles["KCl"],
        "Cl-": formula_moles["NaCl"] + formula_moles["KCl"],
    }
    explicit_total = sum(explicit.values())
    rows = []
    for component in PAPER_COMPONENTS:
        rows.append(
            {
                "basis": "formula",
                "species": component,
                "mass_fraction": FEED_MASS_FRACTIONS[component],
                "moles_on_1kg_feed": formula_moles[component],
                "mole_fraction": formula_moles[component] / formula_total,
            }
        )
    for species in SPECIES:
        rows.append(
            {
                "basis": "explicit_ion",
                "species": species,
                "mass_fraction": "",
                "moles_on_1kg_feed": explicit[species],
                "mole_fraction": explicit[species] / explicit_total,
            }
        )
    return rows


def load_source_rows() -> list[dict[str, str]]:
    with SOURCE_CSV.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_rows(path: Path, fieldnames: list[str], rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


def write_normalized_source(rows: list[dict[str, str]]) -> None:
    mapping = {
        "$x_{water}^{(org)}$": ("org", "H2O"),
        "$x_{butanol}^{(org)}$": ("org", "Butanol"),
        "$x_{NaCl}^{(org)}$": ("org", "NaCl"),
        "$x_{KCl}^{(org)}$": ("org", "KCl"),
        "$x_{water}^{(aq)}$": ("aq", "H2O"),
        "$x_{butanol}^{(aq)}$": ("aq", "Butanol"),
        "$x_{NaCl}^{(aq)}$": ("aq", "NaCl"),
        "$x_{KCl}^{(aq)}$": ("aq", "KCl"),
    }
    phase_rows = []
    for row in rows:
        mapped = mapping.get(row["quantity"])
        if mapped is None:
            continue
        phase, component = mapped
        phase_rows.append(
            {
                "phase": phase,
                "component": component,
                "paper_mole_fraction": row["paper"],
                "model_2020": row["model_2020"],
                "model_2025_num": row["model_2025_num"],
            }
        )
    write_rows(
        NORMALIZED_SOURCE_CSV,
        ["phase", "component", "paper_mole_fraction", "model_2020", "model_2025_num"],
        phase_rows,
    )
    write_rows(
        FEED_CONVERSION_CSV,
        ["basis", "species", "mass_fraction", "moles_on_1kg_feed", "mole_fraction"],
        feed_conversion_rows(),
    )


def _current_result():
    apply_to_current_process()
    import epcsaft
    from epcsaft import ePCSAFTMixture

    mix = ePCSAFTMixture.from_dataset("2022_Ascani", SPECIES, FEED, TEMPERATURE_K)
    options = epcsaft.EquilibriumOptions(max_iterations=500, tolerance=1.0e-8, min_composition=1.0e-12)
    result = mix.equilibrium(kind="electrolyte_lle", T=TEMPERATURE_K, P=PRESSURE_PA, z=FEED, options=options)
    return mix, result, epcsaft.runtime_build_info()["native_dependencies"]["ipopt"]


def solve_payload() -> tuple[bool, dict[str, Any], Any, Any]:
    apply_to_current_process()
    import epcsaft

    runtime_ipopt = epcsaft.runtime_build_info()["native_dependencies"]["ipopt"]
    try:
        mix, result, runtime_ipopt = _current_result()
    except epcsaft.SolutionError as exc:
        diagnostics = dict(getattr(exc, "diagnostics", {}) or {})
        return (
            False,
            {
                "accepted": False,
                "runtime_ipopt": runtime_ipopt,
                "exception_type": type(exc).__name__,
                "exception_message": str(exc),
                "diagnostics": diagnostics,
                "blocker": _blocker_from_diagnostics(diagnostics),
            },
            None,
            None,
        )

    diagnostics = dict(getattr(result, "diagnostics", {}) or {})
    accepted, blocker = _shared_gate_acceptance(diagnostics)
    return (
        accepted,
        {
            "accepted": accepted,
            "runtime_ipopt": runtime_ipopt,
            "solver_backend": diagnostics.get("solver_backend", diagnostics.get("backend", "ipopt")),
            "diagnostics": diagnostics,
            "blocker": blocker,
        },
        mix,
        result,
    )


def _blocker_from_diagnostics(diagnostics: dict[str, Any]) -> dict[str, Any]:
    tpdf = diagnostics.get("tpdf_stability", {})
    if isinstance(tpdf, dict) and tpdf.get("status") == "failed_gate":
        return {"kind": "failed_gate", "reason": tpdf.get("failure_reason", "tpdf_stability")}
    if isinstance(tpdf, dict) and tpdf.get("status") == "blocked_solver":
        return {"kind": "blocked_solver", "reason": "tpdf_stability_solver"}
    route_status = diagnostics.get("route_status")
    solver_status = diagnostics.get("solver_status")
    if route_status or solver_status:
        return {"kind": "blocked_solver", "route_status": route_status, "solver_status": solver_status}
    return {"kind": "blocked_capability", "reason": "no_accepted_public_native_ipopt_result"}


def _shared_gate_acceptance(diagnostics: dict[str, Any]) -> tuple[bool, dict[str, Any]]:
    checks = {
        "solver_backend": diagnostics.get("solver_backend") == "ipopt",
        "derivative_backend": diagnostics.get("derivative_backend") == "cppad_implicit",
        "hessian_approximation": diagnostics.get("hessian_approximation") == "limited-memory",
        "density_backend": diagnostics.get("density_backend") == "liquid_pressure_root",
        "material_balance": float(diagnostics.get("material_balance_norm", math.inf)) <= RESOLVED_TOLERANCES["material_balance_abs"],
        "charge_balance": float(diagnostics.get("charge_balance_norm", math.inf)) <= RESOLVED_TOLERANCES["charge_balance_abs"],
        "neutral_fugacity": float(diagnostics.get("neutral_fugacity_residual_norm", math.inf)) <= RESOLVED_TOLERANCES["neutral_fugacity_abs"],
        "salt_pair_fugacity": float(diagnostics.get("salt_pair_fugacity_residual_norm", math.inf)) <= RESOLVED_TOLERANCES["salt_pair_fugacity_abs"],
        "phase_distance": float(diagnostics.get("phase_distance", -math.inf)) >= RESOLVED_TOLERANCES["phase_distance_min"],
        "minimum_phase_fraction": float(diagnostics.get("minimum_phase_fraction", -math.inf)) >= RESOLVED_TOLERANCES["minimum_phase_fraction_min"],
        "ghat_delta": float(diagnostics.get("gibbs_delta", math.inf)) < RESOLVED_TOLERANCES["ghat_delta_max"],
        "tpdf": bool(diagnostics.get("tpdf_stability", {}).get("accepted", False)),
    }
    density_errors = diagnostics.get("density_recompute_relative_errors", [])
    checks["density_recompute"] = bool(density_errors) and all(
        float(row["relative_error"]) <= RESOLVED_TOLERANCES["density_recompute_rel"] for row in density_errors
    )
    checks["density_min"] = all(
        float(value) >= RESOLVED_TOLERANCES["density_min_mol_m3"]
        for value in diagnostics.get("phase_densities", [])
    )
    failed = [name for name, passed in checks.items() if not passed]
    if failed:
        return False, {"kind": "failed_gate", "failed_checks": failed}
    return True, {}


def current_solution() -> tuple[Any, Any]:
    mix, result, _runtime = _current_result()
    return mix, result


def paper_phase_formula_rows() -> list[dict[str, Any]]:
    source_rows = load_source_rows()
    values: dict[tuple[str, str], float] = {}
    label_map = {
        "$x_{water}^{(org)}$": ("org", "H2O"),
        "$x_{butanol}^{(org)}$": ("org", "Butanol"),
        "$x_{NaCl}^{(org)}$": ("org", "NaCl"),
        "$x_{KCl}^{(org)}$": ("org", "KCl"),
        "$x_{water}^{(aq)}$": ("aq", "H2O"),
        "$x_{butanol}^{(aq)}$": ("aq", "Butanol"),
        "$x_{NaCl}^{(aq)}$": ("aq", "NaCl"),
        "$x_{KCl}^{(aq)}$": ("aq", "KCl"),
    }
    for row in source_rows:
        mapped = label_map.get(row["quantity"])
        if mapped is not None:
            values[mapped] = float(row["paper"])
    out = []
    for phase in ("org", "aq"):
        total = sum(values[(phase, component)] for component in PAPER_COMPONENTS)
        for component in PAPER_COMPONENTS:
            out.append(
                {
                    "source": "Ascani_2022_Case2_paper",
                    "phase": phase,
                    "component": component,
                    "formula_mole_fraction": values[(phase, component)] / total,
                }
            )
    return out


def current_phase_formula_rows(result: Any) -> list[dict[str, Any]]:
    out: list[dict[str, Any]] = []
    for phase in result.phases:
        x = np.asarray(phase.composition, dtype=float)
        formula = {
            "H2O": float(x[0]),
            "Butanol": float(x[1]),
            "NaCl": float(x[2]),
            "KCl": float(x[3]),
        }
        total = sum(formula.values())
        for component in PAPER_COMPONENTS:
            out.append(
                {
                    "source": "current_native_ipopt_full_case2",
                    "phase": phase.label,
                    "component": component,
                    "formula_mole_fraction": formula[component] / total,
                    "phase_fraction": float(phase.phase_fraction),
                    "density_mol_m3": float(phase.density),
                }
            )
    return out


def current_feed_formula_rows() -> list[dict[str, Any]]:
    formula = {
        "H2O": FEED[0],
        "Butanol": FEED[1],
        "NaCl": FEED[2],
        "KCl": FEED[3],
    }
    total = sum(formula.values())
    return [
        {
            "source": "current_native_ipopt_full_case2",
            "phase": "feed",
            "component": component,
            "formula_mole_fraction": formula[component] / total,
        }
        for component in PAPER_COMPONENTS
    ]


def formula_rows_to_phase_map(rows: list[dict[str, Any]]) -> dict[str, dict[str, float]]:
    phase_map: dict[str, dict[str, float]] = {}
    for row in rows:
        phase_map.setdefault(str(row["phase"]), {})[str(row["component"])] = float(row["formula_mole_fraction"])
    return phase_map


def formula_to_mass_fraction(formula: dict[str, float]) -> dict[str, float]:
    masses = {
        component: float(formula.get(component, 0.0)) * MW_FORMULA_KG_PER_MOL[component]
        for component in PAPER_COMPONENTS
    }
    total = sum(masses.values())
    if total <= 0.0:
        raise ValueError("formula composition has no positive mass.")
    return {component: value / total for component, value in masses.items()}


def pseudo_ternary_mass_fraction(formula: dict[str, float]) -> dict[str, float]:
    mass = formula_to_mass_fraction(formula)
    return {
        "H2O": mass["H2O"],
        "Butanol": mass["Butanol"],
        "total_salt": mass["NaCl"] + mass["KCl"],
    }


def ternary_xy(pseudo: dict[str, float]) -> tuple[float, float]:
    butanol = float(pseudo["Butanol"])
    salt = float(pseudo["total_salt"])
    return butanol + 0.5 * salt, (math.sqrt(3.0) / 2.0) * salt


def phase_diagram_rows(result: Any | None = None) -> list[dict[str, Any]]:
    if result is None:
        _mix, result = current_solution()
    current_rows = current_phase_formula_rows(result)
    current_feed = current_feed_formula_rows()
    paper_rows = paper_phase_formula_rows()
    rows: list[dict[str, Any]] = []
    for series, source_rows in (
        ("paper_case2", paper_rows),
        ("current_ipopt", current_rows),
        ("current_feed", current_feed),
    ):
        phase_map = formula_rows_to_phase_map(source_rows)
        for phase, formula in phase_map.items():
            pseudo = pseudo_ternary_mass_fraction(formula)
            x_plot, y_plot = ternary_xy(pseudo)
            rows.append(
                {
                    "series": series,
                    "phase": phase,
                    "w_water": pseudo["H2O"],
                    "w_butanol": pseudo["Butanol"],
                    "w_total_salt": pseudo["total_salt"],
                    "x_plot": x_plot,
                    "y_plot": y_plot,
                }
            )
    return rows


def _state_ln_fugacity_bar(mix: Any, composition: np.ndarray, rho_guess: float | None = None) -> tuple[np.ndarray, float]:
    state = mix.state(TEMPERATURE_K, composition, P=PRESSURE_PA, phase="liq", rho_guess=rho_guess)
    ln_phi = np.asarray(state.fugacity_coefficient(natural_log=True), dtype=float)
    ln_f = np.log(np.maximum(np.asarray(composition, dtype=float), 1.0e-300)) + ln_phi + math.log(PRESSURE_BAR)
    return ln_f, float(state.molar_density())


def fugacity_comparison_rows(mix: Any | None = None, result: Any | None = None) -> list[dict[str, Any]]:
    if mix is None or result is None:
        mix, result = current_solution()
    source_rows = {row["quantity"]: row for row in load_source_rows()}
    paper = {
        "H2O": float(source_rows["$\\ln(f_{water}/bar)$"]["paper"]),
        "Butanol": float(source_rows["$\\ln(f_{butanol}/bar)$"]["paper"]),
        "NaCl": float(source_rows["$\\ln(f_{\\pm,NaCl}/bar)$"]["paper"]),
        "KCl": float(source_rows["$\\ln(f_{\\pm,KCl}/bar)$"]["paper"]),
    }
    current_by_phase: dict[str, dict[str, float]] = {}
    for phase in result.phases:
        ln_f, density = _state_ln_fugacity_bar(mix, np.asarray(phase.composition, dtype=float), float(phase.density))
        current_by_phase[phase.label] = {
            "H2O": float(ln_f[0]),
            "Butanol": float(ln_f[1]),
            "NaCl": 0.5 * (float(ln_f[2]) + float(ln_f[4])),
            "KCl": 0.5 * (float(ln_f[3]) + float(ln_f[4])),
            "density_mol_m3": density,
        }
    rows: list[dict[str, Any]] = []
    for component in PAPER_COMPONENTS:
        for phase in ("org", "aq"):
            current = current_by_phase.get(phase, {}).get(component)
            rows.append(
                {
                    "quantity": f"ln_f_{component}_bar",
                    "component": component,
                    "phase": phase,
                    "paper_ln_f_bar": paper[component],
                    "current_ln_f_bar": "" if current is None else current,
                    "current_minus_paper": "" if current is None else current - paper[component],
                    "current_basis": f"mean_ionic_{component}" if component in {"NaCl", "KCl"} else "component",
                    "note": "",
                }
            )
    return rows


def gibbs_rows(mix: Any | None = None, result: Any | None = None) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    if mix is None or result is None:
        mix, result = current_solution()
    g_feed = float(result.diagnostics["gibbs_feed"]) * R_GAS * TEMPERATURE_K
    g_eq = float(result.diagnostics["gibbs_split"]) * R_GAS * TEMPERATURE_K
    g_delta = g_eq - g_feed
    comparison = [
        {
            "quantity": "g_hat_feed_j_per_mol",
            "paper": PAPER_GIBBS["g_hat_feed_j_per_mol"],
            "current_native_ln_fugacity_basis_j_per_mol": g_feed,
            "current_minus_paper": g_feed - PAPER_GIBBS["g_hat_feed_j_per_mol"],
        },
        {
            "quantity": "g_hat_eq_j_per_mol",
            "paper": PAPER_GIBBS["g_hat_eq_j_per_mol"],
            "current_native_ln_fugacity_basis_j_per_mol": g_eq,
            "current_minus_paper": g_eq - PAPER_GIBBS["g_hat_eq_j_per_mol"],
        },
        {
            "quantity": "delta_g_hat_j_per_mol",
            "paper": PAPER_GIBBS["delta_g_hat_j_per_mol"],
            "current_native_ln_fugacity_basis_j_per_mol": g_delta,
            "current_minus_paper": g_delta - PAPER_GIBBS["delta_g_hat_j_per_mol"],
        },
    ]
    phases = [
        {
            "phase": phase.label,
            "phase_fraction": float(phase.phase_fraction),
            "density_mol_m3": float(phase.density),
        }
        for phase in result.phases
    ]
    return comparison, phases


def stage_c_seed_payloads(result: Any) -> list[dict[str, Any]]:
    feed = np.asarray(FEED, dtype=float)
    phases = {phase.label: np.asarray(phase.composition, dtype=float) for phase in result.phases}
    aq = phases.get("aq", feed)
    org = phases.get("org", feed)
    cation_ratio = feed[2:4] / float(np.sum(feed[2:4]))
    stage_a_mapped = np.asarray([0.37006036297653244, 0.6214918558517049, 0.00844778117176257 * cation_ratio[0], 0.00844778117176257 * cation_ratio[1], 0.00844778117176257])
    seeds = [
        ("source_expected_or_table_seed", aq, org, "source_table"),
        ("water_rich_salt_rich_seed", np.asarray([0.98, 0.01, 0.002, 0.0035, 0.0055]), org, "deterministic"),
        ("butanol_rich_trace_salt_seed", aq, np.asarray([0.05, 0.949, 0.0002, 0.0003, 0.0005]), "deterministic"),
        ("balanced_feed_perturbation_seed", feed * np.asarray([1.02, 0.98, 1.0, 1.0, 1.0]), feed, "deterministic"),
        ("stage_a_accepted_mapped_seed", aq, stage_a_mapped, "stage_a_mapped_by_feed_cation_ratio"),
    ]
    out = []
    for name, aqueous, organic, source in seeds:
        aqueous = _normalize_charge_neutral_seed(aqueous)
        organic = _normalize_charge_neutral_seed(organic)
        out.append(
            {
                "name": name,
                "source": source,
                "aqueous_phase_composition": aqueous.tolist(),
                "organic_phase_composition": organic.tolist(),
            }
        )
    return out


def _normalize_charge_neutral_seed(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float).copy()
    values = np.maximum(values, 1.0e-12)
    values[4] = values[2] + values[3]
    return values / float(np.sum(values))


def summary_payload(accepted: bool, solve: dict[str, Any], result: Any | None = None) -> dict[str, Any]:
    diagnostics = dict(solve.get("diagnostics", {}) or {})
    status = STAGE_STATUS_ACCEPTED if accepted else _strict_status_from_blocker(solve.get("blocker", {}))
    retained_outputs = [
        rel(SUMMARY_JSON),
        rel(FEED_CONVERSION_CSV),
        rel(NORMALIZED_SOURCE_CSV),
        "analyses/paper_validation/native/2022_ascani/results/electrolyte_lle/phase_split.csv",
        "analyses/paper_validation/native/2022_ascani/figures/stability_summary/stability_summary.csv",
        "analyses/paper_validation/native/2022_ascani/figures/density_summary/density_summary.csv",
        "analyses/paper_validation/native/2022_ascani/figures/residual_summary/residual_summary.csv",
    ]
    return {
        "schema_version": 1,
        "stage": "A-C",
        "lane_id": "ascani_2022_distributed_ion_lle",
        "status": status,
        "status_reason": "all shared native Ipopt gates and hard TPDF certification passed" if accepted else "strict gate failed",
        "source_assets": [
            rel(SOURCE_CSV),
            "docs/papers/md/Ascani, Sadowski, Held - 2022 - Calculation of Multiphase Equilibria Containing Mixed Solvents and M.md",
            "docs/papers/md/Ascani, Sadowski, Held - 2022 - Supporting Information for Calculation of multiphase equilibria containing mixed solvents and mixed..md",
        ],
        "command_list": COMMAND_LIST,
        "resolved_tolerances": dict(RESOLVED_TOLERANCES),
        "feed": {
            "species": SPECIES,
            "mass_fractions": dict(FEED_MASS_FRACTIONS),
            "mole_fractions": FEED,
            "temperature_K": TEMPERATURE_K,
            "pressure_Pa": PRESSURE_PA,
        },
        "attempt_matrix": [
            {"basis": "paper_era_no_ssm_ds_born", "attempted": True, "accepted": accepted},
            {"basis": "ssm_ds_born_comparison", "attempted": False, "accepted": False, "reason": "primary basis accepted" if accepted else "not_started"},
        ],
        "deterministic_seed_payloads": [] if result is None else stage_c_seed_payloads(result),
        "ipopt_runtime_diagnostics": {
            "runtime_ipopt": solve.get("runtime_ipopt"),
            "solver_backend": solve.get("solver_backend"),
            "solver_status": diagnostics.get("solver_status"),
            "application_status": diagnostics.get("application_status"),
            "last_callback_exception": diagnostics.get("last_callback_exception"),
        },
        "derivative_diagnostics": {
            "derivative_backend": diagnostics.get("derivative_backend"),
            "gradient_approximation": diagnostics.get("gradient_approximation"),
            "jacobian_approximation": diagnostics.get("jacobian_approximation"),
            "exact_gradient_required": diagnostics.get("exact_gradient_required"),
            "exact_jacobian_required": diagnostics.get("exact_jacobian_required"),
        },
        "hessian_approximation_diagnostics": {
            "hessian_approximation": diagnostics.get("hessian_approximation"),
            "exact_hessian_required": False,
        },
        "density_diagnostics": {
            "density_backend": diagnostics.get("density_backend"),
            "phase_densities": diagnostics.get("phase_densities"),
            "density_recompute_relative_errors": diagnostics.get("density_recompute_relative_errors"),
        },
        "material_charge_fugacity_residuals": {
            "material_balance_norm": diagnostics.get("material_balance_norm"),
            "charge_balance_norm": diagnostics.get("charge_balance_norm"),
            "neutral_fugacity_residual_norm": diagnostics.get("neutral_fugacity_residual_norm"),
            "salt_pair_fugacity_residual_norm": diagnostics.get("salt_pair_fugacity_residual_norm"),
            "neutral_log_fugacity_residuals_raw": diagnostics.get("neutral_log_fugacity_residuals_raw"),
            "salt_pair_log_fugacity_residuals_raw": diagnostics.get("salt_pair_log_fugacity_residuals_raw"),
        },
        "tpdf_stability_results": diagnostics.get("tpdf_stability"),
        "ghat_feed": diagnostics.get("gibbs_feed"),
        "ghat_split": diagnostics.get("gibbs_split"),
        "ghat_delta": diagnostics.get("gibbs_delta"),
        "retained_outputs": retained_outputs,
        "blockers": [] if accepted else [solve.get("blocker", {})],
        "claim_boundary": {
            "accepted_route": accepted,
            "paper_match_claim": "not_claimed",
            "note": "This proves the public native Ipopt liquid-root route for Ascani 2022 Case Study 2; exact paper matching remains a separate source-backed comparison claim.",
        },
        "expected": {
            "status": STAGE_STATUS_ACCEPTED,
            "accepted": True,
            "solver_backend": "ipopt",
            **RESOLVED_TOLERANCES,
        },
        "solve": solve,
    }


def _strict_status_from_blocker(blocker: dict[str, Any]) -> str:
    kind = str(blocker.get("kind", ""))
    if kind in {"blocked_source_data", "blocked_solver", "blocked_capability", "failed_gate"}:
        return kind
    return "failed_gate"


def write_retained_outputs(summary: dict[str, Any], mix: Any | None, result: Any | None) -> None:
    if result is None:
        return
    phase_rows = current_phase_formula_rows(result)
    write_rows(
        RESULTS_DIR / "phase_split.csv",
        ["source", "phase", "component", "formula_mole_fraction", "phase_fraction", "density_mol_m3"],
        phase_rows,
    )
    write_rows(
        ANALYSIS_DIR / "figures" / "figure_4" / "output" / "data" / "current_phase_compositions.csv",
        ["source", "phase", "component", "formula_mole_fraction", "phase_fraction", "density_mol_m3"],
        phase_rows,
    )
    fugacity_rows = fugacity_comparison_rows(mix, result)
    write_rows(
        ANALYSIS_DIR / "figures" / "table_5_fugacity" / "output" / "data" / "table_5_fugacity_comparison.csv",
        ["quantity", "component", "phase", "paper_ln_f_bar", "current_ln_f_bar", "current_minus_paper", "current_basis", "note"],
        fugacity_rows,
    )
    diagnostics = summary["solve"]["diagnostics"]
    stability_rows = []
    for row in diagnostics.get("tpdf_stability", {}).get("trials", []):
        stability_rows.append(
            {
                "parent_phase": row.get("parent_phase"),
                "seed_name": row.get("seed_name"),
                "status": row.get("status"),
                "min_tpd": row.get("min_tpd"),
                "tolerance": row.get("tolerance"),
                "failure_reason": row.get("failure_reason"),
            }
        )
    write_rows(
        ANALYSIS_DIR / "figures" / "stability_summary" / "stability_summary.csv",
        ["parent_phase", "seed_name", "status", "min_tpd", "tolerance", "failure_reason"],
        stability_rows,
    )
    density_rows = diagnostics.get("density_recompute_relative_errors", [])
    write_rows(
        ANALYSIS_DIR / "figures" / "density_summary" / "density_summary.csv",
        ["phase", "reported_density_mol_m3", "recomputed_density_mol_m3", "relative_error"],
        density_rows,
    )
    residual = summary["material_charge_fugacity_residuals"]
    write_rows(
        ANALYSIS_DIR / "figures" / "residual_summary" / "residual_summary.csv",
        ["name", "value"],
        [{"name": key, "value": value} for key, value in residual.items() if not isinstance(value, dict)],
    )
