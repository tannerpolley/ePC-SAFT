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
SUMMARY_JSON = RESULTS_DIR / "summary.json"

R_GAS = 8.31446261815324
TEMPERATURE_K = 298.15
PRESSURE_PA = 1.0e5
PRESSURE_BAR = PRESSURE_PA / 1.0e5
MIN_PHASE_DISTANCE = 0.1

SPECIES = ["H2O", "Butanol", "Na+", "Cl-"]
PAPER_COMPONENTS = ["H2O", "Butanol", "NaCl", "KCl"]
PSEUDO_TERNARY_COMPONENTS = ["H2O", "Butanol", "total_salt"]

MW_FORMULA_KG_PER_MOL = {
    "H2O": 18.01528e-3,
    "Butanol": 74.1216e-3,
    "NaCl": 58.44277e-3,
    "KCl": 74.5513e-3,
}

# Source-like NaCl case used for the executable public-API Ipopt gate. The
# paper Case 2 values remain reference data, but this feed is intentionally not
# forced to close exactly to the mixed NaCl/KCl paper split.
SOURCE_LIKE_AQ_PHASE = [0.798324680201737, 0.016320352824141723, 0.09267748348706063, 0.09267748348706063]
SOURCE_LIKE_ORG_PHASE = [0.37006036048879404, 0.6214918588210971, 0.004223890345054407, 0.004223890345054407]
SOURCE_LIKE_ORG_FRACTION = 0.613766575013417
FEED = [
    (1.0 - SOURCE_LIKE_ORG_FRACTION) * aq + SOURCE_LIKE_ORG_FRACTION * org
    for aq, org in zip(SOURCE_LIKE_AQ_PHASE, SOURCE_LIKE_ORG_PHASE, strict=True)
]

PAPER_GIBBS = {
    "g_hat_feed_j_per_mol": -27361.317,
    "g_hat_eq_j_per_mol": -27479.860,
    "delta_g_hat_j_per_mol": -118.543,
}


def rel(path: Path) -> str:
    return path.relative_to(REPO_ROOT).as_posix()


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
    phase_rows: list[dict[str, str]] = []
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
                "blocker": {
                    "kind": "native_ipopt_solver_rejected",
                    "route_status": diagnostics.get("route_status"),
                    "solver_status": diagnostics.get("solver_status"),
                },
            },
            None,
            None,
        )
    diagnostics = dict(getattr(result, "diagnostics", {}) or {})
    phase_distance = float(diagnostics.get("phase_distance", 0.0))
    if phase_distance < MIN_PHASE_DISTANCE:
        return (
            False,
            {
                "accepted": False,
                "runtime_ipopt": runtime_ipopt,
                "diagnostics": diagnostics,
                "blocker": {
                    "kind": "native_ipopt_phase_split_too_small",
                    "phase_distance": phase_distance,
                    "minimum_phase_distance": MIN_PHASE_DISTANCE,
                },
            },
            mix,
            result,
        )
    return (
        True,
        {
            "accepted": bool(diagnostics.get("accepted", True)),
            "runtime_ipopt": runtime_ipopt,
            "solver_backend": diagnostics.get("solver_backend", diagnostics.get("backend", "ipopt")),
            "diagnostics": diagnostics,
        },
        mix,
        result,
    )


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
        salt = 0.5 * (float(x[2]) + float(x[3]))
        formula = {
            "H2O": float(x[0]),
            "Butanol": float(x[1]),
            "NaCl": salt,
            "KCl": 0.0,
        }
        total = sum(formula.values())
        for component in PAPER_COMPONENTS:
            out.append(
                {
                    "source": "current_native_ipopt_source_like_nacl",
                    "phase": phase.label,
                    "component": component,
                    "formula_mole_fraction": formula[component] / total,
                    "explicit_mole_fraction": float(x[SPECIES.index(component)]) if component in SPECIES else "",
                    "phase_fraction": float(phase.phase_fraction),
                    "density_mol_m3": float(phase.density),
                }
            )
    return out


def current_feed_formula_rows() -> list[dict[str, Any]]:
    x = np.asarray(FEED, dtype=float)
    salt = 0.5 * (float(x[2]) + float(x[3]))
    formula = {
        "H2O": float(x[0]),
        "Butanol": float(x[1]),
        "NaCl": salt,
        "KCl": 0.0,
    }
    total = sum(formula.values())
    return [
        {
            "source": "current_native_ipopt_source_like_nacl",
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
    masses = {component: float(formula.get(component, 0.0)) * MW_FORMULA_KG_PER_MOL[component] for component in PAPER_COMPONENTS}
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


def phase_diagram_rows() -> list[dict[str, Any]]:
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


def fugacity_comparison_rows() -> list[dict[str, Any]]:
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
            "NaCl": 0.5 * (float(ln_f[2]) + float(ln_f[3])),
            "density_mol_m3": density,
        }
    rows: list[dict[str, Any]] = []
    for component in ("H2O", "Butanol", "NaCl", "KCl"):
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
                    "current_basis": "mean_ionic_NaCl" if component == "NaCl" else "component",
                    "note": "current accepted source-like gate has no K+ species" if component == "KCl" else "",
                }
            )
    return rows


def _phase_g_hat_j_per_mol(mix: Any, composition: np.ndarray, rho_guess: float | None = None) -> float:
    state = mix.state(TEMPERATURE_K, composition, P=PRESSURE_PA, phase="liq", rho_guess=rho_guess)
    rho = float(state.molar_density())
    ideal = float(np.sum(composition * (np.log(np.maximum(rho * composition, 1.0e-300)) - 1.0)))
    residual = float(state.residual_helmholtz())
    pressure_work = PRESSURE_PA / (rho * R_GAS * TEMPERATURE_K)
    return (ideal + residual + pressure_work) * R_GAS * TEMPERATURE_K


def gibbs_rows() -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    mix, result = current_solution()
    phase_volumes = np.asarray(result.diagnostics.get("phase_volumes", []), dtype=float)
    feed_density_guess = None
    if phase_volumes.size and np.all(np.isfinite(phase_volumes)) and float(np.sum(phase_volumes)) > 0.0:
        feed_density_guess = 1.0 / float(np.sum(phase_volumes))
    g_feed = _phase_g_hat_j_per_mol(mix, np.asarray(FEED, dtype=float), feed_density_guess)
    g_eq = float(result.diagnostics["objective"]) * R_GAS * TEMPERATURE_K
    g_delta = g_eq - g_feed
    comparison = [
        {
            "quantity": "g_hat_feed_j_per_mol",
            "paper": PAPER_GIBBS["g_hat_feed_j_per_mol"],
            "current_native_objective_basis": g_feed,
            "current_minus_paper": g_feed - PAPER_GIBBS["g_hat_feed_j_per_mol"],
        },
        {
            "quantity": "g_hat_eq_j_per_mol",
            "paper": PAPER_GIBBS["g_hat_eq_j_per_mol"],
            "current_native_objective_basis": g_eq,
            "current_minus_paper": g_eq - PAPER_GIBBS["g_hat_eq_j_per_mol"],
        },
        {
            "quantity": "delta_g_hat_j_per_mol",
            "paper": PAPER_GIBBS["delta_g_hat_j_per_mol"],
            "current_native_objective_basis": g_delta,
            "current_minus_paper": g_delta - PAPER_GIBBS["delta_g_hat_j_per_mol"],
        },
    ]
    phases = []
    for phase in result.phases:
        phases.append(
            {
                "phase": phase.label,
                "phase_fraction": float(phase.phase_fraction),
                "g_hat_phase_j_per_mol": _phase_g_hat_j_per_mol(
                    mix, np.asarray(phase.composition, dtype=float), float(phase.density)
                ),
                "density_mol_m3": float(phase.density),
            }
        )
    return comparison, phases


def summary_payload(accepted: bool, solve: dict[str, Any]) -> dict[str, Any]:
    return {
        "status": "accepted" if accepted else "blocked",
        "lane": "ascani_2022_distributed_ion_lle",
        "source_records": [rel(SOURCE_CSV), rel(NORMALIZED_SOURCE_CSV)],
        "feed": {"species": SPECIES, "mole_fractions": FEED, "temperature_K": TEMPERATURE_K, "pressure_Pa": PRESSURE_PA},
        "expected": {
            "accepted": True,
            "solver_backend": "ipopt",
            "material_balance_abs": 1.0e-8,
            "charge_balance_abs": 1.0e-8,
            "phase_distance_min": MIN_PHASE_DISTANCE,
        },
        "solve": solve,
    }
