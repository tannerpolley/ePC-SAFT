from __future__ import annotations

import csv
import json
import math
import sys
from pathlib import Path
from typing import Any

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[4]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts.dev.native_runtime_env import apply_to_current_process

apply_to_current_process()

import epcsaft
from epcsaft.equilibrium_core.native_requests import build_reactive_speciation_native_request

SCRIPT_DIR = Path(__file__).resolve().parent
ANALYSIS_DIR = SCRIPT_DIR.parent
INPUT_DIR = ANALYSIS_DIR / "data" / "input"
PROCESSED_DIR = ANALYSIS_DIR / "data" / "processed"
RESULTS_DIR = ANALYSIS_DIR / "results" / "pressure_speciation"
SUMMARY_JSON = RESULTS_DIR / "summary.json"

DATASET_DIR = INPUT_DIR / "MEA_CO2_H2O_phase2"
JOU_VLE_CSV = INPUT_DIR / "Jou_1995_VLE.csv"
PROBLEM_JSON = INPUT_DIR / "phase2_activity_speciation_problem.json"
REACTION_CANDIDATES_CSV = INPUT_DIR / "phase2_activity_constant_candidates.csv"
REACTION_BASIS_CSV = INPUT_DIR / "phase2_reaction_constant_basis.csv"
SOURCE_VERIFICATION_CSV = INPUT_DIR / "phase2_reaction_constant_source_verification.csv"
SOURCE_MANIFEST_CSV = INPUT_DIR / "phase2_parameter_artifact_manifest.csv"

SPECIES = ["CO2", "MEA", "H2O", "MEAH+", "MEACOO-", "HCO3-", "CO3^2-", "H3O+", "OH-"]
MW = {"CO2": 0.04401, "MEA": 0.06108, "H2O": 0.01801528}
ELEMENTS = {
    "C": {"CO2": 1.0, "MEA": 2.0, "MEAH+": 2.0, "MEACOO-": 3.0, "HCO3-": 1.0, "CO3^2-": 1.0},
    "N": {"MEA": 1.0, "MEAH+": 1.0, "MEACOO-": 1.0},
    "H": {"MEA": 7.0, "H2O": 2.0, "MEAH+": 8.0, "MEACOO-": 6.0, "HCO3-": 1.0, "H3O+": 3.0, "OH-": 1.0},
    "O": {
        "CO2": 2.0,
        "MEA": 1.0,
        "H2O": 1.0,
        "MEAH+": 1.0,
        "MEACOO-": 3.0,
        "HCO3-": 3.0,
        "CO3^2-": 3.0,
        "H3O+": 1.0,
        "OH-": 1.0,
    },
}
REACTION_STOICHIOMETRY = {
    "R1": {"H2O": -2.0, "H3O+": 1.0, "OH-": 1.0},
    "R2": {"CO2": -1.0, "H2O": -2.0, "HCO3-": 1.0, "H3O+": 1.0},
    "R3": {"HCO3-": -1.0, "H2O": -1.0, "CO3^2-": 1.0, "H3O+": 1.0},
    "R4": {"MEACOO-": -1.0, "H2O": -1.0, "HCO3-": 1.0, "MEA": 1.0},
    "R5": {"MEAH+": -1.0, "H2O": -1.0, "MEA": 1.0, "H3O+": 1.0},
}


def _rel(path: Path) -> str:
    return path.relative_to(REPO_ROOT).as_posix()


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(_json_like(payload), indent=2) + "\n", encoding="utf-8")


def _json_like(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(k): _json_like(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_like(v) for v in value]
    if isinstance(value, np.ndarray):
        return [_json_like(v) for v in value.tolist()]
    if isinstance(value, (np.bool_, bool)):
        return bool(value)
    if isinstance(value, (np.floating, float)):
        if math.isfinite(float(value)):
            return float(value)
        return str(value)
    if isinstance(value, (np.integer, int)):
        return int(value)
    return value


def _normalized_user_options() -> dict[str, Any]:
    source = json.loads((DATASET_DIR / "user_options.json").read_text(encoding="utf-8"))
    normalized = json.loads(json.dumps(source))
    rel_perm = normalized.get("elec_model", {}).get("rel_perm", {})
    born = normalized.get("elec_model", {}).get("born_model", {})
    if rel_perm.get("differential_mode") == "numerical":
        rel_perm["differential_mode"] = "auto"
    mu_born = born.get("mu_born_model", {})
    if mu_born.get("differential_mode") == "numerical":
        mu_born["differential_mode"] = "auto"
    _write_json(
        PROCESSED_DIR / "normalized_user_options.json",
        {
            "source_user_options": source,
            "runtime_user_options": normalized,
            "normalization_policy": (
                "The copied Phase 2 artifact records differential_mode='numerical'. The current package "
                "accepts exact modes only, so the analysis uses 'auto' to request the package-owned exact route."
            ),
        },
    )
    return normalized


def _selected_jou_row() -> dict[str, str]:
    rows = _read_csv(JOU_VLE_CSV)
    if not rows:
        raise ValueError(f"{_rel(JOU_VLE_CSV)} contains no rows.")
    selected = next(
        (
            row
            for row in rows
            if row.get("MEA_weight_fraction") == "0.3"
            and row.get("temperature") == "25"
            and row.get("total_pressure") == "200"
        ),
        rows[0],
    )
    _write_json(
        PROCESSED_DIR / "selected_jou_row.json",
        {
            "selection_policy": "First 30 wt% MEA, 25 C, 200 kPa Jou 1995 row in the copied table.",
            "row": selected,
            "source_path": _rel(JOU_VLE_CSV),
        },
    )
    return selected


def _apparent_totals(row: dict[str, str]) -> dict[str, float]:
    mea_weight_fraction = float(row["MEA_weight_fraction"])
    co2_loading = float(row["CO2_loading"])
    n_mea = 1.0
    mass_mea = n_mea * MW["MEA"]
    mass_water = mass_mea * (1.0 - mea_weight_fraction) / mea_weight_fraction
    n_water = mass_water / MW["H2O"]
    n_co2 = co2_loading * n_mea
    return {
        "C": 2.0 * n_mea + n_co2,
        "N": n_mea,
        "H": 7.0 * n_mea + 2.0 * n_water,
        "O": n_mea + n_water + 2.0 * n_co2,
    }


def _initial_composition(row: dict[str, str]) -> list[float]:
    totals = _apparent_totals(row)
    loading = float(row["CO2_loading"])
    carbon_from_co2 = max(loading, 1.0e-12)
    meacoo = min(0.65 * carbon_from_co2, 0.45)
    hco3 = min(0.05 * carbon_from_co2, max(carbon_from_co2 - meacoo, 0.0))
    co3 = 1.0e-10
    co2 = max(carbon_from_co2 - meacoo - hco3 - co3, 1.0e-10)
    h3o = 1.0e-10
    oh = 1.0e-10
    meah = min(meacoo + hco3 + 2.0 * co3 + oh - h3o, 0.45)
    mea = max(1.0 - meacoo - meah, 1.0e-10)
    h2o = max((totals["H"] - 7.0 * mea - 8.0 * meah - 6.0 * meacoo - hco3 - 3.0 * h3o - oh) / 2.0, 1.0e-8)
    amounts = {
        "CO2": co2,
        "MEA": mea,
        "H2O": h2o,
        "MEAH+": meah,
        "MEACOO-": meacoo,
        "HCO3-": hco3,
        "CO3^2-": co3,
        "H3O+": h3o,
        "OH-": oh,
    }
    total = sum(amounts.values())
    return [amounts[label] / total for label in SPECIES]


def _balances() -> dict[str, dict[str, float]]:
    return {element: {label: coeff for label, coeff in species.items()} for element, species in ELEMENTS.items()}


def _reaction_log_k(row: dict[str, str], temperature_k: float) -> float:
    return (
        float(row["c1"])
        + float(row["c2"]) / temperature_k
        + float(row["c3"]) * math.log(temperature_k)
        + float(row["c4"]) * temperature_k
    )


def _reactions(temperature_k: float) -> list[epcsaft.ReactionDefinition]:
    rows = _read_csv(REACTION_CANDIDATES_CSV)
    by_id = {row["reaction_id"]: row for row in rows}
    reactions: list[epcsaft.ReactionDefinition] = []
    retained: list[dict[str, Any]] = []
    for reaction_id in ("R1", "R2", "R3", "R4", "R5"):
        row = by_id[reaction_id]
        log_k = _reaction_log_k(row, temperature_k)
        reactions.append(
            epcsaft.ReactionDefinition.from_literature_constant(
                REACTION_STOICHIOMETRY[reaction_id],
                log_equilibrium_constant=log_k,
                name=reaction_id,
                standard_state="mole_fraction_activity",
                source=row["source_files"],
                metadata={
                    "candidate_source": row["candidate_source"],
                    "source_basis": row["source_basis"],
                    "validation_role": row["validation_role"],
                },
            )
        )
        retained.append(
            {
                "reaction_id": reaction_id,
                "log_equilibrium_constant": log_k,
                "stoichiometry": REACTION_STOICHIOMETRY[reaction_id],
                "source_row": row,
            }
        )
    _write_json(
        PROCESSED_DIR / "reaction_constants.json",
        {
            "temperature_K": temperature_k,
            "standard_state": "mole_fraction_activity",
            "reactions": retained,
            "source_path": _rel(REACTION_CANDIDATES_CSV),
        },
    )
    return reactions


def _attempt_public_solve(row: dict[str, str], user_options: dict[str, Any]) -> dict[str, Any]:
    temperature_k = float(row["temperature"]) + 273.15
    pressure_pa = float(row["total_pressure"]) * 1000.0
    initial_x = _initial_composition(row)
    totals = _apparent_totals(row)
    reactions = _reactions(temperature_k)

    def mixture_factory(x: Any, T: float, P: float) -> epcsaft.ePCSAFTMixture:
        _ = P
        return epcsaft.ePCSAFTMixture.from_dataset(DATASET_DIR, SPECIES, np.asarray(x, dtype=float), T, user_options=user_options)

    result = epcsaft.solve_reactive_speciation(
        species=SPECIES,
        mixture_factory=mixture_factory,
        T=temperature_k,
        P=pressure_pa,
        balances=_balances(),
        totals=totals,
        reactions=reactions,
        initial_x=initial_x,
        options=epcsaft.ReactiveSpeciationOptions(
            solver_backend="ipopt",
            jacobian_backend="cppad",
            tolerance=1.0e-8,
            max_iterations=100,
            activity_output="always",
        ),
    )
    composition = np.asarray([result.x[label] for label in SPECIES], dtype=float)
    mix = mixture_factory(composition, temperature_k, pressure_pa)
    state = mix.state(T=temperature_k, P=pressure_pa, x=composition, phase="liq")
    ln_phi = np.asarray(state.fugacity_coefficient(natural_log=True), dtype=float)
    pco2_kpa = composition[SPECIES.index("CO2")] * math.exp(float(ln_phi[SPECIES.index("CO2")])) * pressure_pa / 1000.0
    target_pco2_kpa = float(row["CO2_pressure"])
    return {
        "accepted_native_ipopt_speciation": bool(result.success),
        "composition": result.x,
        "activity_coefficients": result.activity_coefficients,
        "diagnostics": result.diagnostics,
        "ln_fugacity_coefficients": {label: float(value) for label, value in zip(SPECIES, ln_phi)},
        "density_mol_m3": float(state.molar_density()),
        "pressure_comparison": {
            "source_CO2_pressure_kPa": target_pco2_kpa,
            "model_CO2_pressure_kPa": pco2_kpa,
            "log10_error": math.log10(max(pco2_kpa, 1.0e-300) / max(target_pco2_kpa, 1.0e-300)),
            "pressure_model": "liquid_fugacity_with_ideal_vapor_side",
        },
    }


def _diagnose_native_residual(row: dict[str, str], user_options: dict[str, Any]) -> dict[str, Any]:
    temperature_k = float(row["temperature"]) + 273.15
    pressure_pa = float(row["total_pressure"]) * 1000.0
    initial_x = np.asarray(_initial_composition(row), dtype=float)
    reactions = _reactions(temperature_k)
    mix = epcsaft.ePCSAFTMixture.from_dataset(DATASET_DIR, SPECIES, initial_x, temperature_k, user_options=user_options)
    request = build_reactive_speciation_native_request(
        T=temperature_k,
        P=pressure_pa,
        initial_x=initial_x,
        balance_matrix=np.asarray([[ELEMENTS[element].get(label, 0.0) for label in SPECIES] for element in ELEMENTS], dtype=float),
        total_vector=np.asarray([_apparent_totals(row)[element] for element in ELEMENTS], dtype=float),
        species=SPECIES,
        reactions=reactions,
        options=epcsaft.ReactiveSpeciationOptions(
            solver_backend="ipopt",
            jacobian_backend="cppad",
            tolerance=1.0e-8,
            max_iterations=100,
            activity_output="always",
        ),
    )
    try:
        payload = epcsaft._core._evaluate_chemical_equilibrium_residual_native(mix._native, request)
    except Exception as exc:
        return {
            "diagnostic_entrypoint": "_evaluate_chemical_equilibrium_residual_native",
            "exception_type": type(exc).__name__,
            "exception_message": str(exc),
        }
    return {
        "diagnostic_entrypoint": "_evaluate_chemical_equilibrium_residual_native",
        "residual_norm": float(payload.get("residual_norm", math.nan)),
        "diagnostics": payload.get("diagnostics", {}),
    }


def _summary(
    *,
    status: str,
    row: dict[str, str],
    solve_payload: dict[str, Any] | None,
    error: Exception | None,
    diagnostic: dict[str, Any] | None,
) -> dict[str, Any]:
    source_paths = [
        DATASET_DIR / "pure" / "any_solvent.csv",
        DATASET_DIR / "mixed" / "binary_interaction" / "k_ij.csv",
        DATASET_DIR / "mixed" / "binary_interaction" / "k_hb_ij.csv",
        DATASET_DIR / "mixed" / "binary_interaction" / "l_ij.csv",
        DATASET_DIR / "mixed" / "rel_perm" / "parameters.csv",
        DATASET_DIR / "user_options.json",
        JOU_VLE_CSV,
        PROBLEM_JSON,
        REACTION_CANDIDATES_CSV,
        REACTION_BASIS_CSV,
        SOURCE_VERIFICATION_CSV,
        SOURCE_MANIFEST_CSV,
    ]
    summary = {
        "lane": "mea_co2_pressure_speciation",
        "status": status,
        "source_paths": [_rel(path) for path in source_paths],
        "selected_source_row": row,
        "public_api": "epcsaft.solve_reactive_speciation",
        "required_solver_backend": "ipopt",
        "required_derivative_backend": "cppad_implicit",
        "accepted_native_ipopt_speciation": bool(solve_payload and solve_payload["accepted_native_ipopt_speciation"]),
        "retained_outputs": [
            _rel(PROCESSED_DIR / "selected_jou_row.json"),
            _rel(PROCESSED_DIR / "normalized_user_options.json"),
            _rel(PROCESSED_DIR / "reaction_constants.json"),
            _rel(SUMMARY_JSON),
        ],
        "solve": solve_payload,
        "diagnostic": diagnostic,
        "vapor_side_assumption": "ideal_fugacity",
    }
    if error is not None:
        message = str(error)
        blocker = message
        diagnostic_message = str(diagnostic.get("exception_message", "")) if diagnostic else ""
        if "Born model 2 phase-state composition sensitivity" in diagnostic_message:
            blocker = "SSM/DS Born composition sensitivity is not implemented for phase-state residuals."
        summary.update(
            {
                "blocker": blocker,
                "exception_type": type(error).__name__,
                "exception_message": message,
            }
        )
    return summary


def main() -> int:
    row = _selected_jou_row()
    user_options = _normalized_user_options()
    error: Exception | None = None
    solve_payload: dict[str, Any] | None = None
    try:
        solve_payload = _attempt_public_solve(row, user_options)
    except Exception as exc:
        error = exc
    diagnostic = None
    if error is not None:
        diagnostic = _diagnose_native_residual(row, user_options)
    status = "executable" if solve_payload and solve_payload["accepted_native_ipopt_speciation"] else "blocked"
    summary = _summary(status=status, row=row, solve_payload=solve_payload, error=error, diagnostic=diagnostic)
    _write_json(SUMMARY_JSON, summary)
    print(json.dumps(_json_like(summary), indent=2))
    return 0 if status == "executable" else 1


if __name__ == "__main__":
    raise SystemExit(main())
