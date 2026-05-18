"""Try the downstream Rezaee calibration through current package routes.

This is a diagnostic candidate, not a promoted validation gate.  It answers
whether the paper-constants refit calibration materially improves the native
residuals and whether the current public reactive electrolyte LLE route can
produce an accepted Ipopt solve from a source-backed Rezaee feed.
"""

from __future__ import annotations

import json
import math
import sys
from copy import deepcopy
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


SCRIPT_DIR = Path(__file__).resolve().parent
ANALYSIS_DIR = SCRIPT_DIR.parent
REPO_ROOT = ANALYSIS_DIR.parents[3]
sys.path.insert(0, str(SCRIPT_DIR))
sys.path.insert(0, str(REPO_ROOT))

from scripts.dev.native_runtime_env import apply_to_current_process  # noqa: E402

apply_to_current_process(cache_path=REPO_ROOT / "build" / "dev" / "CMakeCache.txt")

import epcsaft  # noqa: E402
from epcsaft.equilibrium import SolutionError  # noqa: E402

import rezaee_reactive_equilibrium_replay as replay  # noqa: E402


CALIBRATION_JSON = (
    ANALYSIS_DIR
    / "data"
    / "processed"
    / "rezaee_2026_paper_constants_refit_calibration.json"
)
ROWS_CSV = (
    ANALYSIS_DIR
    / "data"
    / "processed"
    / "rezaee_2026_calibrated_native_ipopt_attempt_rows.csv"
)
SEPARATE_PHASE_ROWS_CSV = (
    ANALYSIS_DIR
    / "data"
    / "processed"
    / "rezaee_2026_calibrated_separate_phase_residual_rows.csv"
)
NATIVE_LOGK_FIT_ROWS_CSV = (
    ANALYSIS_DIR
    / "data"
    / "processed"
    / "rezaee_2026_native_logk_fit_rows.csv"
)
NATIVE_LOGK_FIT_PUBLIC_ROUTE_ROWS_CSV = (
    ANALYSIS_DIR
    / "data"
    / "processed"
    / "rezaee_2026_native_logk_fit_public_route_rows.csv"
)
SUMMARY_JSON = (
    ANALYSIS_DIR
    / "results"
    / "reaction_equilibrium"
    / "rezaee_2026_calibrated_native_ipopt_attempt_summary.json"
)
REPORT_MD = (
    ANALYSIS_DIR
    / "results"
    / "reaction_equilibrium"
    / "rezaee_2026_calibrated_native_ipopt_attempt.md"
)

SPECIES_LABELS = replay.AQ_LABELS + replay.ORG_LABELS
CHARGES = np.array([0.0, 1.0, 1.0, -1.0, 1.0, -1.0, 1.0, 0.0, 0.0, 0.0, 0.0])
CL_INDEX = SPECIES_LABELS.index("Cl-")
H_INDEX = SPECIES_LABELS.index("H+")
FIT_EXPERIMENTS = tuple(range(1, 17)) + (25, 26)
HOLDOUT_EXPERIMENTS = tuple(range(17, 25))
TARGET_RESPONSES_CSV = (
    ANALYSIS_DIR
    / "data"
    / "input"
    / "rezaee_2025_doe_extraction_responses.csv"
)


@dataclass(frozen=True)
class PublicAttemptBasis:
    feed: np.ndarray
    neutralized_feed: np.ndarray
    source_charge: float
    neutralized_charge: float
    neutralization_species: str
    neutralization_delta: float
    balances: dict[str, dict[str, float]]
    totals: dict[str, float]


def _load_calibration() -> dict[str, Any]:
    with CALIBRATION_JSON.open(encoding="utf-8") as handle:
        return json.load(handle)


def _set_symmetric(matrix: Any, labels: list[str], left: str, right: str, value: float) -> None:
    i = labels.index(left)
    j = labels.index(right)
    matrix[i][j] = value
    matrix[j][i] = value


def _calibrated_organic_params(calibration: dict[str, Any]) -> dict[str, Any]:
    _, params, _ = replay._organic_mixture()
    params = deepcopy(params)

    for label in ("RLi", "RNa"):
        species_index = replay.ORG_LABELS.index(label)
        species_params = calibration["parameters"][label]
        params["m"][species_index] = species_params["m"]
        params["s"][species_index] = species_params["sigma_A"]
        params["e"][species_index] = species_params["epsilon_over_k_K"]

    pair_values = calibration["parameters"]["organic_binary_interactions"]
    pair_labels = {
        "DES_RLi": ("DES", "RLi"),
        "DES_RNa": ("DES", "RNa"),
        "DES_TOPO": ("DES", "TOPO"),
        "RLi_RNa": ("RLi", "RNa"),
        "TOPO_RLi": ("TOPO", "RLi"),
        "TOPO_RNa": ("TOPO", "RNa"),
    }
    for key, (left, right) in pair_labels.items():
        _set_symmetric(params["k_ij"], replay.ORG_LABELS, left, right, float(pair_values[key]))
    return params


def _calibrated_combined_mixture(calibration: dict[str, Any]) -> epcsaft.ePCSAFTMixture:
    aqueous_params = replay._aqueous_parameter_payload()
    organic_params = _calibrated_organic_params(calibration)
    return replay._combined_rezaee_mixture(aqueous_params, organic_params)


def _calibrated_organic_mixture(calibration: dict[str, Any]) -> tuple[epcsaft.ePCSAFTMixture, np.ndarray]:
    params = _calibrated_organic_params(calibration)
    mix = epcsaft.ePCSAFTMixture.from_params(params, species=replay.ORG_LABELS)
    pure_ln_phi = []
    for i, label in enumerate(replay.ORG_LABELS):
        pure_params: dict[str, Any] = {}
        for key, value in params.items():
            if key == "assoc_scheme":
                pure_params[key] = [value[i]]
            elif key == "k_ij":
                pure_params[key] = np.zeros((1, 1), dtype=float)
            else:
                pure_params[key] = np.asarray([value[i]], dtype=float)
        pure_mix = epcsaft.ePCSAFTMixture.from_params(pure_params, species=[label])
        pure_state = pure_mix.state(T=replay.TEMPERATURE_K, x=np.asarray([1.0]), P=replay.PRESSURE_PA)
        pure_ln_phi.append(float(pure_state.fugacity_coefficient()[0]))
    return mix, np.asarray(pure_ln_phi, dtype=float)


def _phase_vectors(row: pd.Series, min_composition: float = 1.0e-14) -> tuple[np.ndarray, np.ndarray]:
    aqueous_x = replay._aqueous_x(row)
    organic_x = replay._organic_x(row)
    phase1 = np.concatenate([aqueous_x, np.full(len(replay.ORG_LABELS), min_composition)])
    phase2 = np.concatenate([np.full(len(replay.AQ_LABELS), min_composition), organic_x])
    phase1 = phase1 / phase1.sum()
    phase2 = phase2 / phase2.sum()
    return phase1, phase2


def _neutralize_feed(feed: np.ndarray) -> tuple[np.ndarray, str, float]:
    charge = float(np.dot(feed, CHARGES))
    adjusted = feed.copy()
    if abs(charge) <= 1.0e-12:
        return adjusted / adjusted.sum(), "none", 0.0
    if charge > 0.0:
        adjusted[CL_INDEX] += charge
        species = "Cl-"
        delta = charge
    else:
        adjusted[H_INDEX] += -charge
        species = "H+"
        delta = -charge
    if adjusted.min() <= 0.0:
        raise ValueError("neutralized feed created a nonpositive composition")
    adjusted = adjusted / adjusted.sum()
    return adjusted, species, delta


def _balance_definitions() -> dict[str, dict[str, float]]:
    return {
        "aqueous_oxygen_pool": {"H2O": 1.0, "OH-": 1.0},
        "lithium_total": {"Li+": 1.0, "RLi": 1.0},
        "sodium_total": {"Na+": 1.0, "RNa": 1.0},
        "chloride_total": {"Cl-": 1.0},
        "proton_total": {"H+": 1.0},
        "ammonium_total": {"NH4+": 1.0},
        "organic_carrier_total": {"DES": 1.0, "RLi": 1.0, "RNa": 1.0},
        "topo_total": {"TOPO": 1.0},
    }


def _public_attempt_basis(row: pd.Series) -> PublicAttemptBasis:
    phase1, phase2 = _phase_vectors(row)
    feed = 0.5 * phase1 + 0.5 * phase2
    neutralized_feed, species, delta = _neutralize_feed(feed)
    balances = _balance_definitions()
    totals = {
        name: float(sum(coeff * neutralized_feed[SPECIES_LABELS.index(species)] for species, coeff in row_def.items()))
        for name, row_def in balances.items()
    }
    return PublicAttemptBasis(
        feed=feed,
        neutralized_feed=neutralized_feed,
        source_charge=float(np.dot(feed, CHARGES)),
        neutralized_charge=float(np.dot(neutralized_feed, CHARGES)),
        neutralization_species=species,
        neutralization_delta=float(delta),
        balances=balances,
        totals=totals,
    )


def _row_residuals(mixture: epcsaft.ePCSAFTMixture, constants: dict[str, float]) -> list[dict[str, Any]]:
    rows = pd.read_csv(replay.EQUILIBRIUM_CSV)
    records: list[dict[str, Any]] = []
    for _, row in rows.iterrows():
        aqueous_x = replay._aqueous_x(row)
        organic_x = replay._organic_x(row)
        cross_phase, diagnostics = replay._package_cross_phase_residuals(
            mixture,
            constants,
            aqueous_x,
            organic_x,
        )
        records.append(
            {
                "experiment_no": int(row["experiment_no"]),
                "row_temperature_K": replay.TEMPERATURE_K,
                "li_ln_residual": float(cross_phase[0]),
                "na_ln_residual": float(cross_phase[1]),
                "max_abs_residual": float(max(abs(float(cross_phase[0])), abs(float(cross_phase[1])))),
                "max_abs_charge_residual": _max_abs(diagnostics.get("charge_residuals", [])),
                "max_abs_phase_equilibrium_residual": _max_abs(
                    diagnostics.get("phase_equilibrium_residuals", [0.0])
                ),
                "phase_distance": float(diagnostics.get("phase_distance", float("nan"))),
            }
        )
    return records


def _max_abs(values: Any) -> float:
    values_list = [abs(float(value)) for value in values]
    return float(max(values_list)) if values_list else 0.0


def _separate_phase_residuals(
    organic_mix: epcsaft.ePCSAFTMixture,
    pure_ln_phi: np.ndarray,
    constants: dict[str, float],
) -> list[dict[str, Any]]:
    aqueous_mix, aqueous_charges = replay._aqueous_mixture()
    rows = pd.read_csv(replay.EQUILIBRIUM_CSV)
    records: list[dict[str, Any]] = []
    for row in rows.itertuples(index=False):
        aqueous_x = replay._aqueous_x(row)
        organic_x = replay._organic_x(row)
        aqueous_state = aqueous_mix.state(T=replay.TEMPERATURE_K, x=aqueous_x, P=replay.PRESSURE_PA)
        aqueous_gamma = aqueous_state.activity_coefficient(species=replay.AQ_LABELS)
        organic_gamma = replay._organic_activity_coefficients(organic_mix, pure_ln_phi, organic_x)
        ln_q_li, ln_q_na = replay._reaction_log_quotients(
            aqueous_x,
            organic_x,
            aqueous_gamma,
            organic_gamma,
        )
        li_residual = float(ln_q_li - math.log(constants["Li"]))
        na_residual = float(ln_q_na - math.log(constants["Na"]))
        records.append(
            {
                "experiment_no": int(row.experiment_no),
                "li_ln_residual": li_residual,
                "na_ln_residual": na_residual,
                "max_abs_residual": max(abs(li_residual), abs(na_residual)),
                "aqueous_charge_residual": float(np.dot(aqueous_x, aqueous_charges)),
                "gamma_DES": float(organic_gamma["DES"]),
                "gamma_RLi": float(organic_gamma["RLi"]),
                "gamma_RNa": float(organic_gamma["RNa"]),
            }
        )
    return records


def _summarize_route_result(result: Any) -> dict[str, Any]:
    diagnostics = dict(getattr(result, "diagnostics", {}) or {})
    phase_compositions = getattr(result, "phase_compositions", None)
    phase_fractions = getattr(result, "phase_fractions", None)
    return {
        "accepted": bool(getattr(result, "accepted", False)),
        "solver_backend": getattr(result, "solver_backend", None),
        "derivative_backend": diagnostics.get("derivative_backend"),
        "density_backend": diagnostics.get("density_backend"),
        "phase_fractions": None if phase_fractions is None else [float(value) for value in phase_fractions],
        "phase_compositions": None
        if phase_compositions is None
        else [[float(value) for value in phase] for phase in phase_compositions],
        "densities_mol_m3": diagnostics.get("densities_mol_m3"),
        "max_material_residual": diagnostics.get("max_material_residual"),
        "max_charge_residual": diagnostics.get("max_charge_residual"),
        "max_reaction_residual": diagnostics.get("max_reaction_residual"),
        "max_phase_equilibrium_residual": diagnostics.get("max_phase_equilibrium_residual"),
        "phase_distance": diagnostics.get("phase_distance"),
        "raw_diagnostics": diagnostics,
    }


def _public_route_attempt(
    mixture: epcsaft.ePCSAFTMixture,
    constants: dict[str, float],
    phase_models: dict[str, epcsaft.ePCSAFTMixture],
) -> dict[str, Any]:
    rows = pd.read_csv(replay.EQUILIBRIUM_CSV)
    row = rows.iloc[0]
    basis = _public_attempt_basis(row)
    options = epcsaft.EquilibriumOptions(
        max_iterations=120,
        tolerance=1.0e-8,
        min_composition=1.0e-14,
        timeout_seconds=45.0,
    )
    payload: dict[str, Any] = {
        "experiment_no": int(row["experiment_no"]),
        "source_charge": basis.source_charge,
        "neutralized_charge": basis.neutralized_charge,
        "neutralization_species": basis.neutralization_species,
        "neutralization_delta": basis.neutralization_delta,
        "phase_models": {
            "aq_species": list(phase_models["aq"].species),
            "org_species": list(phase_models["org"].species),
            "contract": 'phase_models={"aq": aqueous_mix, "org": organic_mix}',
        },
        "balances": basis.balances,
        "totals": basis.totals,
    }
    try:
        result = mixture.equilibrium(
            kind="reactive_electrolyte_lle",
            T=replay.TEMPERATURE_K,
            P=replay.PRESSURE_PA,
            z=basis.neutralized_feed,
            balances=basis.balances,
            totals=basis.totals,
            reactions=replay._rezaee_phase_tagged_reactions(constants),
            phase_models=phase_models,
            options=options,
        )
    except SolutionError as exc:
        payload.update(
            {
                "accepted": False,
                "error_type": type(exc).__name__,
                "error": str(exc),
                "diagnostics": getattr(exc, "diagnostics", None),
            }
        )
    except Exception as exc:  # pragma: no cover - diagnostic script path.
        payload.update(
            {
                "accepted": False,
                "error_type": type(exc).__name__,
                "error": str(exc),
            }
        )
    else:
        payload.update(_summarize_route_result(result))
    return payload


def _native_phase_model_log_quotients(
    aqueous_mix: epcsaft.ePCSAFTMixture,
    organic_mix: epcsaft.ePCSAFTMixture,
    row: pd.Series,
) -> tuple[float, float]:
    aqueous_x = replay._aqueous_x(row)
    organic_x = replay._organic_x(row)
    aqueous_state = aqueous_mix.state(T=replay.TEMPERATURE_K, x=aqueous_x, P=replay.PRESSURE_PA)
    organic_state = organic_mix.state(T=replay.TEMPERATURE_K, x=organic_x, P=replay.PRESSURE_PA)
    aqueous_ln_activity = np.log(np.clip(aqueous_x, 1.0e-300, None)) + np.asarray(
        aqueous_state.fugacity_coefficient(),
        dtype=float,
    )
    organic_ln_activity = np.log(np.clip(organic_x, 1.0e-300, None)) + np.asarray(
        organic_state.fugacity_coefficient(),
        dtype=float,
    )
    li_ln_q = (
        organic_ln_activity[replay.ORG_LABELS.index("RLi")]
        + aqueous_ln_activity[replay.AQ_LABELS.index("H2O")]
        - aqueous_ln_activity[replay.AQ_LABELS.index("Li+")]
        - aqueous_ln_activity[replay.AQ_LABELS.index("OH-")]
        - organic_ln_activity[replay.ORG_LABELS.index("DES")]
    )
    na_ln_q = (
        organic_ln_activity[replay.ORG_LABELS.index("RNa")]
        + aqueous_ln_activity[replay.AQ_LABELS.index("H2O")]
        - aqueous_ln_activity[replay.AQ_LABELS.index("Na+")]
        - aqueous_ln_activity[replay.AQ_LABELS.index("OH-")]
        - organic_ln_activity[replay.ORG_LABELS.index("DES")]
    )
    return float(li_ln_q), float(na_ln_q)


def _fit_native_log_constants(
    aqueous_mix: epcsaft.ePCSAFTMixture,
    organic_mix: epcsaft.ePCSAFTMixture,
) -> tuple[dict[str, float], dict[str, float], pd.DataFrame]:
    rows = pd.read_csv(replay.EQUILIBRIUM_CSV)
    records: list[dict[str, Any]] = []
    for _, row in rows.iterrows():
        experiment_no = int(row["experiment_no"])
        li_ln_q, na_ln_q = _native_phase_model_log_quotients(aqueous_mix, organic_mix, row)
        records.append(
            {
                "experiment_no": experiment_no,
                "split": "fit" if experiment_no in FIT_EXPERIMENTS else "holdout",
                "lnQ_Li_native": li_ln_q,
                "lnQ_Na_native": na_ln_q,
            }
        )
    fit_df = pd.DataFrame.from_records(records)
    training = fit_df.loc[fit_df["split"] == "fit"]
    log_constants = {
        "Li": float(training["lnQ_Li_native"].mean()),
        "Na": float(training["lnQ_Na_native"].mean()),
    }
    constants = {key: math.exp(value) for key, value in log_constants.items()}
    fit_df["lnK_Li_fitted"] = log_constants["Li"]
    fit_df["lnK_Na_fitted"] = log_constants["Na"]
    fit_df["Li_ln_residual_to_fit"] = fit_df["lnQ_Li_native"] - log_constants["Li"]
    fit_df["Na_ln_residual_to_fit"] = fit_df["lnQ_Na_native"] - log_constants["Na"]
    return constants, log_constants, fit_df


def _target_response_by_experiment() -> dict[int, Any]:
    targets = pd.read_csv(TARGET_RESPONSES_CSV)
    return {int(row.experiment_no): row for row in targets.itertuples(index=False)}


def _route_extraction_metrics(result: Any, experiment_no: int, target_by_experiment: dict[int, Any]) -> dict[str, float]:
    phases = np.asarray([phase.composition for phase in result.phases], dtype=float)
    phase_totals = np.asarray(
        [float(phase.diagnostics["amount_total"]) for phase in result.phases],
        dtype=float,
    )
    amounts = phases * phase_totals[:, None]
    li_total = float(amounts[:, SPECIES_LABELS.index("Li+")].sum() + amounts[:, SPECIES_LABELS.index("RLi")].sum())
    na_total = float(amounts[:, SPECIES_LABELS.index("Na+")].sum() + amounts[:, SPECIES_LABELS.index("RNa")].sum())
    li_extraction = 100.0 * float(amounts[1, SPECIES_LABELS.index("RLi")]) / li_total
    na_extraction = 100.0 * float(amounts[1, SPECIES_LABELS.index("RNa")]) / na_total
    selectivity = li_extraction / na_extraction
    target = target_by_experiment[experiment_no]
    li_exp = float(target.li_extraction_pct_exp)
    na_exp = float(target.na_extraction_pct_exp)
    selectivity_exp = float(target.li_na_selectivity_exp)
    return {
        "li_extraction_pct_calc": li_extraction,
        "na_extraction_pct_calc": na_extraction,
        "selectivity_calc": selectivity,
        "li_extraction_pct_exp": li_exp,
        "na_extraction_pct_exp": na_exp,
        "selectivity_exp": selectivity_exp,
        "li_extraction_abs_pp": abs(li_extraction - li_exp),
        "na_extraction_abs_pp": abs(na_extraction - na_exp),
        "selectivity_abs_pct": 100.0 * abs(selectivity - selectivity_exp) / max(abs(selectivity_exp), 1.0e-300),
    }


def _fitted_native_logk_public_route(
    mixture: epcsaft.ePCSAFTMixture,
    aqueous_mix: epcsaft.ePCSAFTMixture,
    organic_mix: epcsaft.ePCSAFTMixture,
) -> dict[str, Any]:
    constants, log_constants, fit_df = _fit_native_log_constants(aqueous_mix, organic_mix)
    NATIVE_LOGK_FIT_ROWS_CSV.parent.mkdir(parents=True, exist_ok=True)
    fit_df.to_csv(NATIVE_LOGK_FIT_ROWS_CSV, index=False)
    target_by_experiment = _target_response_by_experiment()
    options = epcsaft.EquilibriumOptions(
        max_iterations=180,
        tolerance=1.0e-8,
        min_composition=1.0e-14,
        timeout_seconds=60.0,
    )
    rows = pd.read_csv(replay.EQUILIBRIUM_CSV)
    records: list[dict[str, Any]] = []
    for _, row in rows.iterrows():
        experiment_no = int(row["experiment_no"])
        split = "fit" if experiment_no in FIT_EXPERIMENTS else "holdout"
        basis = _public_attempt_basis(row)
        record: dict[str, Any] = {
            "experiment_no": experiment_no,
            "split": split,
            "accepted": False,
        }
        try:
            result = mixture.equilibrium(
                kind="reactive_electrolyte_lle",
                T=replay.TEMPERATURE_K,
                P=replay.PRESSURE_PA,
                z=basis.neutralized_feed,
                balances=basis.balances,
                totals=basis.totals,
                reactions=replay._rezaee_phase_tagged_reactions(constants),
                phase_models={"aq": aqueous_mix, "org": organic_mix},
                options=options,
            )
        except SolutionError as exc:
            diagnostics = getattr(exc, "diagnostics", None) or {}
            record.update(
                {
                    "error_type": type(exc).__name__,
                    "solver_status": diagnostics.get("solver_status"),
                    "route_status": diagnostics.get("reactive_route_status"),
                    "rejection_reason": diagnostics.get("rejection_reason"),
                    "reaction_stationarity_norm": diagnostics.get("reaction_stationarity_norm"),
                    "charge_balance_norm": diagnostics.get("charge_balance_norm"),
                }
            )
        except Exception as exc:  # pragma: no cover - retained diagnostic script path.
            record.update(
                {
                    "error_type": type(exc).__name__,
                    "error": str(exc),
                }
            )
        else:
            diagnostics = dict(getattr(result, "diagnostics", {}) or {})
            record.update(
                {
                    "accepted": True,
                    "solver_backend": getattr(result, "backend", None),
                    "derivative_backend": diagnostics.get("derivative_backend"),
                    "density_backend": diagnostics.get("density_backend"),
                    "solver_status": diagnostics.get("solver_status"),
                    "route_status": diagnostics.get("reactive_route_status"),
                    "rejection_reason": diagnostics.get("rejection_reason"),
                    "reaction_stationarity_norm": diagnostics.get("reaction_stationarity_norm"),
                    "charge_balance_norm": diagnostics.get("charge_balance_norm"),
                    "phase_distance": diagnostics.get("phase_distance"),
                }
            )
            record.update(_route_extraction_metrics(result, experiment_no, target_by_experiment))
        records.append(record)

    route_df = pd.DataFrame.from_records(records)
    route_df.to_csv(NATIVE_LOGK_FIT_PUBLIC_ROUTE_ROWS_CSV, index=False)
    holdout = route_df.loc[route_df["split"] == "holdout"]
    all_holdout_rows_solve = bool(holdout["accepted"].all())
    metrics: dict[str, Any] = {}
    if all_holdout_rows_solve:
        metrics = {
            "li_extraction_aard_pp": float(holdout["li_extraction_abs_pp"].mean()),
            "na_extraction_aard_pp": float(holdout["na_extraction_abs_pp"].mean()),
            "selectivity_aard_pct": float(holdout["selectivity_abs_pct"].mean()),
        }
    gate_passed = (
        all_holdout_rows_solve
        and metrics["li_extraction_aard_pp"] <= 20.0
        and metrics["na_extraction_aard_pp"] <= 20.0
        and metrics["selectivity_aard_pct"] <= 25.0
    )
    if gate_passed:
        status = "accepted_public_native_ipopt"
    elif all_holdout_rows_solve:
        status = "failed_gate"
    else:
        status = "blocked_solver"
    return {
        "status": status,
        "fit_method": "mean logK on native x*phi phase-model standard state",
        "constant_source": "fit_rows_1_16_25_26",
        "fit_rows": list(FIT_EXPERIMENTS),
        "holdout_rows": list(HOLDOUT_EXPERIMENTS),
        "fit_row_count": int((route_df["split"] == "fit").sum()),
        "holdout_row_count": int(len(holdout)),
        "all_fit_rows_solve": bool(route_df.loc[route_df["split"] == "fit", "accepted"].all()),
        "all_holdout_rows_solve": all_holdout_rows_solve,
        "log_equilibrium_constants": log_constants,
        "equilibrium_constants": constants,
        "fit_ln_residual_summary": {
            "Li_std": float(fit_df.loc[fit_df["split"] == "fit", "Li_ln_residual_to_fit"].std(ddof=0)),
            "Na_std": float(fit_df.loc[fit_df["split"] == "fit", "Na_ln_residual_to_fit"].std(ddof=0)),
            "holdout_Li_max_abs": float(fit_df.loc[fit_df["split"] == "holdout", "Li_ln_residual_to_fit"].abs().max()),
            "holdout_Na_max_abs": float(fit_df.loc[fit_df["split"] == "holdout", "Na_ln_residual_to_fit"].abs().max()),
        },
        "holdout_metrics": metrics,
        "holdout_requirement": {
            "all_holdout_rows_solve": True,
            "li_extraction_aard_pp_max": 20.0,
            "na_extraction_aard_pp_max": 20.0,
            "selectivity_aard_pct_max": 25.0,
        },
        "public_api": 'mix.equilibrium(kind="reactive_electrolyte_lle")',
        "phase_models": {
            "aq_species": list(aqueous_mix.species),
            "org_species": list(organic_mix.species),
            "contract": 'phase_models={"aq": aqueous_mix, "org": organic_mix}',
        },
        "outputs": {
            "logk_fit_rows_csv": str(NATIVE_LOGK_FIT_ROWS_CSV.relative_to(ANALYSIS_DIR)),
            "public_route_rows_csv": str(NATIVE_LOGK_FIT_PUBLIC_ROUTE_ROWS_CSV.relative_to(ANALYSIS_DIR)),
        },
    }


def _write_report(summary: dict[str, Any]) -> None:
    attempt = summary["public_route_attempt"]
    fitted = summary["fitted_native_logk_public_route"]
    lines = [
        "# Rezaee 2026 calibrated native Ipopt attempt",
        "",
        f"This diagnostic uses `{summary['calibration']['status']}`.",
        "It does not promote the Rezaee lane to direct closure.",
        "",
        "## Source-phase residuals",
        "",
        "Separate-phase replay convention:",
        f"- rows: {summary['separate_phase_residuals']['row_count']}",
        f"- median abs ln residual: {summary['separate_phase_residuals']['median_abs_ln_residual']:.6g}",
        f"- mean abs ln residual: {summary['separate_phase_residuals']['mean_abs_ln_residual']:.6g}",
        f"- max abs ln residual: {summary['separate_phase_residuals']['max_abs_ln_residual']:.6g}",
        "",
        "Native combined-mixture phase-tagged residual:",
        f"- rows: {summary['source_phase_residuals']['row_count']}",
        f"- median abs ln residual: {summary['source_phase_residuals']['median_abs_ln_residual']:.6g}",
        f"- mean abs ln residual: {summary['source_phase_residuals']['mean_abs_ln_residual']:.6g}",
        f"- max abs ln residual: {summary['source_phase_residuals']['max_abs_ln_residual']:.6g}",
        "",
        "## Public route attempt",
        "",
        f"- experiment_no: {attempt['experiment_no']}",
        f"- accepted: {attempt.get('accepted')}",
        f"- error_type: {attempt.get('error_type')}",
        f"- solver_backend: {attempt.get('solver_backend')}",
        f"- derivative_backend: {attempt.get('derivative_backend')}",
        f"- density_backend: {attempt.get('density_backend')}",
        f"- source_charge: {attempt['source_charge']:.6e}",
        f"- neutralized_charge: {attempt['neutralized_charge']:.6e}",
        f"- neutralization: {attempt['neutralization_delta']:.6e} added to {attempt['neutralization_species']}",
        "",
        "## Fitted public route gate",
        "",
        f"- status: {fitted['status']}",
        f"- fit_method: {fitted['fit_method']}",
        f"- all_fit_rows_solve: {fitted['all_fit_rows_solve']}",
        f"- all_holdout_rows_solve: {fitted['all_holdout_rows_solve']}",
        f"- holdout Li extraction AARD / percentage points: {fitted['holdout_metrics'].get('li_extraction_aard_pp')}",
        f"- holdout Na extraction AARD / percentage points: {fitted['holdout_metrics'].get('na_extraction_aard_pp')}",
        f"- holdout selectivity AARD / percent: {fitted['holdout_metrics'].get('selectivity_aard_pct')}",
    ]
    if attempt.get("error"):
        lines.extend(["", "### Error", "", attempt["error"]])
    SUMMARY_JSON.parent.mkdir(parents=True, exist_ok=True)
    REPORT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    calibration = _load_calibration()
    constants = {key: float(value) for key, value in calibration["equilibrium_constants"].items()}
    mixture = _calibrated_combined_mixture(calibration)
    aqueous_mix, _aqueous_charges = replay._aqueous_mixture()
    calibrated_organic_mix, calibrated_pure_ln_phi = _calibrated_organic_mixture(calibration)

    residual_records = _row_residuals(mixture, constants)
    rows_df = pd.DataFrame.from_records(residual_records)
    ROWS_CSV.parent.mkdir(parents=True, exist_ok=True)
    rows_df.to_csv(ROWS_CSV, index=False)

    abs_values = rows_df[["li_ln_residual", "na_ln_residual"]].abs().to_numpy().ravel()
    separate_records = _separate_phase_residuals(calibrated_organic_mix, calibrated_pure_ln_phi, constants)
    separate_df = pd.DataFrame.from_records(separate_records)
    separate_df.to_csv(SEPARATE_PHASE_ROWS_CSV, index=False)
    separate_abs_values = separate_df[["li_ln_residual", "na_ln_residual"]].abs().to_numpy().ravel()
    public_attempt = _public_route_attempt(
        mixture,
        constants,
        {"aq": aqueous_mix, "org": calibrated_organic_mix},
    )
    fitted_native_logk_public_route = _fitted_native_logk_public_route(
        mixture,
        aqueous_mix,
        calibrated_organic_mix,
    )
    summary = {
        "status": "diagnostic_attempt",
        "calibration": {
            "status": calibration["status"],
            "source_project": calibration["source_project"],
            "source_file": calibration["source_file"],
            "source_metrics": calibration["metrics"],
        },
        "source_phase_residuals": {
            "row_count": int(len(rows_df)),
            "median_abs_ln_residual": float(np.median(abs_values)),
            "mean_abs_ln_residual": float(np.mean(abs_values)),
            "max_abs_ln_residual": float(np.max(abs_values)),
            "rows_csv": str(ROWS_CSV.relative_to(ANALYSIS_DIR)),
        },
        "separate_phase_residuals": {
            "row_count": int(len(separate_df)),
            "median_abs_ln_residual": float(np.median(separate_abs_values)),
            "mean_abs_ln_residual": float(np.mean(separate_abs_values)),
            "max_abs_ln_residual": float(np.max(separate_abs_values)),
            "rows_csv": str(SEPARATE_PHASE_ROWS_CSV.relative_to(ANALYSIS_DIR)),
        },
        "public_route_attempt": public_attempt,
        "fitted_native_logk_public_route": fitted_native_logk_public_route,
        "outputs": {
            "summary_json": str(SUMMARY_JSON.relative_to(ANALYSIS_DIR)),
            "report_md": str(REPORT_MD.relative_to(ANALYSIS_DIR)),
            "rows_csv": str(ROWS_CSV.relative_to(ANALYSIS_DIR)),
            "separate_phase_rows_csv": str(SEPARATE_PHASE_ROWS_CSV.relative_to(ANALYSIS_DIR)),
            "native_logk_fit_rows_csv": str(NATIVE_LOGK_FIT_ROWS_CSV.relative_to(ANALYSIS_DIR)),
            "native_logk_fit_public_route_rows_csv": str(
                NATIVE_LOGK_FIT_PUBLIC_ROUTE_ROWS_CSV.relative_to(ANALYSIS_DIR)
            ),
        },
    }
    SUMMARY_JSON.parent.mkdir(parents=True, exist_ok=True)
    SUMMARY_JSON.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    _write_report(summary)

    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
