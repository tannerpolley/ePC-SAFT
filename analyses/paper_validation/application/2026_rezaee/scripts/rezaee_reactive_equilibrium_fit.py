from __future__ import annotations

import json
import math
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy.optimize import least_squares

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import _paths  # noqa: F401,E402
import rezaee_reactive_equilibrium_replay as replay  # noqa: E402
from epcsaft import ePCSAFTMixture  # noqa: E402

ANALYSIS_DIR = Path(__file__).resolve().parents[1]
PROCESSED_DIR = ANALYSIS_DIR / "data" / "processed"
RESULTS_DIR = ANALYSIS_DIR / "results" / "reaction_equilibrium"

FIT_CSV = PROCESSED_DIR / "rezaee_2026_reactive_equilibrium_fit.csv"
FIT_SUMMARY_JSON = RESULTS_DIR / "rezaee_2026_reactive_equilibrium_fit_summary.json"
FIT_REPORT_MD = RESULTS_DIR / "rezaee_2026_reactive_equilibrium_fit.md"

PAIRS = ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3))
BASE_VECTOR = np.asarray(
    [
        11.1223,
        3.8327,
        392.09,
        11.1504,
        4.0254,
        427.44,
        0.0623,
        0.0104,
        0.0158,
        0.0115,
        -0.0139,
        -0.0127,
    ],
    dtype=float,
)
LOWER = np.asarray([1.0, 2.0, 50.0, 1.0, 2.0, 50.0, -0.8, -0.8, -0.8, -0.8, -0.8, -0.8])
UPPER = np.asarray([30.0, 6.0, 1000.0, 30.0, 6.0, 1000.0, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8])


def _jsonable(value: Any) -> Any:
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, dict):
        return {str(k): _jsonable(v) for k, v in value.items()}
    if isinstance(value, list):
        return [_jsonable(v) for v in value]
    return value


def _build_organic(theta: np.ndarray, base_params: dict[str, Any]) -> tuple[ePCSAFTMixture, np.ndarray]:
    params: dict[str, Any] = {}
    for key, value in base_params.items():
        if key in {"m", "s", "e"}:
            arr = np.asarray(value, dtype=float).copy()
            if key == "m":
                arr[2], arr[3] = theta[0], theta[3]
            elif key == "s":
                arr[2], arr[3] = theta[1], theta[4]
            else:
                arr[2], arr[3] = theta[2], theta[5]
            params[key] = arr
        elif key == "k_ij":
            kij = np.zeros((4, 4), dtype=float)
            for val, (i, j) in zip(theta[6:12], PAIRS):
                kij[i, j] = kij[j, i] = float(val)
            params[key] = kij
        else:
            params[key] = value

    mix = ePCSAFTMixture.from_params(params, species=replay.ORG_LABELS)
    pure_ln_phi: list[float] = []
    for i, label in enumerate(replay.ORG_LABELS):
        pure_params: dict[str, Any] = {}
        for key, value in params.items():
            if key == "assoc_scheme":
                pure_params[key] = [value[i]]
            elif key == "k_ij":
                pure_params[key] = np.zeros((1, 1), dtype=float)
            else:
                pure_params[key] = np.asarray([value[i]], dtype=float)
        pure_state = ePCSAFTMixture.from_params(pure_params, species=[label]).state(
            T=replay.TEMPERATURE_K,
            x=np.asarray([1.0]),
            P=replay.PRESSURE_PA,
        )
        pure_ln_phi.append(float(pure_state.fugacity_coefficient()[0]))
    return mix, np.asarray(pure_ln_phi, dtype=float)


def _rows(aqueous_mix: ePCSAFTMixture) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for row in pd.read_csv(replay.EQUILIBRIUM_CSV).itertuples(index=False):
        aqueous_x = replay._aqueous_x(row)
        organic_x = replay._organic_x(row)
        aqueous_gamma = aqueous_mix.state(
            T=replay.TEMPERATURE_K,
            x=aqueous_x,
            P=replay.PRESSURE_PA,
        ).activity_coefficient(species=replay.AQ_LABELS)
        rows.append(
            {
                "experiment_no": int(row.experiment_no),
                "aqueous_x": aqueous_x,
                "organic_x": organic_x,
                "aqueous_gamma": aqueous_gamma,
                "source": row.source,
            }
        )
    return rows


def _lnq_rows(
    rows: list[dict[str, Any]],
    organic_mix: ePCSAFTMixture,
    pure_ln_phi: np.ndarray,
    ln_k: tuple[float, float],
) -> list[dict[str, float]]:
    out: list[dict[str, float]] = []
    for row in rows:
        organic_gamma = replay._organic_activity_coefficients(
            organic_mix,
            pure_ln_phi,
            row["organic_x"],
        )
        ln_q_li, ln_q_na = replay._reaction_log_quotients(
            row["aqueous_x"],
            row["organic_x"],
            row["aqueous_gamma"],
            organic_gamma,
        )
        out.append(
            {
                "experiment_no": row["experiment_no"],
                "lnQ_Li": float(ln_q_li),
                "lnQ_Na": float(ln_q_na),
                "lnQ_minus_lnK_Li": float(ln_q_li - ln_k[0]),
                "lnQ_minus_lnK_Na": float(ln_q_na - ln_k[1]),
            }
        )
    return out


def _residuals(
    theta: np.ndarray,
    rows: list[dict[str, Any]],
    base_params: dict[str, Any],
    ln_k: tuple[float, float] | None,
) -> np.ndarray:
    try:
        organic_mix, pure_ln_phi = _build_organic(theta[:12], base_params)
        active_ln_k = (float(theta[12]), float(theta[13])) if ln_k is None else ln_k
        values = _lnq_rows(rows, organic_mix, pure_ln_phi, active_ln_k)
        residual: list[float] = []
        for row in values:
            residual.extend([row["lnQ_minus_lnK_Li"], row["lnQ_minus_lnK_Na"]])
        return np.asarray(residual, dtype=float)
    except Exception:  # noqa: BLE001
        return np.full(len(rows) * 2, 1.0e6, dtype=float)


def _metrics(residual: np.ndarray) -> dict[str, float]:
    values = np.asarray(residual, dtype=float)
    abs_values = np.abs(values)
    return {
        "mean_abs_ln_residual": float(np.mean(abs_values)),
        "median_abs_ln_residual": float(np.median(abs_values)),
        "max_abs_ln_residual": float(np.max(abs_values)),
        "rms_ln_residual": float(math.sqrt(float(np.mean(values * values)))),
    }


def _required_gamma_diagnostics(
    rows: list[dict[str, Any]],
    organic_mix: ePCSAFTMixture,
    pure_ln_phi: np.ndarray,
    constants: dict[str, float],
) -> dict[str, Any]:
    required_rows: list[dict[str, float]] = []
    for row in rows:
        aqueous_x = np.asarray(row["aqueous_x"], dtype=float)
        organic_x = np.asarray(row["organic_x"], dtype=float)
        aqueous_gamma = row["aqueous_gamma"]
        organic_gamma = replay._organic_activity_coefficients(
            organic_mix,
            pure_ln_phi,
            organic_x,
        )
        required_rli = constants["Li"] * (
            (aqueous_x[1] * aqueous_gamma["Li+"])
            * (aqueous_x[5] * aqueous_gamma["OH-"])
            * (organic_x[0] * organic_gamma["DES"])
        ) / ((organic_x[2]) * (aqueous_x[0] * aqueous_gamma["H2O"]))
        required_rna = constants["Na"] * (
            (aqueous_x[2] * aqueous_gamma["Na+"])
            * (aqueous_x[5] * aqueous_gamma["OH-"])
            * (organic_x[0] * organic_gamma["DES"])
        ) / ((organic_x[3]) * (aqueous_x[0] * aqueous_gamma["H2O"]))
        required_rows.append(
            {
                "experiment_no": float(row["experiment_no"]),
                "gamma_RLi_package": float(organic_gamma["RLi"]),
                "gamma_RLi_required": float(required_rli),
                "ln_gamma_RLi_package_over_required": float(
                    math.log(max(organic_gamma["RLi"], 1.0e-300))
                    - math.log(max(required_rli, 1.0e-300))
                ),
                "gamma_RNa_package": float(organic_gamma["RNa"]),
                "gamma_RNa_required": float(required_rna),
                "ln_gamma_RNa_package_over_required": float(
                    math.log(max(organic_gamma["RNa"], 1.0e-300))
                    - math.log(max(required_rna, 1.0e-300))
                ),
            }
        )

    frame = pd.DataFrame(required_rows)
    return {
        "median_gamma_RLi_package": float(frame["gamma_RLi_package"].median()),
        "median_gamma_RLi_required": float(frame["gamma_RLi_required"].median()),
        "median_ln_gamma_RLi_package_over_required": float(
            frame["ln_gamma_RLi_package_over_required"].median()
        ),
        "median_gamma_RNa_package": float(frame["gamma_RNa_package"].median()),
        "median_gamma_RNa_required": float(frame["gamma_RNa_required"].median()),
        "median_ln_gamma_RNa_package_over_required": float(
            frame["ln_gamma_RNa_package_over_required"].median()
        ),
        "rows": required_rows,
    }


def _fit(
    rows: list[dict[str, Any]],
    base_params: dict[str, Any],
    *,
    ln_k: tuple[float, float] | None,
    initial_ln_k: tuple[float, float],
) -> Any:
    if ln_k is None:
        x0 = np.concatenate([BASE_VECTOR, np.asarray(initial_ln_k, dtype=float)])
        lower = np.concatenate([LOWER, np.asarray([-40.0, -40.0])])
        upper = np.concatenate([UPPER, np.asarray([40.0, 40.0])])
    else:
        x0 = BASE_VECTOR.copy()
        lower = LOWER
        upper = UPPER
    return least_squares(
        lambda theta: _residuals(theta, rows, base_params, ln_k),
        x0,
        bounds=(lower, upper),
        max_nfev=80,
        xtol=1.0e-7,
        ftol=1.0e-7,
        gtol=1.0e-7,
    )


def _parameter_payload(theta: np.ndarray) -> dict[str, Any]:
    return {
        "RLi": {"m": float(theta[0]), "sigma_A": float(theta[1]), "epsilon_over_k_K": float(theta[2])},
        "RNa": {"m": float(theta[3]), "sigma_A": float(theta[4]), "epsilon_over_k_K": float(theta[5])},
        "organic_binary_interactions": {
            "DES_TOPO": float(theta[6]),
            "DES_RLi": float(theta[7]),
            "DES_RNa": float(theta[8]),
            "TOPO_RLi": float(theta[9]),
            "TOPO_RNa": float(theta[10]),
            "RLi_RNa": float(theta[11]),
        },
    }


def _write_report(summary: dict[str, Any]) -> None:
    lines = [
        "# Rezaee 2026 Reactive Equilibrium Fit Diagnostics",
        "",
        "## Source Basis",
        "",
        "- 2025 main paper and supporting information supply the experimental phase-composition rows.",
        "- 2026 main paper and supporting information supply the reaction equations, Gibbs-energy constants, organic parameters, and binary interactions.",
        "- This script checks Rezaee Eqs. 5/6 directly at the experimental SI phase compositions.",
        "",
        "## Diagnostic Results",
        "",
        f"- Published-parameter median abs ln residual: `{summary['published']['metrics']['median_abs_ln_residual']}`.",
        f"- Published-parameter mean abs ln residual: `{summary['published']['metrics']['mean_abs_ln_residual']}`.",
        f"- Refit with paper constants fixed median abs ln residual: `{summary['paper_constants_refit']['metrics']['median_abs_ln_residual']}`.",
        f"- Refit with paper constants fixed mean abs ln residual: `{summary['paper_constants_refit']['metrics']['mean_abs_ln_residual']}`.",
        f"- Diagnostic refit with constants free median abs ln residual: `{summary['free_constants_refit']['metrics']['median_abs_ln_residual']}`.",
        f"- Diagnostic refit with constants free fitted K Li/Na: `{summary['free_constants_refit']['equilibrium_constants']['Li']}`, `{summary['free_constants_refit']['equilibrium_constants']['Na']}`.",
        f"- Median required/package gamma gap for RLi: `{summary['published_required_gamma']['median_ln_gamma_RLi_package_over_required']}` ln units.",
        f"- Median required/package gamma gap for RNa: `{summary['published_required_gamma']['median_ln_gamma_RNa_package_over_required']}` ln units.",
        "",
        "## Interpretation",
        "",
        "The package can evaluate the aqueous ePC-SAFT and organic PC-SAFT activity terms needed by the Rezaee formulation. The published constants and published parameter table do not satisfy the SI phase-composition rows under the package's activity-reference convention. Allowing the organic RLi/RNa parameters and organic binary interactions to refit while keeping the paper constants fixed improves the residual but does not close it to the paper's reported extraction-percent AARD.",
        "",
        "The direct cause is now quantified: matching the SI RLi/RNa mole fractions while holding the published Table 2 constants would require organic complex activity coefficients many orders of magnitude smaller than the package computes from the published Table 8/9 parameters.",
        "",
        "The constants-free fit is diagnostic only: it shows that the remaining gap is a source/reference-state or implementation-convention issue, because the fitted constants move far away from Table 2. Do not present that constants-free fit as the published Rezaee thermodynamic model.",
    ]
    FIT_REPORT_MD.parent.mkdir(parents=True, exist_ok=True)
    FIT_REPORT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    constants = replay._reaction_constants()
    ln_k_paper = (math.log(constants["Li"]), math.log(constants["Na"]))
    aqueous_mix, _aqueous_charges = replay._aqueous_mixture()
    _organic_mix, base_params, pure_ln_phi = replay._organic_mixture()
    rows = _rows(aqueous_mix)

    published_rows = _lnq_rows(rows, _organic_mix, pure_ln_phi, ln_k_paper)
    published_residual = np.asarray(
        [v for row in published_rows for v in (row["lnQ_minus_lnK_Li"], row["lnQ_minus_lnK_Na"])],
        dtype=float,
    )
    required_gamma = _required_gamma_diagnostics(rows, _organic_mix, pure_ln_phi, constants)

    fixed_fit = _fit(rows, base_params, ln_k=ln_k_paper, initial_ln_k=ln_k_paper)
    fixed_theta = fixed_fit.x[:12]
    fixed_mix, fixed_pure = _build_organic(fixed_theta, base_params)
    fixed_rows = _lnq_rows(rows, fixed_mix, fixed_pure, ln_k_paper)
    fixed_residual = np.asarray(
        [v for row in fixed_rows for v in (row["lnQ_minus_lnK_Li"], row["lnQ_minus_lnK_Na"])],
        dtype=float,
    )

    initial_ln_k = (
        float(np.median([row["lnQ_Li"] for row in published_rows])),
        float(np.median([row["lnQ_Na"] for row in published_rows])),
    )
    free_fit = _fit(rows, base_params, ln_k=None, initial_ln_k=initial_ln_k)
    free_theta = free_fit.x[:12]
    free_ln_k = (float(free_fit.x[12]), float(free_fit.x[13]))
    free_mix, free_pure = _build_organic(free_theta, base_params)
    free_rows = _lnq_rows(rows, free_mix, free_pure, free_ln_k)
    free_residual = np.asarray(
        [v for row in free_rows for v in (row["lnQ_minus_lnK_Li"], row["lnQ_minus_lnK_Na"])],
        dtype=float,
    )

    table: list[dict[str, Any]] = []
    for source_row, pub, fixed, free in zip(rows, published_rows, fixed_rows, free_rows):
        table.append(
            {
                "experiment_no": source_row["experiment_no"],
                "published_lnQ_minus_lnK_Li": pub["lnQ_minus_lnK_Li"],
                "published_lnQ_minus_lnK_Na": pub["lnQ_minus_lnK_Na"],
                "paper_constants_refit_lnQ_minus_lnK_Li": fixed["lnQ_minus_lnK_Li"],
                "paper_constants_refit_lnQ_minus_lnK_Na": fixed["lnQ_minus_lnK_Na"],
                "free_constants_refit_lnQ_minus_lnK_Li": free["lnQ_minus_lnK_Li"],
                "free_constants_refit_lnQ_minus_lnK_Na": free["lnQ_minus_lnK_Na"],
                "source": source_row["source"],
            }
        )
    FIT_CSV.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(table).to_csv(FIT_CSV, index=False)

    summary = {
        "row_count": len(rows),
        "published": {
            "equilibrium_constants": constants,
            "parameters": _parameter_payload(BASE_VECTOR),
            "metrics": _metrics(published_residual),
        },
        "published_required_gamma": required_gamma,
        "paper_constants_refit": {
            "success": bool(fixed_fit.success),
            "nfev": int(fixed_fit.nfev),
            "message": str(fixed_fit.message),
            "equilibrium_constants": constants,
            "parameters": _parameter_payload(fixed_theta),
            "metrics": _metrics(fixed_residual),
        },
        "free_constants_refit": {
            "success": bool(free_fit.success),
            "nfev": int(free_fit.nfev),
            "message": str(free_fit.message),
            "equilibrium_constants": {
                "Li": float(math.exp(free_ln_k[0])),
                "Na": float(math.exp(free_ln_k[1])),
                "lnK_Li": free_ln_k[0],
                "lnK_Na": free_ln_k[1],
            },
            "parameters": _parameter_payload(free_theta),
            "metrics": _metrics(free_residual),
        },
        "status": "source_reference_state_gap",
    }
    FIT_SUMMARY_JSON.parent.mkdir(parents=True, exist_ok=True)
    FIT_SUMMARY_JSON.write_text(json.dumps(_jsonable(summary), indent=2, sort_keys=True) + "\n", encoding="utf-8")
    _write_report(summary)
    print(json.dumps(_jsonable(summary), indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
