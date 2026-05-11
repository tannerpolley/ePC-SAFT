from __future__ import annotations

import csv
from functools import lru_cache
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
REPO_ROOT = Path(__file__).resolve().parents[4]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

import _common as common
from scripts._env import require_epcsaft_install

require_epcsaft_install()

import epcsaft
from epcsaft.parameters import get_prop_dict

DATA_PATH = common.analysis_data_path(__file__, "1-butanol-NH4Cl-water-LLE.csv", kind="input")
VALIDATION_PATH = common.analysis_data_path(__file__, "model_validation.csv", kind="processed")
SPECIES = ["H2O", "Butanol", "NH4+", "Cl-"]
MW = np.asarray([18.0153e-3, 74.1216e-3, 18.038e-3, 35.453e-3], dtype=float)
MW_NH4CL = float(MW[2] + MW[3])
IDX = {name: i for i, name in enumerate(SPECIES)}
T_FIGURE = 298.15
P_FIGURE = 101325.0
MODEL_USER_OPTIONS = {
    "elec_model": {
        "rel_perm": {"rule": "linear-massfraction"},
        "include_born_model": False,
    }
}
TABLE_3_4_KIJ = {
    ("NH4+", "Cl-"): -0.566,
    ("NH4+", "H2O"): 0.064,
    ("NH4+", "Butanol"): 0.29,
    ("Cl-", "Butanol"): 0.22,
}
TABLE_4_LIJ = {
    ("NH4+", "Butanol"): 0.140,
    ("Cl-", "Butanol"): 0.245,
}


def _split_nh4cl_mass(total_salt_mass_fraction: float) -> tuple[float, float]:
    total = float(total_salt_mass_fraction)
    return total * float(MW[2] / MW_NH4CL), total * float(MW[3] / MW_NH4CL)


def _mass_to_mole_fraction(w: np.ndarray) -> np.ndarray:
    n = np.asarray(w, dtype=float) / MW
    return n / np.sum(n)


def _mole_to_mass_fraction(x: np.ndarray) -> np.ndarray:
    w = np.asarray(x, dtype=float) * MW
    return w / np.sum(w)


def _phase_mole_fraction_with_trace_salt(
    phase_mass_fraction: np.ndarray, trace_salt_mass_fraction: float
) -> np.ndarray:
    phase = np.asarray(phase_mass_fraction, dtype=float).copy()
    if trace_salt_mass_fraction > 0.0:
        m_nh4, m_cl = _split_nh4cl_mass(trace_salt_mass_fraction)
        phase[IDX["H2O"]] = max(phase[IDX["H2O"]] - trace_salt_mass_fraction, 1.0e-12)
        phase[IDX["NH4+"]] += m_nh4
        phase[IDX["Cl-"]] += m_cl
    return _mass_to_mole_fraction(phase)


def _set_symmetric_pair(
    params: dict, species: list[str], matrix_name: str, left: str, right: str, value: float
) -> None:
    if left not in species or right not in species:
        return
    matrix = params.get(matrix_name)
    if matrix is None:
        return
    i = species.index(left)
    j = species.index(right)
    matrix[i, j] = float(value)
    matrix[j, i] = float(value)


def _apply_held_figure6_overrides(params: dict, species: list[str]) -> dict:
    for (left, right), value in TABLE_3_4_KIJ.items():
        _set_symmetric_pair(params, species, "k_ij", left, right, value)
    for (left, right), value in TABLE_4_LIJ.items():
        _set_symmetric_pair(params, species, "l_ij", left, right, value)
    return params


def _held_figure6_params(feed_x: np.ndarray, user_options: dict = MODEL_USER_OPTIONS) -> dict:
    params = get_prop_dict("2014_Held", SPECIES, feed_x, T_FIGURE, user_options=user_options)
    return _apply_held_figure6_overrides(params, SPECIES)


def _solve_lle(
    aqueous_hint_mass_fraction: np.ndarray,
    organic_hint_mass_fraction: np.ndarray,
    *,
    beta_organic: float = 0.5,
) -> dict | None:
    if float(aqueous_hint_mass_fraction[IDX["NH4+"]] + aqueous_hint_mass_fraction[IDX["Cl-"]]) <= 1.0e-12:
        mw_neutral = MW[:2]
        aq0_neutral = _mass_to_mole_fraction(
            np.asarray([aqueous_hint_mass_fraction[0], aqueous_hint_mass_fraction[1], 0.0, 0.0])
        )
        org0_neutral = _mass_to_mole_fraction(
            np.asarray([organic_hint_mass_fraction[0], organic_hint_mass_fraction[1], 0.0, 0.0])
        )
        aq0 = aq0_neutral[:2] / np.sum(aq0_neutral[:2])
        org0 = org0_neutral[:2] / np.sum(org0_neutral[:2])
        beta = float(beta_organic)
        feed = (1.0 - beta) * aq0 + beta * org0
        neutral_species = ["H2O", "Butanol"]
        params = get_prop_dict("2014_Held", neutral_species, feed, T_FIGURE, user_options=MODEL_USER_OPTIONS)
        mixture = epcsaft.ePCSAFTMixture.from_params(
            _apply_held_figure6_overrides(params, neutral_species),
            species=neutral_species,
        )
        result = mixture.equilibrium(
            kind="lle_flash",
            T=T_FIGURE,
            P=P_FIGURE,
            z=feed,
            initial_phases={"liq1": aq0, "liq2": org0, "phase_fraction": beta},
            options=epcsaft.EquilibriumOptions(
                max_iterations=400,
                tolerance=1.0e-10,
                damping=0.5,
                include_phase_diagnostics=True,
                return_best_effort=True,
            ),
        )
        if not result.split_detected or len(result.phases) != 2:
            return {
                "aqueous_mass_fraction": None,
                "organic_mass_fraction": None,
                "diagnostics": dict(result.diagnostics),
            }
        phases = sorted(result.phases, key=lambda phase: float(phase.composition[0]), reverse=True)
        aq_w = phases[0].composition * mw_neutral
        org_w = phases[1].composition * mw_neutral
        aqueous = np.asarray([*(aq_w / np.sum(aq_w)), 0.0, 0.0], dtype=float)
        organic = np.asarray([*(org_w / np.sum(org_w)), 0.0, 0.0], dtype=float)
        return {
            "aqueous_mass_fraction": aqueous,
            "organic_mass_fraction": organic,
            "beta_aqueous": float(phases[0].phase_fraction),
            "beta_organic": float(phases[1].phase_fraction),
            "diagnostics": dict(result.diagnostics),
            "dielectric_rule": "linear-massfraction",
        }

    trace_salt = 1.0e-8

    aq0 = _phase_mole_fraction_with_trace_salt(aqueous_hint_mass_fraction, 0.0 if trace_salt < 1.0e-7 else trace_salt)
    org0 = _phase_mole_fraction_with_trace_salt(organic_hint_mass_fraction, trace_salt)
    beta = float(beta_organic)
    feed_x = (1.0 - beta) * aq0 + beta * org0
    feed_x = feed_x / np.sum(feed_x)

    mixture = epcsaft.ePCSAFTMixture.from_params(_held_figure6_params(feed_x), species=SPECIES)
    try:
        result = mixture.equilibrium(
            kind="electrolyte_lle",
            T=T_FIGURE,
            P=P_FIGURE,
            z=feed_x,
            initial_phases={"aq": aq0, "org": org0, "phase_fraction": beta},
            options=epcsaft.EquilibriumOptions(
                max_iterations=400,
                tolerance=1.0e-8,
                damping=0.5,
                include_phase_diagnostics=True,
                return_best_effort=True,
            ),
        )
    except epcsaft.SolutionError as exc:
        diagnostics = dict(getattr(exc, "diagnostics", {}) or {})
        result = None
    if result is None:
        return {
            "aqueous_mass_fraction": None,
            "organic_mass_fraction": None,
            "diagnostics": diagnostics,
        }
    if not result.split_detected or len(result.phases) != 2:
        return {
            "aqueous_mass_fraction": None,
            "organic_mass_fraction": None,
            "diagnostics": dict(result.diagnostics),
        }

    phases = {phase.label: phase for phase in result.phases}
    aqueous = phases.get("aq")
    organic = phases.get("org")
    if aqueous is None or organic is None:
        return None
    return {
        "aqueous_mass_fraction": _mole_to_mass_fraction(aqueous.composition),
        "organic_mass_fraction": _mole_to_mass_fraction(organic.composition),
        "beta_aqueous": float(aqueous.phase_fraction),
        "beta_organic": float(organic.phase_fraction),
        "diagnostics": dict(result.diagnostics),
        "dielectric_rule": "linear-massfraction",
    }


@lru_cache(maxsize=1)
def load_experimental_rows() -> tuple[dict, ...]:
    rows: list[dict] = []
    with DATA_PATH.open("r", newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            try:
                w_nh4cl_aq = float(row["w%_NH4Cl in aqueous phase"]) / 100.0
                w_buoh_org = float(row["w%_1-butanol in 1-butanol phase"]) / 100.0
                w_buoh_aq = float(row["w%_1-butanol in aqueous phase"]) / 100.0
            except (KeyError, TypeError, ValueError):
                continue

            m_nh4, m_cl = _split_nh4cl_mass(w_nh4cl_aq)
            aqueous = np.asarray(
                [
                    1.0 - w_nh4cl_aq - w_buoh_aq,
                    w_buoh_aq,
                    m_nh4,
                    m_cl,
                ],
                dtype=float,
            )
            # The provided Figure 6 data only reports NH4Cl in the aqueous phase and 1-butanol in both phases.
            # For plotting the experimental ternary points, the remaining organic fraction is treated as water.
            organic = np.asarray([1.0 - w_buoh_org, w_buoh_org, 0.0, 0.0], dtype=float)

            rows.append(
                {
                    "w_nh4cl_aq": w_nh4cl_aq,
                    "w_buoh_org": w_buoh_org,
                    "w_buoh_aq": w_buoh_aq,
                    "aqueous_mass_fraction": aqueous,
                    "organic_mass_fraction": organic,
                }
            )
    return tuple(rows)


def _objective(target: dict, solved: dict) -> float:
    aq = np.asarray(solved["aqueous_mass_fraction"], dtype=float)
    org = np.asarray(solved["organic_mass_fraction"], dtype=float)
    target_aq = np.asarray(target["aqueous_mass_fraction"], dtype=float)
    target_org = np.asarray(target["organic_mass_fraction"], dtype=float)

    aq_salt = float(aq[IDX["NH4+"]] + aq[IDX["Cl-"]])
    org_salt = float(org[IDX["NH4+"]] + org[IDX["Cl-"]])
    err = 0.0
    err += 3.0 * (aq_salt - float(target["w_nh4cl_aq"])) ** 2
    err += 1.5 * (float(aq[IDX["Butanol"]]) - float(target_aq[IDX["Butanol"]])) ** 2
    err += 0.8 * (float(org[IDX["Butanol"]]) - float(target_org[IDX["Butanol"]])) ** 2
    err += 0.2 * org_salt**2
    if float(org[IDX["Butanol"]]) < float(aq[IDX["Butanol"]]):
        err += 10.0
    return float(err)


@lru_cache(maxsize=1)
def solve_model_rows() -> tuple[dict, ...]:
    solved_rows: list[dict] = []
    for idx, row in enumerate(load_experimental_rows()):
        aq = np.asarray(row["aqueous_mass_fraction"], dtype=float)
        org = np.asarray(row["organic_mass_fraction"], dtype=float)

        best: tuple[float, float, dict] | None = None
        for beta_organic in (0.35, 0.50, 0.65):
            solved = _solve_lle(aq, org, beta_organic=beta_organic)
            if solved is None or solved["aqueous_mass_fraction"] is None or solved["organic_mass_fraction"] is None:
                continue
            score = _objective(row, solved)
            if best is None or score < best[0]:
                best = (score, float(beta_organic), solved)

        if best is None:
            continue

        score, beta_organic, solved = best
        aq_model = np.asarray(solved["aqueous_mass_fraction"], dtype=float)
        org_model = np.asarray(solved["organic_mass_fraction"], dtype=float)
        solved_rows.append(
            {
                "row_id": idx,
                "feed_beta_organic": beta_organic,
                "objective": score,
                "aqueous_mass_fraction": aq_model,
                "organic_mass_fraction": org_model,
                "w_nh4cl_aq": float(aq_model[IDX["NH4+"]] + aq_model[IDX["Cl-"]]),
                "w_buoh_org": float(org_model[IDX["Butanol"]]),
                "w_buoh_aq": float(aq_model[IDX["Butanol"]]),
                "beta_organic": float(solved["beta_organic"]),
                "beta_aqueous": float(solved["beta_aqueous"]),
                "solver_residual_norm": float(solved["diagnostics"].get("solver_residual_norm", np.nan)),
                "gibbs_delta": float(solved["diagnostics"].get("gibbs_delta", np.nan)),
                "dielectric_rule": str(solved.get("dielectric_rule", "linear-massfraction")),
            }
        )

    solved_rows.sort(key=lambda item: item["w_nh4cl_aq"])
    return tuple(solved_rows)


def write_model_validation_table() -> Path:
    rows: list[dict[str, float | int | str]] = []
    experimental = load_experimental_rows()
    for model in solve_model_rows():
        exp = experimental[int(model["row_id"])]
        rows.append(
            {
                "row_id": int(model["row_id"]),
                "dielectric_rule": str(model["dielectric_rule"]),
                "w_NH4Cl_aq_exp_wt_pct": 100.0 * float(exp["w_nh4cl_aq"]),
                "w_NH4Cl_aq_model_wt_pct": 100.0 * float(model["w_nh4cl_aq"]),
                "w_1butanol_org_exp_wt_pct": 100.0 * float(exp["w_buoh_org"]),
                "w_1butanol_org_model_wt_pct": 100.0 * float(model["w_buoh_org"]),
                "w_1butanol_aq_exp_wt_pct": 100.0 * float(exp["w_buoh_aq"]),
                "w_1butanol_aq_model_wt_pct": 100.0 * float(model["w_buoh_aq"]),
                "abs_error_1butanol_org_wt_pct": abs(100.0 * (float(model["w_buoh_org"]) - float(exp["w_buoh_org"]))),
                "abs_error_1butanol_aq_wt_pct": abs(100.0 * (float(model["w_buoh_aq"]) - float(exp["w_buoh_aq"]))),
                "solver_residual_norm": float(model["solver_residual_norm"]),
                "gibbs_delta": float(model["gibbs_delta"]),
            }
        )
    VALIDATION_PATH.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(rows[0]) if rows else []
    with VALIDATION_PATH.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    return VALIDATION_PATH
