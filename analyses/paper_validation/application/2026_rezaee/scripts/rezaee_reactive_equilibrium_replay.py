from __future__ import annotations

import json
import math
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from _paths import ANALYSIS_DIR, REPO_ROOT

if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from _epcsaft_properties import get_prop_dict  # noqa: E402
from epcsaft import ReactionDefinition, _core, ePCSAFTMixture  # noqa: E402
from epcsaft.equilibrium import _reaction_phase_stoichiometry_matrix  # noqa: E402

INPUT_DIR = ANALYSIS_DIR / "data" / "input"
PROCESSED_DIR = ANALYSIS_DIR / "data" / "processed"
RESULTS_DIR = ANALYSIS_DIR / "results" / "reaction_equilibrium"

EQUILIBRIUM_CSV = INPUT_DIR / "rezaee_2025_extraction_equilibrium_mole_fractions.csv"
REACTION_CONSTANTS_CSV = INPUT_DIR / "rezaee_2026_reaction_constants.csv"
ORGANIC_PARAMS_CSV = INPUT_DIR / "rezaee_2026_organic_pcsaft_parameters.csv"
ORGANIC_KIJ_CSV = INPUT_DIR / "rezaee_2026_organic_binary_interactions.csv"

REPLAY_CSV = PROCESSED_DIR / "rezaee_2026_reactive_equilibrium_replay.csv"
SUMMARY_JSON = RESULTS_DIR / "rezaee_2026_reactive_equilibrium_replay_summary.json"
REPORT_MD = RESULTS_DIR / "rezaee_2026_reactive_equilibrium_replay.md"

TEMPERATURE_K = 298.15
PRESSURE_PA = 101325.0

AQ_PARAM_KEYS = ["H2O-Salt-2001", "Li+", "Na+", "Cl-", "H3O+", "OH-", "NH4+"]
AQ_LABELS = ["H2O", "Li+", "Na+", "Cl-", "H+", "OH-", "NH4+"]
ORG_LABELS = ["DES", "TOPO", "RLi", "RNa"]

UNSUPPORTED_FLAT_ELECTROLYTE_KEYS = {
    "dipm",
    "dip_num",
    "elec_model",
    "elec_model_preset",
    "bjeruum_treatment",
    "born_model",
    "born_diff_mode",
    "born_eps_mode",
    "DH_model",
    "dielc_rule",
    "dielc_diff_mode",
    "debug",
}


def _safe_log(value: float) -> float:
    return math.log(max(float(value), 1.0e-300))


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


def _reaction_constants() -> dict[str, float]:
    constants = pd.read_csv(REACTION_CONSTANTS_CSV)
    out: dict[str, float] = {}
    for row in constants.itertuples(index=False):
        reaction = str(row.reaction)
        if reaction.startswith("Li+"):
            out["Li"] = float(row.equilibrium_constant_298K)
        elif reaction.startswith("Na+"):
            out["Na"] = float(row.equilibrium_constant_298K)
    if set(out) != {"Li", "Na"}:
        raise ValueError(f"Expected Li and Na reaction constants in {REACTION_CONSTANTS_CSV}")
    return out


def _aqueous_mixture() -> tuple[ePCSAFTMixture, np.ndarray]:
    params = _aqueous_parameter_payload()
    return ePCSAFTMixture.from_params(params, species=AQ_LABELS), np.asarray(params["z"], dtype=float)


def _aqueous_parameter_payload() -> dict[str, Any]:
    # Rezaee's SI rows contain trace NH4+. The package reference catalog has
    # an NH4+ ion parameter in the 2022 Ascani family; this local table does not,
    # so we carry the same simple cation form here to keep charge accounting explicit.
    nh4_params = {
        "MW": 18.038e-3,
        "m": 1.0,
        "s": 3.5740,
        "e": 230.0,
        "e_assoc": 0.0,
        "vol_a": 0.0,
        "assoc_scheme": None,
        "dipm": 0.0,
        "dip_num": 1,
        "z": 1.0,
        "dielc": 8.0,
        "d_born": 3.0,
        "f_solv": 1.0,
    }
    params = get_prop_dict(
        AQ_PARAM_KEYS,
        np.full(len(AQ_PARAM_KEYS), 1.0 / len(AQ_PARAM_KEYS)),
        TEMPERATURE_K,
        user_params={"NH4+": nh4_params},
    )
    for key in UNSUPPORTED_FLAT_ELECTROLYTE_KEYS:
        params.pop(key, None)
    return params


def _organic_mixture() -> tuple[ePCSAFTMixture, dict[str, Any], np.ndarray]:
    params_df = pd.read_csv(ORGANIC_PARAMS_CSV)
    species = list(params_df["component"])
    if species != ORG_LABELS:
        raise ValueError(f"Expected organic species {ORG_LABELS}, got {species}")

    # Rezaee reports association site counts rather than the package's named
    # site schemes. Keep the mapping explicit so one-site TOPO is not treated
    # as a two-site 2B component.
    assoc_scheme_by_count = {0: None, 1: "1", 2: "2B"}
    assoc_scheme = [
        assoc_scheme_by_count.get(int(sites), "2B")
        for sites in params_df["association_sites"]
    ]
    params: dict[str, Any] = {
        "m": params_df["m"].to_numpy(dtype=float),
        "s": params_df["sigma_A"].to_numpy(dtype=float),
        "e": params_df["epsilon_over_k_K"].to_numpy(dtype=float),
        "e_assoc": params_df["epsilon_assoc_over_k_K"].to_numpy(dtype=float),
        "vol_a": params_df["kappa_assoc"].to_numpy(dtype=float),
        "assoc_scheme": assoc_scheme,
    }

    kij = np.zeros((len(species), len(species)), dtype=float)
    species_index = {label: i for i, label in enumerate(species)}
    for row in pd.read_csv(ORGANIC_KIJ_CSV).itertuples(index=False):
        i = species_index[str(row.component_i)]
        j = species_index[str(row.component_j)]
        kij[i, j] = kij[j, i] = float(row.kij)
    params["k_ij"] = kij

    mix = ePCSAFTMixture.from_params(params, species=species)
    pure_ln_phi = []
    for i, label in enumerate(species):
        pure_params: dict[str, Any] = {}
        for key, value in params.items():
            if key == "assoc_scheme":
                pure_params[key] = [value[i]]
            elif key == "k_ij":
                pure_params[key] = np.zeros((1, 1), dtype=float)
            else:
                pure_params[key] = np.asarray([value[i]], dtype=float)
        pure_mix = ePCSAFTMixture.from_params(pure_params, species=[label])
        pure_state = pure_mix.state(T=TEMPERATURE_K, x=np.asarray([1.0]), P=PRESSURE_PA)
        pure_ln_phi.append(float(pure_state.fugacity_coefficient()[0]))
    return mix, params, np.asarray(pure_ln_phi, dtype=float)


def _combined_rezaee_mixture(aqueous_params: dict[str, Any], organic_params: dict[str, Any]) -> ePCSAFTMixture:
    n_aq = len(AQ_LABELS)
    n_org = len(ORG_LABELS)
    params: dict[str, Any] = {}
    for key in ("m", "s", "e", "e_assoc", "vol_a"):
        params[key] = np.concatenate(
            [np.asarray(aqueous_params[key], dtype=float), np.asarray(organic_params[key], dtype=float)]
        )
    params["assoc_scheme"] = list(aqueous_params.get("assoc_scheme", [None] * n_aq)) + list(
        organic_params["assoc_scheme"]
    )
    params["z"] = np.concatenate([np.asarray(aqueous_params["z"], dtype=float), np.zeros(n_org)])
    params["dielc"] = np.concatenate(
        [np.asarray(aqueous_params.get("dielc", np.full(n_aq, 78.0)), dtype=float), np.full(n_org, 8.0)]
    )
    params["d_born"] = np.concatenate(
        [np.asarray(aqueous_params.get("d_born", np.full(n_aq, 3.0)), dtype=float), np.zeros(n_org)]
    )
    params["f_solv"] = np.concatenate(
        [np.asarray(aqueous_params.get("f_solv", np.ones(n_aq)), dtype=float), np.ones(n_org)]
    )
    if "MW" in aqueous_params:
        params["MW"] = np.concatenate([np.asarray(aqueous_params["MW"], dtype=float), np.full(n_org, 0.3)])
    kij = np.zeros((n_aq + n_org, n_aq + n_org), dtype=float)
    kij[n_aq:, n_aq:] = np.asarray(organic_params["k_ij"], dtype=float)
    params["k_ij"] = kij
    return ePCSAFTMixture.from_params(params, species=AQ_LABELS + ORG_LABELS)


def _rezaee_phase_tagged_reactions(constants: dict[str, float]) -> list[ReactionDefinition]:
    return [
        ReactionDefinition.from_literature_constant(
            {"Li+": -1.0, "OH-": -1.0, "DES": -1.0, "RLi": 1.0, "H2O": 1.0},
            log_equilibrium_constant=math.log(constants["Li"]),
            name="Rezaee_Li_cross_phase",
            standard_state="mole_fraction_activity",
            phase_stoichiometry={
                "aq": {"H2O": 1.0, "Li+": -1.0, "OH-": -1.0},
                "org": {"RLi": 1.0, "DES": -1.0},
            },
            source="Rezaee2026_SI_Table2",
        ),
        ReactionDefinition.from_literature_constant(
            {"Na+": -1.0, "OH-": -1.0, "DES": -1.0, "RNa": 1.0, "H2O": 1.0},
            log_equilibrium_constant=math.log(constants["Na"]),
            name="Rezaee_Na_cross_phase",
            standard_state="mole_fraction_activity",
            phase_stoichiometry={
                "aq": {"H2O": 1.0, "Na+": -1.0, "OH-": -1.0},
                "org": {"RNa": 1.0, "DES": -1.0},
            },
            source="Rezaee2026_SI_Table2",
        ),
    ]


def _package_cross_phase_residuals(
    mixture: ePCSAFTMixture,
    constants: dict[str, float],
    aqueous_x: np.ndarray,
    organic_x: np.ndarray,
) -> tuple[np.ndarray, dict[str, Any]]:
    species = AQ_LABELS + ORG_LABELS
    n_aq = len(AQ_LABELS)
    n_org = len(ORG_LABELS)
    phase1 = np.concatenate([aqueous_x, np.full(n_org, 1.0e-14)])
    phase1 = phase1 / float(np.sum(phase1))
    phase2 = np.concatenate([np.full(n_aq, 1.0e-14), organic_x])
    phase2 = phase2 / float(np.sum(phase2))
    feed = 0.5 * phase1 + 0.5 * phase2
    reactions = _rezaee_phase_tagged_reactions(constants)
    reaction_matrix = np.asarray(
        [[float(reaction.stoichiometry.get(label, 0.0)) for label in species] for reaction in reactions],
        dtype=float,
    )
    phase_matrix, phase_scope = _reaction_phase_stoichiometry_matrix(species, reactions, "electrolyte_lle")
    assert phase_matrix is not None
    request = {
        "T": TEMPERATURE_K,
        "P": PRESSURE_PA,
        "z": feed.tolist(),
        "initial_phases": {"aq": phase1.tolist(), "org": phase2.tolist(), "phase_fraction": 0.5},
        "balance_matrix": np.eye(len(species), dtype=float).reshape(-1).tolist(),
        "balance_rows": len(species),
        "total_vector": feed.tolist(),
        "reaction_stoichiometry": reaction_matrix.reshape(-1).tolist(),
        "reaction_rows": len(reactions),
        "log_equilibrium_constants": [float(reaction.log_equilibrium_constant) for reaction in reactions],
        "reaction_standard_states": [reaction.convention.native_standard_state_code for reaction in reactions],
        "reaction_phase_stoichiometry": phase_matrix.reshape(-1).tolist(),
        "options": {"min_composition": 1.0e-14, "tolerance": 1.0e-8},
    }
    payload = _core._evaluate_reactive_phase_equilibrium_residual_native(mixture._native, request)
    diagnostics = dict(payload["diagnostics"])
    diagnostics["python_phase_scope"] = phase_scope
    return np.asarray(payload["reaction_residuals_cross_phase"], dtype=float), diagnostics


def _aqueous_x(row: Any) -> np.ndarray:
    values = np.asarray(
        [
            row.aqueous_x_H2O,
            row.aqueous_x_Li_plus,
            row.aqueous_x_Na_plus,
            row.aqueous_x_Cl_minus,
            row.aqueous_x_H_plus,
            row.aqueous_x_OH_minus,
            row.aqueous_x_NH4_plus,
        ],
        dtype=float,
    )
    return values / float(np.sum(values))


def _organic_x(row: Any) -> np.ndarray:
    values = np.asarray(
        [row.organic_x_DES, row.organic_x_TOPO, row.organic_x_RLi, row.organic_x_RNa],
        dtype=float,
    )
    return values / float(np.sum(values))


def _organic_activity_coefficients(
    organic_mix: ePCSAFTMixture,
    pure_ln_phi: np.ndarray,
    organic_x: np.ndarray,
) -> dict[str, float]:
    state = organic_mix.state(T=TEMPERATURE_K, x=organic_x, P=PRESSURE_PA)
    ln_gamma = np.asarray(state.fugacity_coefficient(), dtype=float) - pure_ln_phi
    return {label: float(math.exp(value)) for label, value in zip(ORG_LABELS, ln_gamma.tolist())}


def _reaction_log_quotients(
    aqueous_x: np.ndarray,
    organic_x: np.ndarray,
    aqueous_gamma: dict[str, float],
    organic_gamma: dict[str, float],
) -> tuple[float, float]:
    ln_q_li = (
        _safe_log(organic_x[2])
        + _safe_log(organic_gamma["RLi"])
        + _safe_log(aqueous_x[0])
        + _safe_log(aqueous_gamma["H2O"])
        - _safe_log(aqueous_x[1])
        - _safe_log(aqueous_gamma["Li+"])
        - _safe_log(aqueous_x[5])
        - _safe_log(aqueous_gamma["OH-"])
        - _safe_log(organic_x[0])
        - _safe_log(organic_gamma["DES"])
    )
    ln_q_na = (
        _safe_log(organic_x[3])
        + _safe_log(organic_gamma["RNa"])
        + _safe_log(aqueous_x[0])
        + _safe_log(aqueous_gamma["H2O"])
        - _safe_log(aqueous_x[2])
        - _safe_log(aqueous_gamma["Na+"])
        - _safe_log(aqueous_x[5])
        - _safe_log(aqueous_gamma["OH-"])
        - _safe_log(organic_x[0])
        - _safe_log(organic_gamma["DES"])
    )
    return ln_q_li, ln_q_na


def _predict_complexes_from_constants(
    aqueous_x: np.ndarray,
    organic_x_exp: np.ndarray,
    aqueous_gamma: dict[str, float],
    organic_mix: ePCSAFTMixture,
    pure_ln_phi: np.ndarray,
    constants: dict[str, float],
    *,
    max_iter: int = 80,
    tolerance: float = 1.0e-10,
) -> tuple[np.ndarray, int, float]:
    base_ratio = organic_x_exp[:2] / float(np.sum(organic_x_exp[:2]))
    r_li = max(float(organic_x_exp[2]), 1.0e-12)
    r_na = max(float(organic_x_exp[3]), 1.0e-12)
    max_delta = math.inf

    for iteration in range(1, max_iter + 1):
        remainder = max(1.0e-12, 1.0 - r_li - r_na)
        organic_x = np.asarray([base_ratio[0] * remainder, base_ratio[1] * remainder, r_li, r_na], dtype=float)
        organic_gamma = _organic_activity_coefficients(organic_mix, pure_ln_phi, organic_x)
        next_li = constants["Li"] * (
            (aqueous_x[1] * aqueous_gamma["Li+"])
            * (aqueous_x[5] * aqueous_gamma["OH-"])
            * (organic_x[0] * organic_gamma["DES"])
        ) / ((organic_gamma["RLi"]) * (aqueous_x[0] * aqueous_gamma["H2O"]))
        next_na = constants["Na"] * (
            (aqueous_x[2] * aqueous_gamma["Na+"])
            * (aqueous_x[5] * aqueous_gamma["OH-"])
            * (organic_x[0] * organic_gamma["DES"])
        ) / ((organic_gamma["RNa"]) * (aqueous_x[0] * aqueous_gamma["H2O"]))
        next_li = min(max(float(next_li), 1.0e-15), 0.30)
        next_na = min(max(float(next_na), 1.0e-15), 0.30)
        max_delta = max(abs(next_li - r_li), abs(next_na - r_na))
        r_li = 0.5 * r_li + 0.5 * next_li
        r_na = 0.5 * r_na + 0.5 * next_na
        if max_delta < tolerance:
            break

    remainder = max(1.0e-12, 1.0 - r_li - r_na)
    predicted = np.asarray([base_ratio[0] * remainder, base_ratio[1] * remainder, r_li, r_na], dtype=float)
    return predicted, iteration, float(max_delta)


def _summaries(rows: pd.DataFrame, constants: dict[str, float]) -> dict[str, Any]:
    li_abs = np.abs(rows["x_RLi_calc_from_paper_K"] - rows["x_RLi_exp"])
    na_abs = np.abs(rows["x_RNa_calc_from_paper_K"] - rows["x_RNa_exp"])
    li_rel = li_abs / rows["x_RLi_exp"].clip(lower=1.0e-300)
    na_rel = na_abs / rows["x_RNa_exp"].clip(lower=1.0e-300)
    li_lnk_offset = rows["lnQ_Li_at_exp"] - math.log(constants["Li"])
    na_lnk_offset = rows["lnQ_Na_at_exp"] - math.log(constants["Na"])
    native_li = rows["native_cross_phase_residual_Li"]
    native_na = rows["native_cross_phase_residual_Na"]
    return {
        "row_count": int(len(rows)),
        "paper_constants": constants,
        "paper_lnK": {"Li": math.log(constants["Li"]), "Na": math.log(constants["Na"])},
        "median_lnQ_at_experimental_rows": {
            "Li": float(rows["lnQ_Li_at_exp"].median()),
            "Na": float(rows["lnQ_Na_at_exp"].median()),
        },
        "mean_lnQ_at_experimental_rows": {
            "Li": float(rows["lnQ_Li_at_exp"].mean()),
            "Na": float(rows["lnQ_Na_at_exp"].mean()),
        },
        "median_lnQ_minus_lnK": {
            "Li": float(li_lnk_offset.median()),
            "Na": float(na_lnk_offset.median()),
        },
        "package_phase_tagged_cross_phase": {
            "evaluated_rows": int(rows["package_cross_phase_evaluated"].sum()),
            "reaction_phase_scope": "phase_tagged_cross_phase",
            "native_reaction_residual_size": int(rows["native_cross_phase_residual_size"].max()),
            "median_native_cross_phase_residual": {
                "Li": float(native_li.median()),
                "Na": float(native_na.median()),
            },
            "median_abs_native_cross_phase_residual": {
                "Li": float(native_li.abs().median()),
                "Na": float(native_na.abs().median()),
            },
            "max_phase_charge_balance_norm": float(rows["native_phase_charge_balance_norm"].max()),
            "max_element_balance_norm": float(rows["native_element_balance_norm"].max()),
        },
        "max_abs_charge_residual": float(np.max(np.abs(rows["aqueous_charge_residual"]))),
        "median_abs_complex_error_from_paper_K": {
            "RLi": float(li_abs.median()),
            "RNa": float(na_abs.median()),
        },
        "mean_relative_complex_error_from_paper_K": {
            "RLi": float(li_rel.mean()),
            "RNa": float(na_rel.mean()),
        },
        "effective_constants_implied_by_experimental_rows": {
            "Li_median": float(math.exp(rows["lnQ_Li_at_exp"].median())),
            "Na_median": float(math.exp(rows["lnQ_Na_at_exp"].median())),
            "Li_geomean": float(math.exp(rows["lnQ_Li_at_exp"].mean())),
            "Na_geomean": float(math.exp(rows["lnQ_Na_at_exp"].mean())),
        },
        "status": (
            "source_mismatch"
            if float(li_lnk_offset.median()) > 5.0 or float(na_lnk_offset.median()) > 5.0
            else "source_replay_consistent"
        ),
    }


def _write_report(summary: dict[str, Any]) -> None:
    lines = [
        "# Rezaee 2026 Reactive Equilibrium Replay",
        "",
        "## Source Basis",
        "",
        "- 2025 source: `papers/pdf/Rezaee et al. - 2025 - Application of Response Surface Methodology for Selective Extraction of Lithium Using a Hydrophobic DES.pdf`.",
        "- 2025 supporting information: `papers/pdf/Rezaee et al. - 2025 - Supporting information - Application of Response Surface Methodology for Selective Extraction of Lithium.pdf`.",
        "- Main source: `papers/pdf/Rezaee et al. - 2026 - Thermodynamic modeling of lithium extraction from synthetic brine using deep eutectic solvents A PC.pdf`.",
        "- 2026 supporting information: `papers/pdf/Rezaee et al. - 2026 - Supplementary material - Thermodynamic modeling of lithium extraction from synthetic brine using deep eutectic solvents.pdf`.",
        "- Searchable source text: `papers/md/Rezaee et al. - 2026 - Thermodynamic modeling of lithium extraction from synthetic brine using deep eutectic solvents A PC.md`.",
        "- The replay follows the paper's phase-specific reaction-equilibrium formulation rather than a conventional same-species LLE fugacity equality.",
        "- The package check uses `ReactionDefinition.phase_stoichiometry` and the native phase-tagged cross-phase residual block.",
        "- Aqueous phase: H2O, Li+, Na+, Cl-, H+, OH-, NH4+ with ePC-SAFT component activity coefficients.",
        "- Organic phase: DES, TOPO, RLi, RNa with PC-SAFT activity coefficients calculated as mixture fugacity coefficient over pure-component fugacity coefficient.",
        "",
        "## Result",
        "",
        f"- Rows replayed: `{summary['row_count']}`.",
        f"- Status: `{summary['status']}`.",
        f"- Paper lnK Li/Na: `{summary['paper_lnK']['Li']}`, `{summary['paper_lnK']['Na']}`.",
        f"- Median lnQ at experimental rows Li/Na: `{summary['median_lnQ_at_experimental_rows']['Li']}`, `{summary['median_lnQ_at_experimental_rows']['Na']}`.",
        f"- Median lnQ-lnK Li/Na: `{summary['median_lnQ_minus_lnK']['Li']}`, `{summary['median_lnQ_minus_lnK']['Na']}`.",
        f"- Package phase-tagged cross-phase rows evaluated: `{summary['package_phase_tagged_cross_phase']['evaluated_rows']}`.",
        f"- Package native median cross-phase residual Li/Na: `{summary['package_phase_tagged_cross_phase']['median_native_cross_phase_residual']['Li']}`, `{summary['package_phase_tagged_cross_phase']['median_native_cross_phase_residual']['Na']}`.",
        f"- Median absolute RLi/RNa error from paper K replay: `{summary['median_abs_complex_error_from_paper_K']['RLi']}`, `{summary['median_abs_complex_error_from_paper_K']['RNa']}`.",
        f"- Mean relative RLi/RNa error from paper K replay: `{summary['mean_relative_complex_error_from_paper_K']['RLi']}`, `{summary['mean_relative_complex_error_from_paper_K']['RNa']}`.",
        "",
        "## Interpretation",
        "",
        "The package can evaluate the phase-tagged cross-phase reaction residual required by Rezaee's formulation. However, using the paper-reported Table 2 equilibrium constants together with the paper/SI composition rows and the paper-reported organic parameters does not reproduce the reported RLi/RNa complex mole fractions under this activity-reference convention.",
        "",
        "This is not the same as the old four-species fixed-composition LLE smoke. This replay includes the chemical-equilibrium equations that control lithium and sodium extraction. The current blocker is the source/model convention gap exposed by lnQ-lnK, not an omitted call to `electrolyte_lle`.",
        "",
        "Next implementation step: resolve the convention gap by checking the 2026 supporting information/group-contribution worksheet and the 2025 phase-amount basis, then either correct the stored constants/reference-state convention or add a calibrated Rezaee parameter-refit lane that uses these EOS activity calls directly.",
    ]
    REPORT_MD.parent.mkdir(parents=True, exist_ok=True)
    REPORT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    constants = _reaction_constants()
    aqueous_mix, aqueous_charges = _aqueous_mixture()
    organic_mix, _organic_params, pure_ln_phi = _organic_mixture()
    aqueous_params = _aqueous_parameter_payload()
    combined_mixture = _combined_rezaee_mixture(aqueous_params, _organic_params)

    records: list[dict[str, Any]] = []
    for row in pd.read_csv(EQUILIBRIUM_CSV).itertuples(index=False):
        aqueous_x = _aqueous_x(row)
        organic_x_exp = _organic_x(row)
        aqueous_state = aqueous_mix.state(T=TEMPERATURE_K, x=aqueous_x, P=PRESSURE_PA)
        aqueous_gamma = aqueous_state.activity_coefficient(species=AQ_LABELS)
        organic_gamma = _organic_activity_coefficients(organic_mix, pure_ln_phi, organic_x_exp)
        ln_q_li, ln_q_na = _reaction_log_quotients(aqueous_x, organic_x_exp, aqueous_gamma, organic_gamma)
        native_cross_phase_residuals, native_diagnostics = _package_cross_phase_residuals(
            combined_mixture,
            constants,
            aqueous_x,
            organic_x_exp,
        )
        organic_x_calc, iterations, final_delta = _predict_complexes_from_constants(
            aqueous_x,
            organic_x_exp,
            aqueous_gamma,
            organic_mix,
            pure_ln_phi,
            constants,
        )
        records.append(
            {
                "experiment_no": int(row.experiment_no),
                "aqueous_charge_residual": float(np.dot(aqueous_x, aqueous_charges)),
                "x_RLi_exp": float(organic_x_exp[2]),
                "x_RNa_exp": float(organic_x_exp[3]),
                "x_RLi_calc_from_paper_K": float(organic_x_calc[2]),
                "x_RNa_calc_from_paper_K": float(organic_x_calc[3]),
                "x_RLi_abs_error_from_paper_K": float(abs(organic_x_calc[2] - organic_x_exp[2])),
                "x_RNa_abs_error_from_paper_K": float(abs(organic_x_calc[3] - organic_x_exp[3])),
                "lnQ_Li_at_exp": float(ln_q_li),
                "lnQ_Na_at_exp": float(ln_q_na),
                "lnQ_minus_lnK_Li": float(ln_q_li - math.log(constants["Li"])),
                "lnQ_minus_lnK_Na": float(ln_q_na - math.log(constants["Na"])),
                "package_cross_phase_evaluated": bool(native_diagnostics["cross_phase_reaction_residuals"]),
                "native_reaction_phase_scope": str(native_diagnostics["reaction_phase_scope"]),
                "native_cross_phase_residual_size": int(native_diagnostics["cross_phase_reaction_residual_size"]),
                "native_cross_phase_residual_Li": float(native_cross_phase_residuals[0]),
                "native_cross_phase_residual_Na": float(native_cross_phase_residuals[1]),
                "native_phase_charge_balance_norm": float(native_diagnostics["phase_charge_balance_norm"]),
                "native_element_balance_norm": float(native_diagnostics["element_balance_norm"]),
                "gamma_H2O": float(aqueous_gamma["H2O"]),
                "gamma_Li_plus": float(aqueous_gamma["Li+"]),
                "gamma_Na_plus": float(aqueous_gamma["Na+"]),
                "gamma_OH_minus": float(aqueous_gamma["OH-"]),
                "gamma_DES": float(organic_gamma["DES"]),
                "gamma_RLi": float(organic_gamma["RLi"]),
                "gamma_RNa": float(organic_gamma["RNa"]),
                "paper_K_fixed_point_iterations": int(iterations),
                "paper_K_fixed_point_final_delta": float(final_delta),
                "source": row.source,
            }
        )

    replay = pd.DataFrame(records)
    REPLAY_CSV.parent.mkdir(parents=True, exist_ok=True)
    replay.to_csv(REPLAY_CSV, index=False)
    summary = _summaries(replay, constants)
    SUMMARY_JSON.parent.mkdir(parents=True, exist_ok=True)
    SUMMARY_JSON.write_text(json.dumps(_jsonable(summary), indent=2) + "\n", encoding="utf-8")
    _write_report(summary)
    print(json.dumps(_jsonable(summary), indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
