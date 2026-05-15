from __future__ import annotations

import json
import math
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import _paths  # noqa: F401,E402
import rezaee_paper_basis_reaction_coordinate as paper_basis  # noqa: E402
import rezaee_reactive_epcsaft_option_scan as option_scan  # noqa: E402
import rezaee_reactive_equilibrium_replay as replay  # noqa: E402
from epcsaft import ePCSAFTMixture  # noqa: E402

ANALYSIS_DIR = Path(__file__).resolve().parents[1]
INPUT_DIR = ANALYSIS_DIR / "data" / "input"
PROCESSED_DIR = ANALYSIS_DIR / "data" / "processed"
RESULTS_DIR = ANALYSIS_DIR / "results" / "reaction_equilibrium"

DOE_CSV = INPUT_DIR / "rezaee_2025_doe_extraction_responses.csv"
ROWS_CSV = PROCESSED_DIR / "rezaee_2026_section32_equilibrium_replication_rows.csv"
SUMMARY_JSON = RESULTS_DIR / "rezaee_2026_section32_equilibrium_replication_summary.json"
REPORT_MD = RESULTS_DIR / "rezaee_2026_section32_equilibrium_replication.md"

PRESSURE_PA = 101325.0
TEMPERATURE_K = 298.15
EXPERIMENT_COUNT = 26


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


def _organic_mixture(*, use_table9_kij: bool) -> tuple[ePCSAFTMixture, np.ndarray]:
    mix, params, pure_ln_phi = replay._organic_mixture()
    if use_table9_kij:
        return mix, pure_ln_phi

    no_kij = dict(params)
    no_kij["k_ij"] = np.zeros_like(np.asarray(params["k_ij"], dtype=float))
    mix = ePCSAFTMixture.from_params(no_kij, species=replay.ORG_LABELS)
    pure_ln_phi = []
    for i, label in enumerate(replay.ORG_LABELS):
        pure_params: dict[str, Any] = {}
        for key, value in no_kij.items():
            if key == "assoc_scheme":
                pure_params[key] = [value[i]]
            elif key == "k_ij":
                pure_params[key] = np.zeros((1, 1), dtype=float)
            else:
                pure_params[key] = np.asarray([value[i]], dtype=float)
        pure_mix = ePCSAFTMixture.from_params(pure_params, species=[label])
        pure_state = pure_mix.state(T=TEMPERATURE_K, x=np.asarray([1.0]), P=PRESSURE_PA)
        pure_ln_phi.append(float(pure_state.fugacity_coefficient()[0]))
    return mix, np.asarray(pure_ln_phi, dtype=float)


def _activity_coefficients(
    aqueous_mix: ePCSAFTMixture | None,
    organic_mix: ePCSAFTMixture,
    pure_ln_phi: np.ndarray,
    aqueous_x: np.ndarray,
    organic_x: np.ndarray,
    *,
    use_ideal_aqueous: bool,
) -> tuple[dict[str, float], dict[str, float]]:
    if use_ideal_aqueous:
        aqueous_gamma = {label: 1.0 for label in replay.AQ_LABELS}
    else:
        if aqueous_mix is None:
            raise ValueError("aqueous_mix is required unless use_ideal_aqueous=True")
        aqueous_gamma = aqueous_mix.state(
            T=TEMPERATURE_K,
            x=aqueous_x,
            P=PRESSURE_PA,
        ).activity_coefficient(species=replay.AQ_LABELS)
    organic_gamma = replay._organic_activity_coefficients(organic_mix, pure_ln_phi, organic_x)
    return aqueous_gamma, organic_gamma


def _paper_complex_mole_fractions(
    constants: dict[str, float],
    aqueous_x: np.ndarray,
    organic_x: np.ndarray,
    aqueous_gamma: dict[str, float],
    organic_gamma: dict[str, float],
) -> tuple[float, float]:
    x_rli = (
        constants["Li"]
        * (aqueous_x[1] * aqueous_gamma["Li+"])
        * (aqueous_x[5] * aqueous_gamma["OH-"])
        * (organic_x[0] * organic_gamma["DES"])
        / ((organic_gamma["RLi"]) * (aqueous_x[0] * aqueous_gamma["H2O"]))
    )
    x_rna = (
        constants["Na"]
        * (aqueous_x[2] * aqueous_gamma["Na+"])
        * (aqueous_x[5] * aqueous_gamma["OH-"])
        * (organic_x[0] * organic_gamma["DES"])
        / ((organic_gamma["RNa"]) * (aqueous_x[0] * aqueous_gamma["H2O"]))
    )
    return float(max(x_rli, 0.0)), float(max(x_rna, 0.0))


def _solve_section32_row(
    row: Any,
    constants: dict[str, float],
    aqueous_mix: ePCSAFTMixture | None,
    organic_mix: ePCSAFTMixture,
    pure_ln_phi: np.ndarray,
    *,
    case_id: str,
    use_ideal_aqueous: bool,
    oh_inventory_multiplier: float,
    invert_reaction_constants: bool,
    max_iter: int = 120,
    tol: float = 1.0e-12,
) -> dict[str, Any]:
    constants = dict(constants)
    if invert_reaction_constants:
        constants = {"Li": 1.0 / constants["Li"], "Na": 1.0 / constants["Na"]}
    initial = paper_basis._initial_moles(row)
    initial["OH-"] *= oh_inventory_multiplier
    shared_limit = max(min(initial["OH-"], initial["DES"], initial["Li+"] + initial["Na+"]), 1.0e-14)
    xi_li = min(1.0e-12, 0.25 * shared_limit, 0.01 * initial["Li+"])
    xi_na = min(1.0e-12, 0.25 * max(shared_limit - xi_li, 1.0e-14), 0.01 * initial["Na+"])
    converged = False
    max_delta = math.inf
    organic_total_initial = initial["DES"] + initial["TOPO"] + initial["RLi"] + initial["RNa"]

    for iteration in range(1, max_iter + 1):
        moles = paper_basis._extent_to_moles(initial, xi_li, xi_na)
        if min(moles.values()) <= 0.0:
            xi_li *= 0.5
            xi_na *= 0.5
            continue
        aqueous_x = paper_basis._phase_x(moles, replay.AQ_LABELS)
        organic_x = paper_basis._phase_x(moles, replay.ORG_LABELS)
        aqueous_gamma, organic_gamma = _activity_coefficients(
            aqueous_mix,
            organic_mix,
            pure_ln_phi,
            aqueous_x,
            organic_x,
            use_ideal_aqueous=use_ideal_aqueous,
        )
        x_rli, x_rna = _paper_complex_mole_fractions(constants, aqueous_x, organic_x, aqueous_gamma, organic_gamma)
        next_li = min(max(x_rli * organic_total_initial, 0.0), initial["Li+"], initial["DES"], initial["OH-"])
        next_na = min(
            max(x_rna * organic_total_initial, 0.0),
            initial["Na+"],
            max(initial["DES"] - next_li, 0.0),
            max(initial["OH-"] - next_li, 0.0),
        )
        max_delta = max(abs(next_li - xi_li), abs(next_na - xi_na))
        xi_li = 0.65 * xi_li + 0.35 * next_li
        xi_na = 0.65 * xi_na + 0.35 * next_na
        if max_delta <= tol:
            converged = True
            break

    li_extraction = 100.0 * xi_li / initial["Li+"]
    na_extraction = 100.0 * xi_na / initial["Na+"]
    selectivity = li_extraction / max(na_extraction, 1.0e-300)
    return {
        "case_id": case_id,
        "experiment_no": int(row.experiment_no),
        "T_C": float(row.T_C),
        "pH": float(row.pH),
        "topo_wt_pct": float(row.topo_wt_pct),
        "na_li_mass_ratio": float(row.na_li_mass_ratio),
        "oh_inventory_multiplier": float(oh_inventory_multiplier),
        "invert_reaction_constants": bool(invert_reaction_constants),
        "initial_Li_mol": float(initial["Li+"]),
        "initial_Na_mol": float(initial["Na+"]),
        "initial_OH_mol": float(initial["OH-"]),
        "Li_extraction_pct_calc": float(li_extraction),
        "Na_extraction_pct_calc": float(na_extraction),
        "selectivity_calc": float(selectivity),
        "Li_extraction_pct_exp": float(row.li_extraction_pct_exp),
        "Na_extraction_pct_exp": float(row.na_extraction_pct_exp),
        "selectivity_exp": float(row.li_na_selectivity_exp),
        "Li_extraction_AARD_contribution_pct": float(
            100.0 * abs(li_extraction - float(row.li_extraction_pct_exp)) / float(row.li_extraction_pct_exp)
        ),
        "selectivity_AARD_contribution_pct": float(
            100.0 * abs(selectivity - float(row.li_na_selectivity_exp)) / float(row.li_na_selectivity_exp)
        ),
        "converged": bool(converged),
        "iterations": int(iteration),
        "final_extent_delta": float(max_delta),
        "source": row.source,
    }


def _build_cases() -> list[dict[str, Any]]:
    multipliers = [1.0, 1.0e3, 1.0e5, 1.0e6, 1.0e7]
    cases: list[dict[str, Any]] = [
        {
            "case_id": "ideal_aqueous_table9_kij_pH_stoich",
            "option_id": None,
            "use_ideal_aqueous": True,
            "use_table9_kij": True,
            "oh_inventory_multiplier": 1.0,
            "invert_reaction_constants": False,
        },
        {
            "case_id": "held_2014_s2_no_born_table9_kij_pH_stoich",
            "option_id": "held_2014_s2_no_born_constant_eps",
            "use_ideal_aqueous": False,
            "use_table9_kij": True,
            "oh_inventory_multiplier": 1.0,
            "invert_reaction_constants": False,
        },
        {
            "case_id": "held_2014_s2_no_born_no_kij_pH_stoich",
            "option_id": "held_2014_s2_no_born_constant_eps",
            "use_ideal_aqueous": False,
            "use_table9_kij": False,
            "oh_inventory_multiplier": 1.0,
            "invert_reaction_constants": False,
        },
        {
            "case_id": "held_2014_s2_no_born_table9_kij_inverseK_diagnostic",
            "option_id": "held_2014_s2_no_born_constant_eps",
            "use_ideal_aqueous": False,
            "use_table9_kij": True,
            "oh_inventory_multiplier": 1.0,
            "invert_reaction_constants": True,
        },
    ]
    for multiplier in multipliers[1:]:
        cases.append(
            {
                "case_id": f"held_2014_s2_no_born_table9_kij_ohx{multiplier:g}",
                "option_id": "held_2014_s2_no_born_constant_eps",
                "use_ideal_aqueous": False,
                "use_table9_kij": True,
                "oh_inventory_multiplier": multiplier,
                "invert_reaction_constants": False,
            }
        )
    return cases


def _run_case(case: dict[str, Any], rows: pd.DataFrame, constants: dict[str, float]) -> list[dict[str, Any]]:
    aqueous_mix = None
    if not case["use_ideal_aqueous"]:
        option = next(item for item in option_scan.OPTION_SETS if item["option_id"] == case["option_id"])
        aqueous_mix, _resolved, _snapshot = option_scan._build_aqueous_mixture(option)
    organic_mix, pure_ln_phi = _organic_mixture(use_table9_kij=bool(case["use_table9_kij"]))
    return [
        _solve_section32_row(
            row,
            constants,
            aqueous_mix,
            organic_mix,
            pure_ln_phi,
            case_id=str(case["case_id"]),
            use_ideal_aqueous=bool(case["use_ideal_aqueous"]),
            oh_inventory_multiplier=float(case["oh_inventory_multiplier"]),
            invert_reaction_constants=bool(case["invert_reaction_constants"]),
        )
        for row in rows.itertuples(index=False)
    ]


def _summarize(rows: pd.DataFrame) -> dict[str, Any]:
    grouped = rows.groupby("case_id", as_index=False).agg(
        row_count=("experiment_no", "count"),
        converged_rows=("converged", "sum"),
        Li_extraction_AARD_pct=("Li_extraction_AARD_contribution_pct", "mean"),
        selectivity_AARD_pct=("selectivity_AARD_contribution_pct", "mean"),
        median_Li_extraction_pct_calc=("Li_extraction_pct_calc", "median"),
        median_Na_extraction_pct_calc=("Na_extraction_pct_calc", "median"),
        median_selectivity_calc=("selectivity_calc", "median"),
        median_initial_OH_mol=("initial_OH_mol", "median"),
    )
    best = grouped.sort_values(["Li_extraction_AARD_pct", "selectivity_AARD_pct"]).iloc[0].to_dict()
    direct = grouped.loc[grouped["case_id"].eq("held_2014_s2_no_born_table9_kij_pH_stoich")]
    direct_row = direct.iloc[0].to_dict() if not direct.empty else {}
    return {
        "status": "section32_equation_replication_ran",
        "row_count": int(rows["experiment_no"].nunique()),
        "intended_benchmark_row_count": EXPERIMENT_COUNT,
        "paper_reference_AARD_pct": {
            "before_table9_kij": {"Li_extraction": 15.11, "selectivity": 16.73},
            "after_table9_kij": {"Li_extraction": 7.89, "selectivity": 8.63},
        },
        "direct_held2014_table9_pH_stoich": direct_row,
        "best_case": best,
        "case_summary": grouped.to_dict(orient="records"),
        "interpretation": (
            "The script follows Rezaee Section 3.2 after Table 8: Eq. 17 updates reaction-coordinate "
            "mole fractions, ePC-SAFT/PC-SAFT supplies activity coefficients for Eqs. 14 and 15, and "
            "Eqs. 18-20 evaluate Li extraction, Li/Na selectivity, and AARD over the 26 benchmark rows. "
            "The direct Held-2014-like no-Born pH-stoichiometric basis is the first source-aligned case."
        ),
    }


def _write_report(summary: dict[str, Any]) -> None:
    best = summary["best_case"]
    direct = summary["direct_held2014_table9_pH_stoich"]
    lines = [
        "# Rezaee 2026 Section 3.2 Equilibrium Replication",
        "",
        "## Scope",
        "",
        "This script starts at the paper text immediately after Table 8: after all PC-SAFT parameters are determined, reaction equilibrium is calculated to predict equilibrium composition, lithium extraction percentage, and Li/Na selectivity.",
        "",
        "The validation uses the paper equations directly rather than package-owned equilibrium solvers:",
        "",
        "- Eq. 17 supplies mole fractions from reaction coordinates.",
        "- Eqs. 14 and 15 update RLi and RNa mole fractions from activity coefficients.",
        "- Eqs. 18 and 19 compute Li extraction and selectivity.",
        "- Eq. 20 computes AARD over the 26 benchmark rows.",
        "",
        "The package is used only to compute activity coefficients. The first source-aligned package option is Held-2014-S2-like ePC-SAFT with no Born term and constant dielectric behavior for the aqueous phase, plus Table 9 organic binary interactions.",
        "",
        "## Result",
        "",
        f"- Rows evaluated: `{summary['row_count']}`.",
        "- Paper reference AARD after Table 9 kij: `7.89%` for lithium extraction and `8.63%` for selectivity.",
        f"- Direct Held-2014/Table-9/pH-stoichiometric Li extraction AARD: `{direct.get('Li_extraction_AARD_pct')}`.",
        f"- Direct Held-2014/Table-9/pH-stoichiometric selectivity AARD: `{direct.get('selectivity_AARD_pct')}`.",
        f"- Best diagnostic case: `{best['case_id']}`.",
        f"- Best diagnostic Li extraction AARD: `{best['Li_extraction_AARD_pct']}`.",
        f"- Best diagnostic selectivity AARD: `{best['selectivity_AARD_pct']}`.",
        f"- Best diagnostic median calculated Li extraction: `{best['median_Li_extraction_pct_calc']}`.",
        f"- Best diagnostic median calculated selectivity: `{best['median_selectivity_calc']}`.",
        "",
        "## Interpretation",
        "",
        summary["interpretation"],
        "",
        "If the direct pH-stoichiometric case does not approach the paper AARD, the next source task is not a package phase-equilibrium solve. It is identifying the paper's exact initial mole basis for OH-/NH4OH/base inventory behind Eq. 17, because Section 3.2 says those initial moles are inputs but the main-text DOE table only exposes pH.",
        "",
        "## Generated Files",
        "",
        f"- `{ROWS_CSV.relative_to(ANALYSIS_DIR)}`",
        f"- `{SUMMARY_JSON.relative_to(ANALYSIS_DIR)}`",
    ]
    REPORT_MD.parent.mkdir(parents=True, exist_ok=True)
    REPORT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    rows = pd.read_csv(DOE_CSV)
    if rows["experiment_no"].nunique() != EXPERIMENT_COUNT:
        raise ValueError(f"Expected {EXPERIMENT_COUNT} benchmark rows, found {rows['experiment_no'].nunique()}")
    constants = replay._reaction_constants()
    records: list[dict[str, Any]] = []
    for case in _build_cases():
        records.extend(_run_case(case, rows, constants))
    out = pd.DataFrame(records)
    ROWS_CSV.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(ROWS_CSV, index=False)
    summary = _summarize(out)
    SUMMARY_JSON.parent.mkdir(parents=True, exist_ok=True)
    SUMMARY_JSON.write_text(json.dumps(_jsonable(summary), indent=2) + "\n", encoding="utf-8")
    _write_report(summary)
    print(json.dumps(_jsonable(summary), indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
