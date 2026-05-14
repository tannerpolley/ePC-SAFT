from __future__ import annotations

import csv
import sys
from pathlib import Path

import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
FIGURE_DIR = SCRIPT_DIR.parent
ANALYSIS_ROOT = FIGURE_DIR.parents[2]
REPO_ROOT = ANALYSIS_ROOT.parents[1]

if str(ANALYSIS_ROOT) not in sys.path:
    sys.path.insert(0, str(ANALYSIS_ROOT))
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts._env import require_epcsaft_install

require_epcsaft_install()

import _plot_common as common
import _model_overlay as overlay
from epcsaft.parameters import get_prop_dict
from scripts._epcsaft_oop import epcsaft_density, epcsaft_fugacity_coefficient_terms, epcsaft_pressure

DATA_PATH = common.analysis_data_path(FIGURE_DIR, "water_contributions.csv", kind="processed", category="figure_3")
OUTPUT_CSV = common.analysis_runs_path(
    __file__, "figure3_mu_contribution_numerical_comparison.csv", category=("figure_3", "diagnostics")
)
T_REF = overlay.T_REF
P_REF = overlay.P_REF
EPS = overlay.EPS
EPS_INF = overlay.EPS_INF
RT_KJMOL = overlay.R_GAS * T_REF / 1000.0

CONTRIBUTION_MAP = {
    "hc": {"paper_rows": ("hc avg", "hc"), "suffix": "hc"},
    "disp": {"paper_rows": ("disp avg", "disp"), "suffix": "disp"},
    "assoc": {"paper_rows": ("assoc avg", "assoc"), "suffix": "assoc"},
    "dh": {"paper_rows": (), "suffix": "ion"},
    "born": {"paper_rows": ("born avg", "born"), "suffix": "born"},
}

NUMERICAL_USER_OPTIONS = {
    "elec_model": {
        "hc_model": {"dadx_differential_mode": "numerical"},
        "disp_model": {"dadx_differential_mode": "numerical"},
        "assoc_model": {"dadx_differential_mode": "numerical"},
        "polar_model": {"dadx_differential_mode": "numerical"},
    }
}


def _paper_value(frame: common.Table, contribution: str, ion: str) -> float:
    if contribution == "dh":
        return 0.0
    for row_key in CONTRIBUTION_MAP[contribution]["paper_rows"]:
        if row_key in frame.index:
            return frame.scalar(row_key, ion)
    raise KeyError(f"Missing paper row for contribution '{contribution}'.")


def _kjmol(value: float) -> float:
    return float(RT_KJMOL * float(value))


def _terms_for_user_options(ion: str, user_options: dict | None) -> tuple[dict[str, object], int]:
    species = overlay._species_for_ion(ion, "water")
    x = np.asarray([EPS, EPS, 1.0 - 2.0 * EPS], dtype=float)
    params = get_prop_dict("2020_Bulow", species, x, T_REF, user_options=user_options or {})
    rho = epcsaft_density(T_REF, P_REF, x, params, phase="liq")

    z = np.asarray(params.get("z", []), dtype=float)
    idx_ion = np.where(np.abs(z) > 1.0e-12)[0]
    idx_solv = np.where(np.abs(z) <= 1.0e-12)[0]
    x_ref = x.copy()
    x_ref[idx_ion] = 0.0
    solv_sum = float(np.sum(x_ref[idx_solv]))
    if solv_sum > 0.0:
        x_ref[idx_solv] /= solv_sum
    else:
        x_ref[idx_solv] = 1.0 / len(idx_solv)

    p_ref = epcsaft_pressure(T_REF, rho, x_ref, params)
    x_inf = x_ref.copy()
    ion_idx = species.index(ion)
    x_inf[ion_idx] = EPS_INF
    x_inf /= np.sum(x_inf)
    phase = "vap" if rho < 900.0 else "liq"
    rho_inf = epcsaft_density(T_REF, p_ref, x_inf, params, phase=phase)
    return epcsaft_fugacity_coefficient_terms(T_REF, rho_inf, x_inf, params), ion_idx


def _extract_terms(terms: dict[str, object], idx: int, suffix: str) -> dict[str, float]:
    mu_term = float(np.asarray(terms[f"mu_{suffix}"], dtype=float)[idx])
    a_term = float(terms[f"a_{suffix}"])
    z_term = float(terms[f"z_{suffix}"])
    dadx_term = float(np.asarray(terms[f"dadx_{suffix}"], dtype=float)[idx])
    sum_x_dadx_term = float(terms[f"sum_x_dadx_{suffix}"])
    return {
        "mu": _kjmol(mu_term),
        "a": _kjmol(a_term),
        "z": _kjmol(z_term),
        "dadx": _kjmol(dadx_term),
        "sum_xj_dadx": _kjmol(sum_x_dadx_term),
        "manual": _kjmol(a_term + z_term + dadx_term - sum_x_dadx_term),
    }


def _build_rows() -> list[dict[str, object]]:
    frame = common.load_indexed_csv(DATA_PATH)
    ions = list(frame.columns)
    rows: list[dict[str, object]] = []
    contribution_order = {key: idx for idx, key in enumerate(CONTRIBUTION_MAP)}
    ion_order = {ion: idx for idx, ion in enumerate(ions)}

    for contribution, meta in CONTRIBUTION_MAP.items():
        for ion in ions:
            paper_mu = _paper_value(frame, contribution, ion)
            analytical_terms, analytical_idx = _terms_for_user_options(ion, None)
            numerical_terms, numerical_idx = _terms_for_user_options(ion, NUMERICAL_USER_OPTIONS)
            analytical = _extract_terms(analytical_terms, analytical_idx, str(meta["suffix"]))
            numerical = _extract_terms(numerical_terms, numerical_idx, str(meta["suffix"]))
            rows.append(
                {
                    "ion": ion,
                    "contr": contribution,
                    "paper_mu_contr": paper_mu,
                    "epcsaft_mu_contr_analytical": analytical["mu"],
                    "epcsaft_mu_manual_sum_analytical": analytical["manual"],
                    "a_contr_analytical": analytical["a"],
                    "z_contr_analytical": analytical["z"],
                    "dadx_contr_analytical": analytical["dadx"],
                    "sum_xj_dadx_contr_analytical": analytical["sum_xj_dadx"],
                    "epcsaft_mu_contr_numerical": numerical["mu"],
                    "epcsaft_mu_manual_sum_numerical": numerical["manual"],
                    "a_contr_numerical": numerical["a"],
                    "z_contr_numerical": numerical["z"],
                    "dadx_contr_numerical": numerical["dadx"],
                    "sum_xj_dadx_contr_numerical": numerical["sum_xj_dadx"],
                    "num_minus_analytical": numerical["mu"] - analytical["mu"],
                    "paper_minus_analytical": paper_mu - analytical["mu"],
                    "paper_minus_numerical": paper_mu - numerical["mu"],
                }
            )

    rows.sort(key=lambda row: (contribution_order[str(row["contr"])], ion_order[str(row["ion"])]))
    return rows


def main() -> None:
    rows = _build_rows()
    OUTPUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    with OUTPUT_CSV.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "ion",
                "contr",
                "paper_mu_contr",
                "epcsaft_mu_contr_analytical",
                "epcsaft_mu_manual_sum_analytical",
                "a_contr_analytical",
                "z_contr_analytical",
                "dadx_contr_analytical",
                "sum_xj_dadx_contr_analytical",
                "epcsaft_mu_contr_numerical",
                "epcsaft_mu_manual_sum_numerical",
                "a_contr_numerical",
                "z_contr_numerical",
                "dadx_contr_numerical",
                "sum_xj_dadx_contr_numerical",
                "num_minus_analytical",
                "paper_minus_analytical",
                "paper_minus_numerical",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {OUTPUT_CSV}", flush=True)


if __name__ == "__main__":
    main()

