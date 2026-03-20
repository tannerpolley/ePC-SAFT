from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

REPO_ROOT = Path(__file__).resolve().parents[3]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts._env import require_pcsaft_install

require_pcsaft_install()

from pcsaft.parameters import get_prop_dict
from pcsaft import flashTQ, pcsaft_multiphase_lle
from scripts.multiphase_model_analysis import ascani_case2_dataset_comparison as case2


def _max_mean_ionic_gap(out: dict) -> float:
    phases = out["phases"]
    x0 = np.asarray(phases[0]["x"], dtype=float)
    x1 = np.asarray(phases[1]["x"], dtype=float)
    if x0[1] >= x1[1]:
        org = phases[0]
        aq = phases[1]
    else:
        org = phases[1]
        aq = phases[0]
    lnf_delta = np.asarray(org["lnfug"], dtype=float) - np.asarray(aq["lnfug"], dtype=float)
    charged_idx = np.asarray(out["charged_species_indices"], dtype=int)
    e_matrix = np.asarray(out["e_matrix"], dtype=float)
    residual = e_matrix.dot(lnf_delta[charged_idx])
    return float(np.max(np.abs(residual))) if residual.size else 0.0


def _solve_binary_lle(species: list[str], dataset: str, feed: np.ndarray, t: float = 298.15) -> dict:
    params = get_prop_dict(dataset, species, feed, t, user_options={})
    params["debug"] = False
    attempt_options = [
        {
            "tpdf_global_trials": 500,
            "tpdf_local_trials": 200,
            "tpdf_tol": -1e-6,
            "solver_tol": 1e-9,
            "max_nfev": 300,
        },
        {
            "tpdf_global_trials": 1500,
            "tpdf_local_trials": 600,
            "tpdf_tol": -1e-6,
            "solver_tol": 1e-10,
            "max_nfev": 800,
            "charge_weight": 5000.0,
            "solver_accept_norm": 0.5,
            "split_tol": 1e-4,
        },
    ]
    last = None
    for opt in attempt_options:
        out = pcsaft_multiphase_lle(t, 1.0e5, feed, params, species, options=opt)
        last = out
        if bool(out.get("converged", False)) and int(out.get("n_phases", 0)) == 2:
            return out
    if last is None:
        raise RuntimeError("No binary LLE solve was executed.")
    return last


def test_ascani_case2_current_option_sets_are_equilibrated():
    results = [case2._solve_dataset(config) for config in case2._default_model_configs()]
    assert len(results) == 2

    for result in results:
        diag = result["paper_compare"]
        org = result["phases"]["organic"]["x"]
        aq = result["phases"]["aqueous"]["x"]

        assert diag["ghat_delta_j_per_mol"] < 0.0
        assert diag["tpdf_min"] < 0.0
        assert diag["mass_balance_max"] < 1.0e-10
        assert abs(diag["phase_charge_org"]) < 1.0e-6
        assert abs(diag["phase_charge_aq"]) < 1.0e-6
        assert diag["neutral_gap_max"] < 6.0e-2
        assert diag["mean_ionic_gap_max"] < 5.0e-4
        assert diag["e_matrix_rank"] == pytest.approx(diag["mean_ionic_pair_count"], rel=0.0, abs=1.0e-12)

        assert org["Butanol"] > aq["Butanol"]
        assert aq["H2O"] > org["H2O"]
        assert diag["eta_na_to_aq_pct"] > 80.0
        assert diag["eta_k_to_aq_pct"] > 80.0
        assert diag["eta_cl_to_aq_pct"] > 80.0


def test_workbook_vle_methane_ethane_propane_example():
    t = 233.15
    x = np.asarray([0.1, 0.3, 0.6], dtype=float)
    m = np.asarray([1.0000, 1.6069, 2.0020], dtype=float)
    s = np.asarray([3.7039, 3.5206, 3.6184], dtype=float)
    e = np.asarray([150.03, 191.42, 208.11], dtype=float)
    k_ij = np.asarray(
        [
            [0.0, 3.0e-4, 1.15e-2],
            [3.0e-4, 0.0, 5.10e-3],
            [1.15e-2, 5.10e-3, 0.0],
        ],
        dtype=float,
    )
    params = {"m": m, "s": s, "e": e, "k_ij": k_ij}

    pressure, xl, xv = flashTQ(t, 0, x, params)
    expected_pressure = 1276369.47358564
    expected_y = np.asarray([0.7246577940944616, 0.20293469102909495, 0.0724075148764435], dtype=float)

    assert pressure == pytest.approx(expected_pressure, rel=2.0e-5, abs=25.0)
    assert np.allclose(np.asarray(xl, dtype=float), x, atol=1.0e-12, rtol=0.0)
    assert np.allclose(np.asarray(xv, dtype=float), expected_y, atol=5.0e-7, rtol=5.0e-7)


def test_lle_water_c4mim_ntf2_example():
    species = ["H2O", "C4mim+", "NTf2-"]
    feed = np.asarray([0.5, 0.25, 0.25], dtype=float)
    out = _solve_binary_lle(species, "bulow_2019", feed)

    assert bool(out.get("converged", False))
    assert int(out.get("n_phases", 0)) == 2
    assert float(out.get("tpdf_min", 0.0)) < 0.0
    assert float(out.get("residual_norm", np.inf)) < 1.0e-8
    assert _max_mean_ionic_gap(out) < 1.0e-8

    x0 = np.asarray(out["phases"][0]["x"], dtype=float)
    x1 = np.asarray(out["phases"][1]["x"], dtype=float)
    if x0[0] >= x1[0]:
        water_rich, il_rich = x0, x1
        beta_water, beta_il = float(out["phases"][0]["beta"]), float(out["phases"][1]["beta"])
    else:
        water_rich, il_rich = x1, x0
        beta_water, beta_il = float(out["phases"][1]["beta"]), float(out["phases"][0]["beta"])

    assert water_rich[0] > 0.89
    assert il_rich[0] < 1.0e-8
    assert water_rich[1] == pytest.approx(water_rich[2], rel=0.0, abs=1.0e-12)
    assert il_rich[1] == pytest.approx(il_rich[2], rel=0.0, abs=1.0e-12)
    assert beta_il == pytest.approx(0.44235402480674835, rel=1.0e-6, abs=1.0e-9)
    assert beta_water == pytest.approx(0.5576459751932517, rel=1.0e-6, abs=1.0e-9)


