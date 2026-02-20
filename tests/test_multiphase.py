# -*- coding: utf-8 -*-
"""
Tests for multiphase electrolyte LLE helper and solver.
"""
import numpy as np
import pcsaft as pcs
from data.epcsaft_properties import get_prop_dict, molality_to_molefraction


def _ascani_case_feed_molefractions():
    # Ascani 2022 case-study 2 feed (mass fraction basis):
    # water + 1-butanol + NaCl + KCl = (0.8094, 0.1728, 0.0054, 0.0124)
    w_water = 0.8094
    w_but = 0.1728
    w_nacl = 0.0054
    w_kcl = 0.0124

    mw_water = 18.01528e-3
    mw_but = 74.12e-3
    mw_nacl = (22.98976928 + 35.453) * 1e-3
    mw_kcl = (39.0983 + 35.453) * 1e-3

    n_water = w_water / mw_water
    n_but = w_but / mw_but
    n_na = w_nacl / mw_nacl
    n_k = w_kcl / mw_kcl
    n_cl = n_na + n_k

    species = ["H2O-2B-Li", "Butanol", "Na+", "K+", "Cl-"]
    n = np.array([n_water, n_but, n_na, n_k, n_cl], dtype=float)
    return species, n / np.sum(n)


def test_build_e_matrix_rank():
    species = ["water", "butanol", "Na+", "K+", "Cl-"]
    z = np.array([0.0, 0.0, 1.0, 1.0, -1.0], dtype=float)
    z_feed = np.array([0.40, 0.40, 0.09, 0.06, 0.15], dtype=float)

    E, info = pcs._build_e_matrix(z, z_feed, species=species)
    assert E.shape == (2, 3)
    assert np.linalg.matrix_rank(E) == 2
    assert len(info["ion_pair_rows"]) == 2
    assert info["charged_species"] == ["Na+", "K+", "Cl-"]


def test_multiphase_ascani_case_smoke():
    t = 298.15
    p = 1.0e5
    species, z_feed = _ascani_case_feed_molefractions()
    params = get_prop_dict(species, z_feed, t, user_options={"elec_model": "2020", "debug": False, 'dielc_rule': 3})

    result = pcs.pcsaft_multiphase_lle(
        t, p, z_feed, params, species,
        options={
            "tpdf_global_trials": 1200,
            "tpdf_local_trials": 600,
            "solver_tol": 1e-9,
            "max_nfev": 220,
            "debug": False,
        },
    )

    assert result["converged"]
    assert result["n_phases"] == 2
    assert len(result["phases"]) == 2

    z = np.asarray(params["z"], dtype=float)
    x1 = np.asarray(result["phases"][0]["x"], dtype=float)
    x2 = np.asarray(result["phases"][1]["x"], dtype=float)
    beta = float(result["phases"][0]["beta"])

    assert abs(np.sum(x1) - 1.0) < 1e-10
    assert abs(np.sum(x2) - 1.0) < 1e-10
    assert abs(np.dot(z, x1)) < 1e-6
    assert abs(np.dot(z, x2)) < 1e-6
    assert np.max(np.abs(z_feed - (beta*x1 + (1.0 - beta)*x2))) < 1e-8

    # Qualitative split: butanol-rich and aqueous-rich phases are distinct.
    idx_but = species.index("Butanol")
    assert abs(x1[idx_but] - x2[idx_but]) > 0.01

    # Check residual system consistency from returned states.
    neutral_idx = np.where(np.abs(z) <= 1e-12)[0]
    charged_idx = np.where(np.abs(z) > 1e-12)[0]
    E = np.asarray(result["e_matrix"], dtype=float)
    dlnf = np.asarray(result["phases"][0]["lnfug"]) - np.asarray(result["phases"][1]["lnfug"])
    res = np.concatenate([dlnf[neutral_idx], E.dot(dlnf[charged_idx]), np.array([np.dot(z, x1)])])
    assert np.linalg.norm(res) < 2e-1


def test_multiphase_stable_dilute_returns_one_phase():
    t = 298.15
    p = 1.0e5
    species = ["H2O-2B-Li", "Na+", "Cl-"]
    z_feed = molality_to_molefraction(1e-4, species=species, solvent="H2O-2B-Li")
    params = get_prop_dict(species, z_feed, t, user_options={"elec_model": "2020", "debug": False})

    result = pcs.pcsaft_multiphase_lle(
        t, p, z_feed, params, species,
        options={
            "tpdf_global_trials": 500,
            "tpdf_local_trials": 200,
            "tpdf_tol": -1e-6,
            "solver_tol": 1e-9,
            "max_nfev": 150,
        },
    )

    assert result["converged"]
    assert result["n_phases"] == 1
    assert len(result["phases"]) == 1
