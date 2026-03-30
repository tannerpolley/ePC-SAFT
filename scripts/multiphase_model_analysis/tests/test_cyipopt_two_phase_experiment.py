from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pytest

pytest.importorskip("cyipopt")

REPO_ROOT = Path(__file__).resolve().parents[3]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts._env import require_pcsaft_install

require_pcsaft_install()

from pcsaft.parameters import _resolve_runtime_options
from scripts.multiphase_model_analysis import ascani_case2_dataset_comparison as case2
from scripts.multiphase_model_analysis.cyipopt_two_phase_experiment import solve_two_phase_lle_cyipopt


def _dataset_file(*parts: str) -> Path:
    return REPO_ROOT / "src" / "pcsaft" / "data" / "pcsaft_parameters" / Path(*parts)


def _runtime_to_elec_model(runtime: dict) -> dict:
    return {
        "include_ion_term": True,
        "zeta_mode": 1,
        "include_dipole_term": False,
        "alpha_profile": "constant",
        "alpha0": 0.0,
        "f_solv_rule": 0,
        "rel_perm": {"rule": int(runtime["dielc_rule"]), "differential_mode": int(runtime["dielc_diff_mode"])},
        "include_born_model": False,
        "born_model": {"d_Born_mode": 0, "solvation_shell_model": False, "dielectric_saturation": False},
    }


def _canonical_case_inputs() -> tuple[float, float, list[str], np.ndarray, dict]:
    t = 298.15
    p = 1.0e5
    species = ["H2O-2B-Li", "Na+", "Cl-"]
    canonical = json.loads(_dataset_file("2020_Bulow", "user_options.json").read_text(encoding="utf-8"))
    runtime = _resolve_runtime_options(canonical)["runtime"]
    runtime["dielc_rule"] = 1
    runtime["dielc_diff_mode"] = 0
    s_water = 2.7927 + 10.11 * np.exp(-0.01775 * t) - 1.417 * np.exp(-0.01146 * t)
    params = {
        "MW": np.asarray([18.01528e-3, 22.98e-3, 35.45e-3]),
        "m": np.asarray([1.2047, 1.0, 1.0]),
        "s": np.asarray([s_water, 2.8232, 2.7560]),
        "e": np.asarray([353.95, 230.0, 170.0]),
        "e_assoc": np.asarray([2425.7, 0.0, 0.0]),
        "vol_a": np.asarray([0.04509, 0.0, 0.0]),
        "assoc_scheme": ["2B", None, None],
        "dipm": np.asarray([0.0, 0.0, 0.0]),
        "dip_num": np.asarray([1.0, 1.0, 1.0]),
        "z": np.asarray([0.0, 1.0, -1.0]),
        "dielc": np.asarray([78.09, 8.0, 8.0]),
        "d_born": np.asarray([0.0, 3.445, 4.1]),
        "f_solv": np.asarray([1.5, 1.0, 1.0]),
        "k_ij": np.asarray([[0.0, 0.0045, -0.25], [0.0045, 0.0, 0.317], [-0.25, 0.317, 0.0]]),
        "l_ij": np.zeros((3, 3)),
        "k_hb": np.zeros((3, 3)),
        "elec_model": _runtime_to_elec_model(runtime),
        "debug": False,
    }
    n = np.asarray([1.0 / 0.01801528, 1.0e-4, 1.0e-4])
    z_feed = n / np.sum(n)
    return t, p, species, z_feed, params


def _max_mean_ionic_gap(out: dict) -> float:
    if int(out.get("n_phases", 0)) != 2:
        return 0.0
    phases = out["phases"]
    lnf_delta = np.asarray(phases[0]["lnfug"], dtype=float) - np.asarray(phases[1]["lnfug"], dtype=float)
    charged_idx = np.asarray(out["charged_species_indices"], dtype=int)
    e_matrix = np.asarray(out["e_matrix"], dtype=float)
    residual = e_matrix.dot(lnf_delta[charged_idx]) if charged_idx.size else np.zeros(0, dtype=float)
    return float(np.max(np.abs(residual))) if residual.size else 0.0


def _phase_roles(out: dict, species: list[str]) -> tuple[dict, dict]:
    butanol_name = "Butanol" if "Butanol" in species else species[0]
    idx = species.index(butanol_name)
    ph0, ph1 = out["phases"]
    org = ph0 if ph0["x"][idx] >= ph1["x"][idx] else ph1
    aq = ph1 if org is ph0 else ph0
    return org, aq


def test_cyipopt_experiment_returns_expected_structure_for_canonical_case():
    t, p, species, z_feed, params = _canonical_case_inputs()
    out = solve_two_phase_lle_cyipopt(
        t,
        p,
        z_feed,
        params,
        species,
        options={"tpdf_global_trials": 300, "tpdf_local_trials": 120, "tpdf_tol": -1.0e-6, "solver_tol": 1.0e-8, "max_nfev": 250, "fd_eps": 1.0e-7},
    )
    assert "n_phases" in out
    assert "phases" in out
    assert "e_matrix" in out
    assert "solver_info" in out
    assert "tpdf_seed_x" in out
    if int(out["n_phases"]) == 2:
        x0 = np.asarray(out["phases"][0]["x"], dtype=float)
        x1 = np.asarray(out["phases"][1]["x"], dtype=float)
        beta0 = float(out["phases"][0]["beta"])
        beta1 = float(out["phases"][1]["beta"])
        z = np.asarray(params["z"], dtype=float)
        assert np.max(np.abs(z_feed - (beta0 * x0 + beta1 * x1))) <= 1.0e-8
        assert abs(float(np.dot(z, x0))) <= 1.0e-6
        assert abs(float(np.dot(z, x1))) <= 1.0e-6
        assert np.max(np.abs(x0 - x1)) > 1.0e-3
        assert _max_mean_ionic_gap(out) <= 1.0e-3


def test_cyipopt_experiment_tracks_current_solver_on_ascani_case2():
    config = case2._default_model_configs()[0]
    species, z_feed, _ = case2._case2_feed()
    params = case2._build_params_for_config(config, species, z_feed)
    current = case2._solve_lle_with_retries(case2.T_REF, case2.P_REF, z_feed, params, species)
    assert bool(current.get("converged", False))
    assert int(current.get("n_phases", 0)) == 2
    cy_options = dict(current.get("_solve_options", {}))
    cy_options["fd_eps"] = 1.0e-7
    trial = solve_two_phase_lle_cyipopt(case2.T_REF, case2.P_REF, z_feed, params, species, options=cy_options)
    assert bool(trial.get("converged", False))
    assert int(trial.get("n_phases", 0)) == 2
    trial_org, trial_aq = _phase_roles(trial, species)
    idx_water = species.index("H2O")
    idx_butanol = species.index("Butanol")
    idx_na = species.index("Na+")
    idx_k = species.index("K+")
    idx_cl = species.index("Cl-")
    assert trial_org["x"][idx_butanol] > trial_aq["x"][idx_butanol]
    assert trial_aq["x"][idx_water] > trial_org["x"][idx_water]
    assert trial_aq["x"][idx_na] > trial_org["x"][idx_na]
    assert trial_aq["x"][idx_k] > trial_org["x"][idx_k]
    assert trial_aq["x"][idx_cl] > trial_org["x"][idx_cl]
    assert _max_mean_ionic_gap(trial) <= max(1.0e-3, 10.0 * max(_max_mean_ionic_gap(current), 1.0e-8))
    assert float(trial["residual_norm"]) <= max(1.0, 10.0 * float(current["residual_norm"]))
