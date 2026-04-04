# -*- coding: utf-8 -*-
"""Object-oriented integration tests for the native PC-SAFT runtime."""

import json
from pathlib import Path

import numpy as np
import pytest

import pcsaft
from pcsaft import DATASET_ROOT
from pcsaft import aly_lee
from pcsaft import PCSAFTMixture
from pcsaft import SolutionError
from pcsaft.parameters import _resolve_runtime_options


def _dataset_file(*parts: str) -> Path:
    return DATASET_ROOT.joinpath(*parts)


def _runtime_to_elec_model(runtime):
    radius_to_d_born = {1: 0, 2: 1, 3: 2, 4: 3, 5: 3}
    born_radius_model = int(runtime.get("born_radius_model", 1))
    born_diff_mode = int(runtime.get("born_diff_mode", 0))
    born_model = int(runtime.get("born_model", 1))
    return {
        "rel_perm": {
            "rule": int(runtime.get("dielc_rule", 1)),
            "differential_mode": int(runtime.get("dielc_diff_mode", 0)),
        },
        "hc_model": {
            "dadx_differential_mode": int(runtime.get("hc_dadx_diff_mode", 0)),
        },
        "disp_model": {
            "dadx_differential_mode": int(runtime.get("disp_dadx_diff_mode", 0)),
        },
        "assoc_model": {
            "dadx_differential_mode": int(runtime.get("assoc_dadx_diff_mode", 0)),
        },
        "polar_model": {
            "dadx_differential_mode": int(runtime.get("polar_dadx_diff_mode", 0)),
        },
        "DH_model": {
            "d_ion_mode": int(runtime.get("d_ion_mode", 1)),
            "bjeruum_treatment": bool(runtime.get("bjeruum_treatment", False)),
            "mu_DH_model": {
                "differential_mode": int(runtime.get("mu_DH_diff_mode", 0)),
                "comp_dep_rel_perm": bool(runtime.get("mu_DH_comp_dep_rel_perm", True)),
                "include_sum_term": bool(runtime.get("mu_DH_include_sum_term", True)),
            },
        },
        "include_born_model": bool(runtime.get("include_born_model", born_model != 0)),
        "born_model": {
            "d_Born_mode": int(runtime.get("d_born_mode", radius_to_d_born.get(born_radius_model, 0))),
            "solvation_shell_model": bool(runtime.get("born_solvation_shell_model", born_model == 2)),
            "dielectric_saturation": bool(runtime.get("born_dielectric_saturation", born_model == 2)),
            "bulk_mode": int(runtime.get("born_bulk_mode", runtime.get("born_eps_mode", 0))),
            "mu_born_model": {
                "differential_mode": int(runtime.get("mu_born_diff_mode", 1 if born_diff_mode == 1 else 0)),
                "comp_dep_rel_perm": bool(runtime.get("mu_born_comp_dep_rel_perm", born_diff_mode != 3)),
                "include_sum_term": bool(runtime.get("mu_born_include_sum_term", born_diff_mode != 2)),
                "comp_dep_delta_d": bool(runtime.get("mu_born_comp_dep_delta_d", False)),
            },
        },
    }


def test_package_exports_are_object_first():
    assert hasattr(pcsaft, "PCSAFTMixture")
    assert hasattr(pcsaft, "PCSAFTState")
    assert not hasattr(pcsaft, "pcsaft_den")
    assert not hasattr(pcsaft, "pcsaft_cp")
    assert not hasattr(pcsaft, "pcsaft_ares")
    assert not hasattr(pcsaft, "pcsaft_lnfugcoef")
    assert not hasattr(pcsaft, "pcsaft_dielectric_eval")
    with pytest.raises(ImportError):
        exec("from pcsaft.pcsaft import pcsaft_ares", {})


def test_state_properties_match_regression_values():
    t = 233.15
    rho = 14330.417110
    x = np.array([0.1, 0.3, 0.6])
    params = {
        "m": np.asarray([1.0000, 1.6069, 2.0020]),
        "s": np.asarray([3.7039, 3.5206, 3.6184]),
        "e": np.asarray([150.03, 191.42, 208.11]),
        "k_ij": np.asarray([
            [0.0, 3.0e-4, 1.15e-2],
            [3.0e-4, 0.0, 5.10e-3],
            [1.15e-2, 5.10e-3, 0.0],
        ]),
    }

    mix = PCSAFTMixture.from_params(params)
    state = mix.state(T=t, x=x, rho=rho)

    assert state.density() == pytest.approx(rho)
    assert state.a_res() == pytest.approx(-3.54988543593195, rel=1e-6)
    assert state.pressure() > 0.0
    assert state.Z() > 0.0


def test_state_roundtrip_hres_for_liquid_and_vapor():
    t = 325.0
    p = 101325.0
    x = np.asarray([1.0])

    tol = 1e-2

    # Toluene
    mix = PCSAFTMixture.from_params({
        "m": np.asarray([2.8149]),
        "s": np.asarray([3.7169]),
        "e": np.asarray([285.69]),
    })

    liq = mix.state(T=t, x=x, P=p, phase="liq")
    vap = mix.state(T=t, x=x, P=p, phase="vap")

    assert liq.density() > vap.density()
    assert liq.h_res() == pytest.approx(-36809.39, rel=tol)
    assert vap.h_res() == pytest.approx(-362.6777, rel=tol)


def test_state_methods_and_breakdown_match_existing_accessors():
    t = 233.15
    rho = 14330.417110
    x = np.array([0.1, 0.3, 0.6])
    params = {
        "m": np.asarray([1.0000, 1.6069, 2.0020]),
        "s": np.asarray([3.7039, 3.5206, 3.6184]),
        "e": np.asarray([150.03, 191.42, 208.11]),
        "k_ij": np.asarray([
            [0.0, 3.0e-4, 1.15e-2],
            [3.0e-4, 0.0, 5.10e-3],
            [1.15e-2, 5.10e-3, 0.0],
        ]),
    }

    mix = PCSAFTMixture.from_params(params)
    state = mix.state(T=t, x=x, rho=rho)

    assert not hasattr(state, "ares")
    assert not hasattr(state, "hres")
    assert not hasattr(state, "sres")
    assert not hasattr(state, "gres")
    assert not hasattr(state, "mures")
    assert not hasattr(state, "fugcoef")
    assert not hasattr(state, "dielc_eval")
    assert np.isfinite(state.a_res())
    assert np.isfinite(state.h_res())
    assert np.isfinite(state.s_res())
    assert np.isfinite(state.g_res())
    assert np.all(np.isfinite(state.mu_res()))
    assert np.all(np.isfinite(state.gamma()))

    breakdown = state.breakdown()
    assert breakdown["miac"] == {}
    assert breakdown["miac_m"] == {}
    assert breakdown["gsolv"] == {}
    assert breakdown["dielectric_eval"] is None
    assert breakdown["osmoticC"] is None
    np.testing.assert_allclose(breakdown["a_res"], state.a_res())
    np.testing.assert_allclose(breakdown["h_res"], state.h_res())
    np.testing.assert_allclose(breakdown["s_res"], state.s_res())
    np.testing.assert_allclose(breakdown["g_res"], state.g_res())
    np.testing.assert_allclose(breakdown["mu_res"], state.mu_res())
    np.testing.assert_allclose(breakdown["lnfugcoef"], state.lnfugcoef())
    np.testing.assert_allclose(breakdown["gamma"], state.gamma())
    assert "fugcoef" not in breakdown

    aly = np.asarray([1.2, 0.8, -0.01, 2.0e-5, -3.0e-8], dtype=float)
    cp_mix = PCSAFTMixture.from_params({
        "m": np.asarray([2.8149]),
        "s": np.asarray([3.7169]),
        "e": np.asarray([285.69]),
    })
    cp_state = cp_mix.state(T=325.0, x=np.asarray([1.0]), P=101325.0)
    cp_expected = aly_lee(cp_state.T, aly)
    cp_expected += (
        cp_mix.state(T=cp_state.T + 0.001, x=np.asarray([1.0]), P=101325.0, phase="liq").h_res()
        - cp_mix.state(T=cp_state.T - 0.001, x=np.asarray([1.0]), P=101325.0, phase="liq").h_res()
    ) / 0.002
    assert cp_state.cp(aly) == pytest.approx(cp_expected, rel=1e-8, abs=1e-8)


def test_flash_and_vaporization_return_structured_results():
    t = 325.0
    p = 101325.0
    x = np.asarray([1.0])
    mix = PCSAFTMixture.from_params({
        "m": np.asarray([2.8149]),
        "s": np.asarray([3.7169]),
        "e": np.asarray([285.69]),
    })
    state = mix.state(T=t, x=x, P=p, phase="liq")

    flash = state.flashTQ(q=0.0)
    hvap = state.Hvap()

    assert flash.kind == "TQ"
    assert flash.pressure is not None
    assert len(flash.phases) == 2
    assert flash.phases[0].x.shape == (1,)
    assert hvap.pressure > 0.0
    assert np.isfinite(hvap.value)


def test_state_constructor_rejects_invalid_pressure_and_density():
    mix = PCSAFTMixture.from_params({
        "m": np.asarray([1.0]),
        "s": np.asarray([3.0]),
        "e": np.asarray([150.0]),
    })

    with pytest.raises(pcsaft.InputError):
        mix.state(T=300.0, x=np.asarray([1.0]), P=0.0)

    with pytest.raises(pcsaft.InputError):
        mix.state(T=300.0, x=np.asarray([1.0]), rho=-1.0)


def test_dataset_state_methods_and_lle_workflow():
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
        "k_ij": np.asarray([
            [0.0, 0.0045, -0.25],
            [0.0045, 0.0, 0.317],
            [-0.25, 0.317, 0.0],
        ]),
        "l_ij": np.zeros((3, 3)),
        "k_hb": np.zeros((3, 3)),
        "elec_model": _runtime_to_elec_model(runtime),
        "debug": bool(runtime["debug"]),
    }

    z_feed = np.asarray([1.0 / 0.01801528, 1e-4, 1e-4], dtype=float)
    z_feed /= np.sum(z_feed)

    mix = PCSAFTMixture.from_params(params, species=species)
    state = mix.state(T=t, x=z_feed, P=p)

    eps, deps = state.dielectric_eval()
    assert eps > 0.0
    assert deps.shape == (3,)
    assert state.osmoticC().shape == (1,)
    miac = state.miac(species=species)
    miac_m = state.miac_m(species=species)
    gsolv = state.gsolv(species=species)
    assert list(miac) == ["Na+Cl-"]
    assert list(miac_m) == ["Na+Cl-"]
    assert list(gsolv) == ["Na+", "Cl-"]
    assert np.isfinite(next(iter(miac.values())))
    assert np.isfinite(next(iter(miac_m.values())))
    assert all(np.isfinite(value) for value in gsolv.values())

    breakdown = state.breakdown(species=species)
    assert breakdown["miac"] == miac
    assert breakdown["miac_m"] == miac_m
    assert breakdown["gsolv"] == gsolv
    np.testing.assert_allclose(breakdown["mu_res"], state.mu_res())
    np.testing.assert_allclose(breakdown["gamma"], state.gamma())
    np.testing.assert_allclose(breakdown["lnfugcoef"], state.lnfugcoef())
    np.testing.assert_allclose(breakdown["dielectric_eval"][0], eps)
    np.testing.assert_allclose(breakdown["dielectric_eval"][1], deps)

    lle = state.multiphase_lle(z_feed, species=species, options={"tpdf_global_trials": 50, "tpdf_local_trials": 20, "tpdf_tol": -1e-6})
    assert lle.n_phases >= 1
    assert lle.phases
    assert lle.e_matrix.shape[1] == len(lle.charged_species)
    assert lle.charged_species


def test_missing_public_procedural_api_is_not_exported():
    with pytest.raises(AttributeError):
        getattr(pcsaft, "pcsaft_den")
