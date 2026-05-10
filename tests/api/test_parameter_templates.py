# -*- coding: utf-8 -*-
"""Tests for user-owned parameter template generation."""

import json

import numpy as np
import pytest

from epcsaft import create_parameter_template
from epcsaft import ePCSAFTMixture
from epcsaft.parameters import _resolve_runtime_options
from epcsaft.parameters import minimize_user_options


def test_create_parameter_template_creates_loadable_scaffold(tmp_path):
    root = create_parameter_template(tmp_path, "2026_User", ["H2O"])

    assert root == tmp_path / "2026_User"
    assert (root / "user_options.json").exists()
    assert json.loads((root / "user_options.json").read_text(encoding="utf-8")) == {}

    pure_path = root / "pure" / "water.csv"
    assert pure_path.exists()
    pure_lines = pure_path.read_text(encoding="utf-8").splitlines()
    assert pure_lines[0] == "component,m,s,e,e_assoc,vol_a,assoc_scheme,z,dielc,d_born,f_solv,MW"
    assert pure_lines[1].startswith("H2O,")

    for filename in ("k_ij.csv", "l_ij.csv", "k_hb_ij.csv"):
        matrix_lines = (root / "mixed" / "binary_interaction" / filename).read_text(encoding="utf-8").splitlines()
        assert matrix_lines[0] == "component,H2O"
        assert matrix_lines[1] == "H2O,0.0"

    rel_perm_path = root / "mixed" / "rel_perm" / "parameters.csv"
    assert rel_perm_path.read_text(encoding="utf-8").splitlines() == ["organic,a,b,c"]

    mixture = ePCSAFTMixture.from_dataset(root, ["H2O"], np.asarray([1.0]), 298.15)

    assert mixture.ncomp == 1
    assert mixture.parameters["z"].size == 0
    np.testing.assert_allclose(mixture.parameters["m"], [1.2047])
    np.testing.assert_allclose(mixture.parameters["k_ij"], np.zeros((1, 1)))


def test_runtime_options_reject_legacy_electrolyte_shorthand():
    with pytest.raises(KeyError, match="Unknown user_options key"):
        _resolve_runtime_options({"dielc_rule": "constant"})

    with pytest.raises(KeyError, match="Unknown user_options key"):
        _resolve_runtime_options({"debug": True})

    with pytest.raises(TypeError, match="elec_model\\['DH_model'\\] must be a dict"):
        _resolve_runtime_options({"elec_model": {"DH_model": 2}})

    with pytest.raises(TypeError, match="elec_model\\['born_model'\\] must be a dict"):
        _resolve_runtime_options({"elec_model": {"born_model": 1}})

    with pytest.raises(KeyError, match="unsupported key"):
        _resolve_runtime_options({"elec_model": {"born_rel_perm": "solvent"}})


def test_runtime_options_accept_autodiff_modes_and_preserve_explicit_overrides():
    user_options = {
        "elec_model": {
            "rel_perm": {"differential_mode": "autodiff"},
            "hc_model": {"dadx_differential_mode": "autodiff"},
            "disp_model": {"dadx_differential_mode": "autodiff"},
            "assoc_model": {"dadx_differential_mode": "autodiff"},
            "DH_model": {"mu_DH_model": {"differential_mode": "autodiff"}},
            "born_model": {"mu_born_model": {"differential_mode": "autodiff"}},
        }
    }

    resolved = _resolve_runtime_options(user_options)
    model = resolved["model"]
    runtime = resolved["runtime"]

    assert model["rel_perm"]["differential_mode"] == 2
    assert model["hc_model"]["dadx_differential_mode"] == 2
    assert model["disp_model"]["dadx_differential_mode"] == 2
    assert model["assoc_model"]["dadx_differential_mode"] == 2
    assert model["DH_model"]["mu_DH_model"]["differential_mode"] == 2
    assert model["born_model"]["mu_born_model"]["differential_mode"] == 2

    assert runtime["dielc_diff_mode"] == 2
    assert runtime["hc_dadx_diff_mode"] == 2
    assert runtime["disp_dadx_diff_mode"] == 2
    assert runtime["assoc_dadx_diff_mode"] == 2
    assert runtime["mu_DH_diff_mode"] == 2
    assert runtime["mu_born_diff_mode"] == 2

    minimized = minimize_user_options(user_options)
    assert minimized == user_options


def test_runtime_options_accept_salt_free_solvent_massfraction_dielectric_rule():
    resolved = _resolve_runtime_options({"elec_model": {"rel_perm": {"rule": "salt-free-massfraction"}}})

    assert resolved["model"]["rel_perm"]["rule"] == 9
    assert resolved["runtime"]["dielc_rule"] == 9


def test_salt_free_solvent_massfraction_dielectric_rule_reaches_native_runtime():
    mixture = ePCSAFTMixture.from_dataset(
        "2014_Held",
        ["H2O", "Butanol"],
        np.asarray([0.8, 0.2]),
        298.15,
        user_options={"elec_model": {"rel_perm": {"rule": "salt-free-massfraction"}}},
    )
    state = mixture.state(T=298.15, rho=1000.0, x=np.asarray([0.8, 0.2]))

    assert mixture.parameters["elec_model"]["rel_perm"]["rule"] == 9
    assert np.isfinite(state.compressibility_factor())


def test_runtime_options_default_to_auto_derivative_policy():
    resolved = _resolve_runtime_options({})
    model = resolved["model"]
    runtime = resolved["runtime"]

    assert model["rel_perm"]["differential_mode"] == 3
    assert model["hc_model"]["dadx_differential_mode"] == 3
    assert model["disp_model"]["dadx_differential_mode"] == 3
    assert model["assoc_model"]["dadx_differential_mode"] == 3
    assert model["DH_model"]["mu_DH_model"]["differential_mode"] == 3
    assert model["born_model"]["mu_born_model"]["differential_mode"] == 3

    assert runtime["dielc_diff_mode"] == 3
    assert runtime["hc_dadx_diff_mode"] == 3
    assert runtime["disp_dadx_diff_mode"] == 3
    assert runtime["assoc_dadx_diff_mode"] == 3
    assert runtime["mu_DH_diff_mode"] == 3
    assert runtime["mu_born_diff_mode"] == 3


def test_runtime_options_reject_removed_polar_model():
    with pytest.raises(KeyError, match="unsupported key"):
        _resolve_runtime_options({"elec_model": {"polar_model": {"dadx_differential_mode": "autodiff"}}})
