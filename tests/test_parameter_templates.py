# -*- coding: utf-8 -*-
"""Tests for user-owned parameter template generation."""

import json

import numpy as np
import pytest

from epcsaft import create_parameter_template
from epcsaft import ePCSAFTMixture
from epcsaft.parameters import _resolve_runtime_options


def test_create_parameter_template_creates_loadable_scaffold(tmp_path):
    root = create_parameter_template(tmp_path, "2026_User", ["H2O"])

    assert root == tmp_path / "2026_User"
    assert (root / "user_options.json").exists()
    assert json.loads((root / "user_options.json").read_text(encoding="utf-8")) == {}

    pure_path = root / "pure" / "water.csv"
    assert pure_path.exists()
    pure_lines = pure_path.read_text(encoding="utf-8").splitlines()
    assert pure_lines[0] == "component,m,s,e,e_assoc,vol_a,assoc_scheme,dipm,dip_num,z,dielc,d_born,f_solv,MW"
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

    with pytest.raises(TypeError, match="elec_model\\['DH_model'\\] must be a dict"):
        _resolve_runtime_options({"elec_model": {"DH_model": 2}})

    with pytest.raises(TypeError, match="elec_model\\['born_model'\\] must be a dict"):
        _resolve_runtime_options({"elec_model": {"born_model": 1}})

    with pytest.raises(KeyError, match="unsupported key"):
        _resolve_runtime_options({"elec_model": {"born_rel_perm": "solvent"}})
