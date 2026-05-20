"""Tests for user-owned parameter template generation."""

import json

import numpy as np
import pytest

from epcsaft import ParameterSet, create_parameter_template, ePCSAFTMixture
from epcsaft.parameters import _resolve_runtime_options, minimize_user_options
from tests.helpers.numeric import assert_allclose


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
    assert_allclose(mixture.parameters["m"], [1.2047])
    assert_allclose(mixture.parameters["k_ij"], np.zeros((1, 1)))


def test_create_parameter_template_accepts_explicit_legacy_schema(tmp_path):
    root = create_parameter_template(tmp_path, "legacy_user", ["H2O", "Ethanol"], schema="legacy")

    assert (root / "pure" / "any_solvent.csv").exists()
    assert (root / "mixed" / "binary_interaction" / "k_ij.csv").exists()
    assert not (root / "parameter_set.json").exists()


def test_create_parameter_template_canonical_schema_creates_json_scaffold(tmp_path):
    root = create_parameter_template(tmp_path, "canonical_user", ["H2O", "Na+", "Cl-"], schema="canonical")

    assert root == tmp_path / "canonical_user"
    assert not (root / "pure").exists()
    assert json.loads((root / "user_options.json").read_text(encoding="utf-8")) == {}

    parameter_set = json.loads((root / "parameter_set.json").read_text(encoding="utf-8"))
    assert parameter_set["schema"] == "canonical"
    assert parameter_set["schema_version"] == 1
    assert parameter_set["components"] == ["H2O", "Na+", "Cl-"]
    assert parameter_set["binary_records"] == []

    pure_records = parameter_set["pure_records"]
    assert [record["component"] for record in pure_records] == ["H2O", "Na+", "Cl-"]
    assert {record["molar_mass_units"] for record in pure_records} == {"kg/mol"}
    assert all(record["molar_mass"] is None for record in pure_records)


def test_canonical_parameter_template_loads_through_dataset_seam(tmp_path):
    root = create_parameter_template(tmp_path, "canonical_runtime", ["H2O", "Ethanol"], schema="canonical")
    path = root / "parameter_set.json"
    payload = json.loads(path.read_text(encoding="utf-8"))
    payload["pure_records"] = [
        {
            "component": "H2O",
            "molar_mass": 18.01528e-3,
            "m": 1.2047,
            "sigma": 2.7927,
            "epsilon_k": 353.95,
            "charge": 0.0,
            "relative_permittivity": 78.3,
        },
        {
            "component": "Ethanol",
            "molar_mass": 46.07e-3,
            "m": 2.0,
            "sigma": 3.5,
            "epsilon_k": 250.0,
            "charge": 0.0,
        },
    ]
    payload["binary_records"] = [
        {
            "components": ["H2O", "Ethanol"],
            "k_ij": 0.021,
        }
    ]
    payload["runtime_options"] = {"source_tag": "canonical-template-test"}
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    (root / "user_options.json").write_text(
        json.dumps(
            {
                "canonical_user_options": {
                    "elec_model": {
                        "preset": "legacy-shorthand",
                        "base": "legacy-base",
                        "rel_perm": {"rule": "constant"},
                    }
                }
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )

    params = ParameterSet.from_dataset(root, ["H2O", "Ethanol"], np.asarray([0.5, 0.5]), 298.15)
    mixture = ePCSAFTMixture.from_dataset(root, ["H2O", "Ethanol"], np.asarray([0.5, 0.5]), 298.15)

    assert params.components == ("H2O", "Ethanol")
    assert params.to_runtime_dict()["k_ij"][0, 1] == pytest.approx(0.021)
    assert params.to_runtime_dict()["source_tag"] == "canonical-template-test"
    assert "preset" not in params.runtime_options["elec_model"]
    assert "base" not in params.runtime_options["elec_model"]
    assert params.runtime_options["elec_model"]["rel_perm"]["rule"] == "constant"
    assert mixture.species == ["H2O", "Ethanol"]
    assert mixture.parameters["k_ij"][0, 1] == pytest.approx(0.021)


def test_create_parameter_template_rejects_unknown_schema(tmp_path):
    with pytest.raises(ValueError, match="Unsupported parameter template schema"):
        create_parameter_template(tmp_path, "bad_schema", ["H2O"], schema="spreadsheet")


def test_runtime_options_reject_removed_electrolyte_shorthand():
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


def test_runtime_options_accept_cppad_modes_and_preserve_explicit_overrides():
    user_options = {
        "elec_model": {
            "rel_perm": {"differential_mode": "cppad"},
            "hc_model": {"dadx_differential_mode": "cppad"},
            "disp_model": {"dadx_differential_mode": "cppad"},
            "assoc_model": {"dadx_differential_mode": "cppad"},
            "DH_model": {"mu_DH_model": {"differential_mode": "cppad"}},
            "born_model": {"mu_born_model": {"differential_mode": "cppad"}},
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


def test_runtime_options_reject_removed_generic_derivative_mode():
    with pytest.raises(ValueError, match="Unknown rule option"):
        _resolve_runtime_options({"elec_model": {"rel_perm": {"differential_mode": "autodiff"}}})


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


def test_salt_free_solvent_massfraction_dielectric_rule_supports_ionic_pressure_state():
    mixture = ePCSAFTMixture.from_dataset(
        "2014_Held",
        ["H2O", "Butanol", "NH4+", "Cl-"],
        np.asarray([0.715, 0.27, 0.0075, 0.0075]),
        298.15,
        user_options={
            "elec_model": {
                "rel_perm": {"rule": "salt-free-massfraction"},
                "include_born_model": False,
            }
        },
    )
    state = mixture.state(T=298.15, P=101325.0, x=np.asarray([0.715, 0.27, 0.0075, 0.0075]))

    assert mixture.parameters["elec_model"]["rel_perm"]["rule"] == 9
    assert state.density() > 0.0
    assert np.isfinite(state.pressure())


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
        _resolve_runtime_options({"elec_model": {"polar_model": {"dadx_differential_mode": "cppad"}}})
