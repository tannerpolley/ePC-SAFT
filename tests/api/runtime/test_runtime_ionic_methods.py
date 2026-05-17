"""Object-oriented integration tests for the pybind11 native ePC-SAFT runtime."""

import numpy as np
import pytest

from tests.helpers.runtime_cases import (
    _assert_array,
    _ionic_state,
    _ionic_state_with_elec_model,
)


def test_ionic_activity_and_solution_methods_return_expected_values():
    state, species = _ionic_state()

    relative_permittivity = state.relative_permittivity()
    osmotic_coefficient = state.osmotic_coefficient()
    component_activity = state.activity_coefficient(species=species)
    mean_ionic_mole = state.activity_coefficient(species=species, mean_ionic_form=True, basis="mole")
    mean_ionic_molality = state.activity_coefficient(species=species, mean_ionic_form=True, basis="molality")
    solvation_free_energy = state.solvation_free_energy(species=species)
    diagnostics = state.state_diagnostics(species=species)

    assert relative_permittivity[0] == pytest.approx(78.075982)
    _assert_array(relative_permittivity[1], [78.09, 8.0, 8.0])
    _assert_array(osmotic_coefficient, [0.9739566103279091], rtol=1e-4)
    assert component_activity == {
        "water": pytest.approx(1.0000051724037697),
        "Na+": pytest.approx(0.9222113778654043),
        "Cl-": pytest.approx(0.9222258090371313),
    }
    assert mean_ionic_mole == {"Na+Cl-": pytest.approx(0.9222185934230398)}
    assert mean_ionic_molality == {"Na+Cl-": pytest.approx(0.9220341497043553)}
    assert solvation_free_energy == {
        "Na+": pytest.approx(-475461.4260703414),
        "Cl-": pytest.approx(-489572.50284416083),
    }

    _assert_array(diagnostics["relative_permittivity"][1], [78.09, 8.0, 8.0])
    _assert_array(diagnostics["osmotic_coefficient"], osmotic_coefficient)
    assert diagnostics["activity_coefficient"] == component_activity
    assert diagnostics["mean_ionic_activity_coefficient_mole"] == mean_ionic_mole
    assert diagnostics["mean_ionic_activity_coefficient_molality"] == mean_ionic_molality
    assert diagnostics["solvation_free_energy"] == solvation_free_energy
    _assert_array(
        diagnostics["fugacity_coefficient"], [0.031479320480733174, 4.651483659012546e-84, 1.5683276992772872e-86]
    )
    _assert_array(
        diagnostics["residual_chemical_potential"], [-10.682420304620588, -199.10395742942775, -204.79630395556683]
    )
    _assert_array(
        np.exp(state.fugacity_coefficient()),
        [0.031479320480733174, 4.651483659012546e-84, 1.5683276992772872e-86],
    )
    _assert_array(state.fugacity_coefficient(), [-3.458424439279275, -191.8799615776576, -197.57230810636238])
    _assert_array(
        state.fugacity_coefficient(natural_log=False),
        [0.031479320480733174, 4.651483659012546e-84, 1.5683276992772872e-86],
    )
    assert state.compressibility_factor() == pytest.approx(0.000728884077611683)
    assert state.residual_helmholtz() == pytest.approx(-9.7214027218058)
    assert state.temperature_derivative_residual_helmholtz() == pytest.approx(0.032388021640507005)
    assert state.residual_enthalpy() == pytest.approx(-26415.160790583413)
    assert state.residual_entropy() == pytest.approx(-59.523895812302186)
    assert state.residual_gibbs() == pytest.approx(-8668.111254145517)
    _assert_array(state.residual_chemical_potential(), [-10.682420304620588, -199.10395742942775, -204.79630395556683])
    assert state.pressure() == pytest.approx(100000.0)
    assert state.density() == pytest.approx(55344.274540081075)
    assert state.density(units="molar") == pytest.approx(55344.274540081075)
    assert state.molar_density() == pytest.approx(55344.274540081075)
    assert state.density(units="mass") == pytest.approx(997.1665703121223)
    assert state.density(units="kg/m^3") == pytest.approx(997.1665703121223)
    assert state.mass_density() == pytest.approx(997.1665703121223)
    _assert_array(
        diagnostics["fugacity_coefficient_terms"]["mu_total"],
        [-10.682420304620588, -199.10395742942775, -204.79630395556683],
    )
    _assert_array(
        diagnostics["fugacity_coefficient_terms"]["lnfugcoef_total"],
        [-3.4584244392944334, -191.87996157767273, -197.5723081063775],
    )

def test_rel_perm_cppad_matches_analytic_density_derivative_usage():
    base_state, _ = _ionic_state()
    ad_state, _ = _ionic_state_with_elec_model(
        {
            "rel_perm": {"differential_mode": "cppad"},
        }
    )

    base_dadx = base_state.dadx()
    ad_dadx = ad_state.dadx()
    _assert_array(ad_dadx["terms"]["ion"], base_dadx["terms"]["ion"], rtol=1e-7, atol=1e-9)
    _assert_array(ad_dadx["terms"]["born"], base_dadx["terms"]["born"], rtol=1e-7, atol=1e-9)

def test_default_dadx_reports_auto_derivative_policy():
    state, _ = _ionic_state()
    dadx = state.dadx()

    assert dadx["derivative_backend"]["hc"] == "analytic"
    assert dadx["derivative_backend"]["disp"] == "analytic"
    assert dadx["derivative_backend"]["ion"] == "analytic"
    assert dadx["derivative_backend"]["born"] == "analytic"
    assert dadx["derivative_backend"]["assoc"] == "analytic_implicit"
    assert dadx["derivative_available"] is True

def test_removed_derivative_backend_names_are_rejected():
    removed_backend = "f" + "d"
    with pytest.raises(ValueError, match="Unknown option value"):
        _ionic_state_with_elec_model(
            {
                "hc_model": {"dadx_differential_mode": removed_backend},
                "disp_model": {"dadx_differential_mode": removed_backend},
                "DH_model": {"mu_DH_model": {"differential_mode": removed_backend}},
                "born_model": {"mu_born_model": {"differential_mode": removed_backend}},
            }
        )

def test_hc_dadx_cppad_matches_analytic_terms():
    base_state, _ = _ionic_state()
    ad_state, _ = _ionic_state_with_elec_model(
        {
            "hc_model": {"dadx_differential_mode": "cppad"},
        }
    )

    _assert_array(ad_state.dadx()["terms"]["hc"], base_state.dadx()["terms"]["hc"], rtol=1e-7, atol=1e-9)

def test_mu_dh_cppad_matches_analytic_ion_terms():
    base_state, _ = _ionic_state()
    ad_state, _ = _ionic_state_with_elec_model(
        {
            "DH_model": {
                "mu_DH_model": {"differential_mode": "cppad"},
            },
        }
    )

    _assert_array(ad_state.dadx()["terms"]["ion"], base_state.dadx()["terms"]["ion"], rtol=1e-7, atol=1e-9)

def test_mu_born_cppad_matches_analytic_born_terms():
    base_state, _ = _ionic_state()
    ad_state, _ = _ionic_state_with_elec_model(
        {
            "born_model": {
                "mu_born_model": {"differential_mode": "cppad"},
            },
        }
    )

    _assert_array(ad_state.dadx()["terms"]["born"], base_state.dadx()["terms"]["born"], rtol=1e-7, atol=1e-9)
