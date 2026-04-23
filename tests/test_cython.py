# -*- coding: utf-8 -*-
"""Object-oriented integration tests for the native ePC-SAFT runtime."""

import numpy as np
import pytest

import epcsaft
from epcsaft import ePCSAFTMixture


def _neutral_state():
    species = ["A", "B", "C"]
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
    mix = ePCSAFTMixture.from_params(params, species=species)
    state = mix.state(T=233.15, x=np.array([0.1, 0.3, 0.6]), rho=14330.417110)
    return state, species


def _ionic_state():
    t = 298.15
    species = ["water", "Na+", "Cl-"]
    params = _ionic_params()
    mix = ePCSAFTMixture.from_params(params, species=species)
    state = mix.state(T=t, x=np.array([0.9998, 1.0e-4, 1.0e-4]), P=1.0e5)
    return state, species


def _ionic_state_with_elec_model(elec_model):
    t = 298.15
    species = ["water", "Na+", "Cl-"]
    params = _ionic_params()
    params["elec_model"] = elec_model
    mix = ePCSAFTMixture.from_params(params, species=species)
    state = mix.state(T=t, x=np.array([0.9998, 1.0e-4, 1.0e-4]), P=1.0e5)
    return state, species


def _ionic_params():
    t = 298.15
    s_water = 2.7927 + 10.11 * np.exp(-0.01775 * t) - 1.417 * np.exp(-0.01146 * t)
    return {
        "MW": np.asarray([18.01528e-3, 22.98e-3, 35.45e-3]),
        "m": np.asarray([1.2047, 1.0, 1.0]),
        "s": np.asarray([s_water, 2.8232, 2.7560]),
        "e": np.asarray([353.95, 230.0, 170.0]),
        "e_assoc": np.asarray([2425.7, 0.0, 0.0]),
        "vol_a": np.asarray([0.04509, 0.0, 0.0]),
        "assoc_scheme": ["2B", None, None],
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
    }


def _assert_array(actual, expected, rtol=1e-8, atol=1e-10):
    np.testing.assert_allclose(np.asarray(actual, dtype=float), np.asarray(expected, dtype=float), rtol=rtol, atol=atol)


def _sum_term_arrays(terms):
    total = None
    for value in terms.values():
        arr = np.asarray(value, dtype=float)
        total = arr.copy() if total is None else total + arr
    return total


def test_package_exports_are_available():
    assert hasattr(epcsaft, "ePCSAFTMixture")
    assert hasattr(epcsaft, "ePCSAFTState")
    assert hasattr(epcsaft, "ActivityCoefficientResult")


def test_from_params_rejects_legacy_electrolyte_keys():
    species = ["water", "Na+", "Cl-"]
    params = _ionic_params()
    params["dielc_rule"] = 1
    with pytest.raises(ValueError, match='Flat electrostatic params are no longer supported'):
        ePCSAFTMixture.from_params(params, species=species)

    params = _ionic_params()
    params["elec_model"] = {"born_rel_perm": "solvent"}
    with pytest.raises(ValueError, match='unsupported key'):
        ePCSAFTMixture.from_params(params, species=species)


def test_neutral_scalar_methods_return_expected_values():
    state, species = _neutral_state()

    assert state.density() == pytest.approx(14330.417110)
    assert state.density(units="molar") == pytest.approx(14330.417110)
    assert state.density(units="mol/m^3") == pytest.approx(14330.417110)
    assert state.molar_density() == pytest.approx(14330.417110)
    with pytest.raises(epcsaft.InputError, match="density units must be"):
        state.density(units="weird")
    with pytest.raises(epcsaft.InputError, match="Mass density requires component molecular weights"):
        state.mass_density()
    with pytest.raises(epcsaft.InputError, match="Mass density requires component molecular weights"):
        state.density(units="mass")
    assert state.pressure() == pytest.approx(1276374.1152948933)
    assert state.compressibility_factor() == pytest.approx(0.04594621208078564)
    assert state.residual_helmholtz() == pytest.approx(-3.54988545131505)
    assert state.temperature_derivative_residual_helmholtz() == pytest.approx(0.03077401856781036)
    assert state.residual_enthalpy() == pytest.approx(-15758.229958475444)
    assert state.residual_entropy() == pytest.approx(-55.751451436621096)
    assert state.residual_gibbs() == pytest.approx(-2759.779056027235)
    _assert_array(state.residual_chemical_potential(), [-1.1478687523834008, -3.6543804288405415, -5.488063725572939])

    fugacity_coefficient = state.fugacity_coefficient()
    fugacity_coefficient_coeff = state.fugacity_coefficient(natural_log=False)
    _assert_array(fugacity_coefficient, [1.9324151168689134, -0.5740965595882255, -2.407779856320623])
    _assert_array(fugacity_coefficient_coeff, [6.906169322700795, 0.5632134688356544, 0.09001491894620331])
    np.testing.assert_allclose(np.exp(fugacity_coefficient), fugacity_coefficient_coeff)


def test_state_method_aliases_match_canonical_methods():
    state, species = _ionic_state()
    aliases = state.method_aliases()
    assert aliases == {
        "pressure": "p",
        "density": "rho",
        "molar_density": "rho_molar",
        "mass_density": "rho_mass",
        "compressibility_factor": "z",
        "residual_helmholtz": "ares",
        "temperature_derivative_residual_helmholtz": "dadt",
        "composition_derivative_residual_helmholtz": "dadx",
        "residual_enthalpy": "hres",
        "residual_entropy": "sres",
        "residual_gibbs": "gres",
        "residual_chemical_potential": "mures",
        "activity_coefficient": "gamma",
        "fugacity_coefficient": "fugcoef",
        "relative_permittivity": "epsr",
        "osmotic_coefficient": "osmotic_coef",
        "state_diagnostics": "diag",
        "solvation_free_energy": "gsolv",
    }

    assert state.p() == pytest.approx(state.pressure())
    assert state.rho() == pytest.approx(state.density())
    assert state.rho_molar() == pytest.approx(state.molar_density())
    assert state.rho_mass() == pytest.approx(state.mass_density())
    assert state.z() == pytest.approx(state.compressibility_factor())
    assert state.ares() == pytest.approx(state.residual_helmholtz())
    assert state.dadt() == pytest.approx(state.temperature_derivative_residual_helmholtz())
    assert state.hres() == pytest.approx(state.residual_enthalpy())
    assert state.sres() == pytest.approx(state.residual_entropy())
    assert state.gres() == pytest.approx(state.residual_gibbs())
    assert state.dadx()["z_total"] == pytest.approx(state.composition_derivative_residual_helmholtz()["z_total"])
    _assert_array(state.mures(), state.residual_chemical_potential())
    _assert_array(state.fugcoef(), state.fugacity_coefficient())
    _assert_array(state.fugcoef(natural_log=False), state.fugacity_coefficient(natural_log=False))
    assert state.epsr()[0] == pytest.approx(state.relative_permittivity()[0])
    _assert_array(state.osmotic_coef(), state.osmotic_coefficient())
    assert state.gamma(species=species) == state.activity_coefficient(species=species)
    assert state.gamma(species=species, mean_ionic_form=True, basis="molality") == state.activity_coefficient(
        species=species, mean_ionic_form=True, basis="molality"
    )
    assert state.diag(species=species)["pressure"] == pytest.approx(state.state_diagnostics(species=species)["pressure"])
    assert state.gsolv(species=species) == state.solvation_free_energy(species=species)
    with pytest.raises(AttributeError):
        getattr(state, "fugacity_coefficient_terms")
    with pytest.raises(AttributeError):
        getattr(state, "lnfug_terms")


def test_state_contribution_term_payloads_match_totals():
    state, species = _ionic_state()

    ares = state.ares(return_contribution_terms=True)
    assert set(ares) == {"total", "terms"}
    assert set(ares["terms"]) == {"hc", "disp", "assoc", "ion", "born"}
    assert ares["total"] == pytest.approx(state.residual_helmholtz())
    assert sum(ares["terms"].values()) == pytest.approx(ares["total"])

    z_terms = state.z(return_contribution_terms=True)
    assert set(z_terms) == {"total", "terms"}
    assert set(z_terms["terms"]) == {"hc", "disp", "assoc", "ion", "born", "ideal"}
    assert z_terms["total"] == pytest.approx(state.compressibility_factor())
    assert sum(z_terms["terms"].values()) == pytest.approx(z_terms["total"])

    dadt = state.dadt(return_contribution_terms=True)
    assert set(dadt) == {"total", "terms"}
    assert set(dadt["terms"]) == {"hc", "disp", "assoc", "ion", "born"}
    assert dadt["total"] == pytest.approx(state.temperature_derivative_residual_helmholtz())
    assert sum(dadt["terms"].values()) == pytest.approx(dadt["total"])

    mures = state.mures(return_contribution_terms=True)
    assert set(mures) == {"total", "terms"}
    assert set(mures["terms"]) == {"hc", "disp", "assoc", "ion", "born"}
    _assert_array(mures["total"], state.residual_chemical_potential())
    _assert_array(_sum_term_arrays(mures["terms"]), mures["total"])

    dadx = state.dadx()
    assert set(dadx) == {"total", "terms", "ares_terms", "sum_x_terms", "z_raw_terms", "z_terms", "z_total"}
    assert set(dadx["terms"]) == {"hc", "disp", "assoc", "ion", "born"}
    assert set(dadx["ares_terms"]) == {"hc", "disp", "assoc", "ion", "born"}
    assert set(dadx["sum_x_terms"]) == {"hc", "disp", "assoc", "ion", "born"}
    assert set(dadx["z_raw_terms"]) == {"hc", "disp", "assoc", "ion", "born"}
    assert set(dadx["z_terms"]) == {"hc", "disp", "assoc", "ion", "born"}
    _assert_array(dadx["total"], _sum_term_arrays(dadx["terms"]))
    assert dadx["z_total"] == pytest.approx(state.compressibility_factor())
    assert dadx["ares_terms"]["hc"] == pytest.approx(state.ares(return_contribution_terms=True)["terms"]["hc"])

    fugcoef = state.fugcoef(return_contribution_terms=True)
    fugcoef_coeff = state.fugcoef(natural_log=False, return_contribution_terms=True)
    assert set(fugcoef) == {"total", "terms", "term_basis", "terms_total_natural_log"}
    assert set(fugcoef["terms"]) == {"hc", "disp", "assoc", "ion", "born"}
    _assert_array(fugcoef["total"], state.fugacity_coefficient())
    _assert_array(fugcoef_coeff["total"], state.fugacity_coefficient(natural_log=False))
    _assert_array(_sum_term_arrays(fugcoef["terms"]), fugcoef["terms_total_natural_log"])
    _assert_array(_sum_term_arrays(fugcoef_coeff["terms"]), fugcoef_coeff["terms_total_natural_log"])
    assert fugcoef["term_basis"] == "natural_log"


def test_dadrho_hierarchy_identities_hold_for_neutral_and_ionic_states():
    for factory in (_neutral_state, _ionic_state):
        state, species = factory()
        z_payload = state.z(return_contribution_terms=True)
        dadx = state.dadx()
        mures = state.mures()
        diagnostics = state.state_diagnostics(species=species)
        lnfug = np.asarray(diagnostics["fugacity_coefficient_terms"]["lnfugcoef_total"], dtype=float)

        z_residual = sum(float(dadx["z_terms"][key]) for key in ("hc", "disp", "assoc", "ion", "born"))
        assert z_residual == pytest.approx(state.z() - 1.0)
        assert z_payload["terms"]["ideal"] + z_residual == pytest.approx(z_payload["total"])

        reconstructed_mu = sum(
            dadx["ares_terms"][key] + dadx["z_raw_terms"][key] + dadx["terms"][key] - dadx["sum_x_terms"][key]
            for key in ("hc", "disp", "assoc", "ion", "born")
        )
        np.testing.assert_allclose(reconstructed_mu, mures)
        np.testing.assert_allclose(lnfug, mures - np.log(state.z()))


def test_neutral_composition_and_fugacity_terms_return_expected_values():
    state, species = _neutral_state()

    dadx = state.composition_derivative_residual_helmholtz()
    terms = state.state_diagnostics(species=species)["fugacity_coefficient_terms"]
    assert set(dadx) == {"total", "terms", "ares_terms", "sum_x_terms", "z_raw_terms", "z_terms", "z_total"}
    for key in ("hc", "disp", "assoc", "ion", "born"):
        np.testing.assert_allclose(dadx["terms"][key], terms[f"dadx_{key}"])
        np.testing.assert_allclose(dadx["ares_terms"][key], terms[f"a_{key}"])
        np.testing.assert_allclose(dadx["sum_x_terms"][key], terms[f"sum_x_dadx_{key}"])
        np.testing.assert_allclose(dadx["z_raw_terms"][key], terms[f"z_raw_{key}"])
        np.testing.assert_allclose(dadx["z_terms"][key], terms[f"z_{key}"])
    _assert_array(dadx["terms"]["hc"], [6.883394758977889, 9.511092421642308, 12.258394188218157])
    _assert_array(dadx["terms"]["disp"], [-7.506335725134188, -12.640545058126808, -17.221530131978034])
    _assert_array(dadx["terms"]["assoc"], [0.0, 0.0, 0.0])
    _assert_array(dadx["terms"]["ion"], [0.0, 0.0, 0.0])
    _assert_array(dadx["terms"]["born"], [0.0, 0.0, 0.0])
    assert dadx["z_total"] == pytest.approx(0.04594621208078564)
    _assert_array(terms["mu_total"], [-1.1478687523834008, -3.6543804288405415, -5.488063725572939])
    _assert_array(terms["lnfugcoef_total"], [1.9324151168689134, -0.5740965595882255, -2.407779856320623])
    assert terms["z_total"] == pytest.approx(0.04594621208078564)
    reconstructed_mu = sum(
        dadx["ares_terms"][key] + dadx["z_raw_terms"][key] + dadx["terms"][key] - dadx["sum_x_terms"][key]
        for key in ("hc", "disp", "assoc", "ion", "born")
    )
    np.testing.assert_allclose(reconstructed_mu, state.residual_chemical_potential())


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
    _assert_array(osmotic_coefficient, [0.9739566103279091], rtol=1e-6)
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
    _assert_array(diagnostics["osmotic_coefficient"], [0.9739566103279091], rtol=1e-6)
    assert diagnostics["activity_coefficient"] == component_activity
    assert diagnostics["mean_ionic_activity_coefficient_mole"] == mean_ionic_mole
    assert diagnostics["mean_ionic_activity_coefficient_molality"] == mean_ionic_molality
    assert diagnostics["solvation_free_energy"] == solvation_free_energy
    _assert_array(diagnostics["fugacity_coefficient"], [0.031479320480733174, 4.651483659012546e-84, 1.5683276992772872e-86])
    _assert_array(diagnostics["residual_chemical_potential"], [-10.682420304620588, -199.10395742942775, -204.79630395556683])
    _assert_array(
        np.exp(state.fugacity_coefficient()),
        [0.031479320480733174, 4.651483659012546e-84, 1.5683276992772872e-86],
    )
    _assert_array(state.fugacity_coefficient(), [-3.458424439279275, -191.8799615776576, -197.57230810636238])
    _assert_array(state.fugacity_coefficient(natural_log=False), [0.031479320480733174, 4.651483659012546e-84, 1.5683276992772872e-86])
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
    _assert_array(diagnostics["fugacity_coefficient_terms"]["mu_total"], [-10.682420304620588, -199.10395742942775, -204.79630395556683])
    _assert_array(diagnostics["fugacity_coefficient_terms"]["lnfugcoef_total"], [-3.4584244392944334, -191.87996157767273, -197.5723081063775])


def test_rel_perm_autodiff_matches_analytic_density_derivative_usage():
    base_state, _ = _ionic_state()
    ad_state, _ = _ionic_state_with_elec_model({
        "rel_perm": {"differential_mode": "autodiff"},
    })

    base_dadx = base_state.dadx()
    ad_dadx = ad_state.dadx()
    _assert_array(ad_dadx["terms"]["ion"], base_dadx["terms"]["ion"], rtol=1e-7, atol=1e-9)
    _assert_array(ad_dadx["terms"]["born"], base_dadx["terms"]["born"], rtol=1e-7, atol=1e-9)


def test_hc_dadx_autodiff_matches_analytic_terms():
    base_state, _ = _ionic_state()
    ad_state, _ = _ionic_state_with_elec_model({
        "hc_model": {"dadx_differential_mode": "autodiff"},
    })

    _assert_array(ad_state.dadx()["terms"]["hc"], base_state.dadx()["terms"]["hc"], rtol=1e-7, atol=1e-9)


def test_mu_dh_autodiff_matches_analytic_ion_terms():
    base_state, _ = _ionic_state()
    ad_state, _ = _ionic_state_with_elec_model({
        "DH_model": {
            "mu_DH_model": {"differential_mode": "autodiff"},
        },
    })

    _assert_array(ad_state.dadx()["terms"]["ion"], base_state.dadx()["terms"]["ion"], rtol=1e-7, atol=1e-9)


def test_mu_born_autodiff_matches_analytic_born_terms():
    base_state, _ = _ionic_state()
    ad_state, _ = _ionic_state_with_elec_model({
        "born_model": {
            "mu_born_model": {"differential_mode": "autodiff"},
        },
    })

    _assert_array(ad_state.dadx()["terms"]["born"], base_state.dadx()["terms"]["born"], rtol=1e-7, atol=1e-9)


def test_state_diagnostics_matches_public_methods():
    neutral_state, neutral_species = _neutral_state()
    neutral_diag = neutral_state.state_diagnostics(species=neutral_species)
    assert neutral_diag["activity_coefficient"] == {}
    assert neutral_diag["mean_ionic_activity_coefficient_mole"] == {}
    assert neutral_diag["mean_ionic_activity_coefficient_molality"] == {}
    assert neutral_diag["solvation_free_energy"] == {}
    assert neutral_diag["relative_permittivity"] is None
    assert neutral_diag["osmotic_coefficient"] is None
    assert neutral_diag["mass_density"] is None
    assert neutral_diag["pressure"] == pytest.approx(neutral_state.pressure())
    assert neutral_diag["density"] == pytest.approx(neutral_state.density())
    assert neutral_diag["density_molar"] == pytest.approx(neutral_state.molar_density())
    assert neutral_diag["compressibility_factor"] == pytest.approx(neutral_state.compressibility_factor())
    _assert_array(neutral_diag["residual_chemical_potential"], neutral_state.residual_chemical_potential())
    _assert_array(neutral_diag["fugacity_coefficient"], neutral_state.fugacity_coefficient(natural_log=False))

    ionic_state, ionic_species = _ionic_state()
    ionic_diag = ionic_state.state_diagnostics(species=ionic_species)
    assert ionic_diag["activity_coefficient"] == ionic_state.activity_coefficient(species=ionic_species)
    assert ionic_diag["density"] == pytest.approx(ionic_state.molar_density())
    assert ionic_diag["density_molar"] == pytest.approx(ionic_state.molar_density())
    assert ionic_diag["mass_density"] == pytest.approx(ionic_state.mass_density())
    assert ionic_diag["mean_ionic_activity_coefficient_mole"] == ionic_state.activity_coefficient(
        species=ionic_species,
        mean_ionic_form=True,
        basis="mole",
    )
    assert ionic_diag["mean_ionic_activity_coefficient_molality"] == ionic_state.activity_coefficient(
        species=ionic_species,
        mean_ionic_form=True,
        basis="molality",
    )
    assert ionic_diag["solvation_free_energy"] == ionic_state.solvation_free_energy(species=ionic_species)
    _assert_array(ionic_diag["relative_permittivity"][1], ionic_state.relative_permittivity()[1])
    _assert_array(ionic_diag["osmotic_coefficient"], ionic_state.osmotic_coefficient())
    _assert_array(ionic_diag["residual_chemical_potential"], ionic_state.residual_chemical_potential())
    _assert_array(ionic_diag["fugacity_coefficient"], ionic_state.fugacity_coefficient(natural_log=False))


def test_state_constructor_rejects_invalid_pressure_and_density():
    mix = ePCSAFTMixture.from_params({
        "m": np.asarray([1.0]),
        "s": np.asarray([3.0]),
        "e": np.asarray([150.0]),
    })

    with pytest.raises(epcsaft.InputError):
        mix.state(T=300.0, x=np.asarray([1.0]), P=0.0)

    with pytest.raises(epcsaft.InputError):
        mix.state(T=300.0, x=np.asarray([1.0]), rho=-1.0)


def test_pressure_based_state_matches_equivalent_density_state():
    mix = ePCSAFTMixture.from_params({
        "m": np.asarray([2.8149]),
        "s": np.asarray([3.7169]),
        "e": np.asarray([285.69]),
    }, species=["Toluene"])

    state_tp = mix.state(T=320.0, x=np.asarray([1.0]), P=101325.0, phase="liq")
    state_trho = mix.state(T=320.0, x=np.asarray([1.0]), rho=state_tp.density(), phase="liq")

    assert state_tp.density() == pytest.approx(state_trho.density())
    assert state_tp.pressure() == pytest.approx(state_trho.pressure())
    assert state_tp.compressibility_factor() == pytest.approx(state_trho.compressibility_factor())
    assert state_tp.residual_helmholtz() == pytest.approx(state_trho.residual_helmholtz())


def test_pressure_based_state_raises_solver_error_during_construction():
    mix = ePCSAFTMixture.from_params({
        "m": np.asarray([1.0]),
        "s": np.asarray([3.0]),
        "e": np.asarray([150.0]),
    }, species=["A"])

    with pytest.raises(epcsaft.SolutionError):
        mix.state(T=300.0, x=np.asarray([1.0]), P=1.0e-12, phase="liq")
