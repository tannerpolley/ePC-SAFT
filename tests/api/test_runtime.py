"""Object-oriented integration tests for the pybind11 native ePC-SAFT runtime."""

import numpy as np
import pytest

import epcsaft
import epcsaft.epcsaft as epcsaft_module
import epcsaft.ipopt_backend as ipopt_backend
from epcsaft import ePCSAFTMixture
from tests.helpers.runtime_cases import (
    _assert_array,
    _ionic_params,
    _ionic_state,
    _ionic_state_with_elec_model,
    _neutral_state,
    _sum_term_arrays,
)


def test_package_exports_are_available():
    assert hasattr(epcsaft, "ePCSAFTMixture")
    assert hasattr(epcsaft, "ePCSAFTState")
    assert hasattr(epcsaft, "ActivityCoefficientResult")
    assert isinstance(epcsaft.__version__, str)
    assert epcsaft.__version__ != "0+unknown"
    assert isinstance(epcsaft.__git_commit__, str)
    assert callable(epcsaft.runtime_build_info)
    assert callable(epcsaft.capabilities)
    assert callable(epcsaft.evaluate_fugacity_coefficients)
    assert callable(epcsaft.evaluate_fugacity_coefficients_batch)
    assert callable(epcsaft.validate_dataset_bundle)


def test_organized_public_import_modules_are_available():
    import epcsaft.diagnostics
    import epcsaft.electrolyte
    import epcsaft.eos
    import epcsaft.reactive

    assert epcsaft.eos.ePCSAFTMixture is epcsaft.ePCSAFTMixture
    assert epcsaft.eos.Mixture is epcsaft.ePCSAFTMixture
    assert epcsaft.eos.State is epcsaft.ePCSAFTState
    assert epcsaft.electrolyte.ElectrolyteLLEProblem is epcsaft.ElectrolyteLLEProblem
    assert epcsaft.reactive.ReactiveSpeciationProblem is epcsaft.ReactiveSpeciationProblem
    assert epcsaft.diagnostics.capabilities is epcsaft.capabilities


def test_runtime_build_info_and_capabilities_are_json_like():
    info = epcsaft.runtime_build_info()
    assert info["package_version"] == epcsaft.__version__
    assert "source_git_commit" in info
    assert info["native_extension_available"] is True

    capabilities = epcsaft.capabilities()
    assert capabilities["native_extension"] is True
    ipopt = capabilities["optimizers"]["ipopt"]
    assert ipopt["backend"] == "cyipopt"
    fit_contract = capabilities["regression"]["reactive_electrolyte_batch_context"]["fit_status_contract"]
    assert fit_contract["available"] is True
    assert "bounded_incomplete" not in fit_contract["statuses"]
    assert fit_contract["public_placeholder_statuses"] == []
    mixed_regression = capabilities["regression"]["reactive_electrolyte_batch_context"][
        "bounded_mixed_pressure_speciation_regression"
    ]
    assert mixed_regression["available"] is True
    assert mixed_regression["status"] == "production"
    assert mixed_regression["supports_pressure_targets"] is True
    assert mixed_regression["supports_speciation_targets"] is True
    assert mixed_regression["supports_bounds"] is True
    assert mixed_regression["native_hot_loop"] is False
    assert ipopt["available"] is ipopt_backend.cyipopt_available()
    assert ipopt["formulations"] == ["bound_constrained_residual_minimization"]
    assert ipopt["full_constrained_nlp_available"] is False
    assert ipopt["default_auto_uses_ipopt"] is False
    assert ipopt["exact_hessian_available"] is False
    assert info["optional_dependencies"]["cyipopt"]["available"] is ipopt["available"]
    assert capabilities["equilibrium"]["neutral_tp_flash"]["available"] is True
    assert capabilities["equilibrium"]["neutral_bubble_dew"] == {
        "available": True,
        "backend": "native_state_fugacity_with_python_scalar_root",
        "methods": ["bubble_p", "bubble_t", "dew_p", "dew_t"],
        "status": "production",
    }
    electrolyte_bubble = capabilities["equilibrium"]["electrolyte_bubble_pressure"]
    assert electrolyte_bubble["available"] is True
    assert electrolyte_bubble["backend"] == "native"
    assert electrolyte_bubble["scope"] == "fixed liquid composition with neutral vapor species; ions remain liquid-only"
    assert capabilities["equilibrium"]["electrolyte_lle"]["default_auto_uses_ipopt"] is False
    assert capabilities["equilibrium"]["electrolyte_lle"]["full_constrained_nlp_available"] is False
    reactive_bubble = capabilities["equilibrium"]["reactive_electrolyte_bubble"]
    assert reactive_bubble["available"] is True
    assert reactive_bubble["backend"] == "native"
    assert (
        reactive_bubble["scope"]
        == "native chemical speciation followed by native fixed-liquid electrolyte bubble pressure"
    )
    assert capabilities["equilibrium"]["reactive_speciation"]["default_auto_uses_ipopt"] is False
    assert capabilities["equilibrium"]["reactive_speciation"]["full_constrained_nlp_available"] is False
    assert (
        capabilities["equilibrium"]["reactive_speciation"]["jacobian_auto_policy"]
        == "analytic_where_available_else_backend_unavailable"
    )
    assert capabilities["equilibrium"]["reactive_speciation"]["derivative_gap_status"] == "backend_unavailable"
    assert capabilities["equilibrium"]["reactive_speciation"]["explicit_autodiff_raises_when_unavailable"] is True
    assert capabilities["regression"]["pure_neutral"]["backend"] == "native"
    reactive_regression = capabilities["regression"]["reactive_electrolyte_residuals"]
    assert reactive_regression["available"] is True
    assert reactive_regression["backend"] == "python_orchestrated_native_solvers"
    assert "downstream-owned" in reactive_regression["scope"]
    batch_context = capabilities["regression"]["reactive_electrolyte_batch_context"]
    assert batch_context["available"] is True
    assert batch_context["backend"] == "python_batched_native_solvers"
    assert "ReactiveElectrolyteRegressionContext" in batch_context["classes"]
    assert batch_context["methods"] == ["evaluate_objective"]
    assert capabilities["equilibrium"]["problem_objects"]["entrypoint"] == "mixture.solve_equilibrium(problem)"
    assert (
        capabilities["equilibrium"]["contribution_maps"]["activity_coefficient_term_decomposition_available"] is False
    )


def test_cyipopt_import_prepares_configured_windows_dll_directory(monkeypatch, tmp_path):
    calls: list[str] = []
    dll_dir = tmp_path / "ipopt-bin"
    dll_dir.mkdir()

    monkeypatch.setenv("EPCSAFT_IPOPT_DLL_DIR", str(dll_dir))
    monkeypatch.setenv("PATH", "base-path")
    monkeypatch.setattr(ipopt_backend.os, "add_dll_directory", lambda path: calls.append(str(path)), raising=False)
    ipopt_backend._prepare_ipopt_dll_search_path()

    assert calls == [str(dll_dir)]


def test_fast_fugacity_helper_matches_state_call_and_reports_density() -> None:
    state, _ = _ionic_state()
    mix = state.mixture
    helper = epcsaft.evaluate_fugacity_coefficients(
        mix,
        T=state.T,
        x=state.x,
        P=state.pressure(),
        phase="liq",
        natural_log=True,
    )

    _assert_array(helper["ln_fugacity_coefficient"], state.fugacity_coefficient(natural_log=True))
    assert helper["density"] == pytest.approx(state.molar_density())
    assert helper["phase"] == "liq"


def test_batch_fugacity_helper_matches_scalar_rows() -> None:
    state, _ = _ionic_state()
    mix = state.mixture
    rows = [
        {"T": state.T, "P": state.pressure(), "x": state.x, "phase": "liq"},
        {"T": state.T, "rho": state.molar_density(), "x": state.x, "phase": "liq"},
    ]

    batch = epcsaft.evaluate_fugacity_coefficients_batch(mix, rows=rows, natural_log=True)

    assert len(batch) == 2
    for payload in batch:
        _assert_array(payload["ln_fugacity_coefficient"], state.fugacity_coefficient(natural_log=True))
        assert payload["density"] == pytest.approx(state.molar_density())


def test_state_accepts_rho_seed_alias() -> None:
    state, _ = _ionic_state()
    mix = state.mixture
    seeded = mix.state(T=state.T, x=state.x, P=state.pressure(), phase="liq", rho_seed=state.molar_density())

    assert seeded.molar_density() == pytest.approx(state.molar_density())


def test_validate_dataset_bundle_reports_reaction_and_charge_errors() -> None:
    report = epcsaft.validate_dataset_bundle(
        {
            "m": np.asarray([1.0, 1.0]),
            "s": np.asarray([3.0, 3.0]),
            "e": np.asarray([200.0, 200.0]),
            "z": np.asarray([1.0, 1.0]),
            "d_born": np.asarray([3.0, np.nan]),
        },
        species=["Na+", "Cl-"],
        reactions=[epcsaft.ReactionDefinition({"Na+": -1.0, "Missing": 1.0}, 0.0)],
    )

    assert report["valid"] is False
    assert any("charge sign mismatch" in error for error in report["errors"])
    assert any("non-finite d_born" in error for error in report["errors"])
    assert any("Unknown species 'Missing'" in error for error in report["errors"])


def test_from_params_rejects_legacy_electrolyte_keys():
    species = ["water", "Na+", "Cl-"]
    params = _ionic_params()
    params["dielc_rule"] = 1
    with pytest.raises(ValueError, match="Flat electrostatic params are no longer supported"):
        ePCSAFTMixture.from_params(params, species=species)

    params = _ionic_params()
    params["elec_model"] = {"born_rel_perm": "solvent"}
    with pytest.raises(ValueError, match="unsupported key"):
        ePCSAFTMixture.from_params(params, species=species)


def test_neutral_scalar_methods_return_expected_values():
    state, _species = _neutral_state()

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
    assert state.diag(species=species)["pressure"] == pytest.approx(
        state.state_diagnostics(species=species)["pressure"]
    )
    assert state.gsolv(species=species) == state.solvation_free_energy(species=species)
    with pytest.raises(AttributeError):
        _ = state.fugacity_coefficient_terms
    with pytest.raises(AttributeError):
        _ = state.lnfug_terms


def test_state_contribution_term_payloads_match_totals():
    state, _species = _ionic_state()

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
    assert set(dadx) == {
        "total",
        "terms",
        "ares_terms",
        "sum_x_terms",
        "z_raw_terms",
        "z_terms",
        "z_total",
        "derivative_backend",
        "derivative_available",
        "backend_unavailable_reason",
    }
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


def test_state_contribution_map_aliases_use_public_family_names():
    state, _species = _ionic_state()

    families = {"hard_chain", "dispersion", "association", "ionic", "born"}
    helmholtz = state.helmholtz_contributions()
    residual_helmholtz = state.residual_helmholtz_contributions()
    pressure = state.pressure_contributions()
    chemical_potential = state.chemical_potential_contributions()
    ln_phi = state.ln_fugacity_coefficient_contributions()

    assert set(helmholtz["terms"]) == families
    assert set(residual_helmholtz["terms"]) == families
    assert families | {"ideal"} == set(pressure["terms"])
    assert set(chemical_potential["terms"]) == families
    assert set(ln_phi["terms"]) == families
    assert helmholtz["term_basis"] == "dimensionless_residual_helmholtz"
    assert pressure["term_basis"] == "pressure_from_compressibility_factor"
    np.testing.assert_allclose(ln_phi["total"], state.fugacity_coefficient())


def test_activity_coefficient_contribution_map_is_explicitly_unsupported():
    state, _species = _ionic_state()

    with pytest.raises(NotImplementedError, match="Activity coefficients are available"):
        state.activity_coefficient_contributions()


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
    assert set(dadx) == {
        "total",
        "terms",
        "ares_terms",
        "sum_x_terms",
        "z_raw_terms",
        "z_terms",
        "z_total",
        "derivative_backend",
        "derivative_available",
        "backend_unavailable_reason",
    }
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


def test_rel_perm_autodiff_matches_analytic_density_derivative_usage():
    base_state, _ = _ionic_state()
    ad_state, _ = _ionic_state_with_elec_model(
        {
            "rel_perm": {"differential_mode": "autodiff"},
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
    assert dadx["derivative_backend"]["assoc"] == "analytic"
    assert dadx["derivative_available"] is True
    assert dadx["backend_unavailable_reason"] == ""


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


def test_hc_dadx_autodiff_matches_analytic_terms():
    base_state, _ = _ionic_state()
    ad_state, _ = _ionic_state_with_elec_model(
        {
            "hc_model": {"dadx_differential_mode": "autodiff"},
        }
    )

    _assert_array(ad_state.dadx()["terms"]["hc"], base_state.dadx()["terms"]["hc"], rtol=1e-7, atol=1e-9)


def test_mu_dh_autodiff_matches_analytic_ion_terms():
    base_state, _ = _ionic_state()
    ad_state, _ = _ionic_state_with_elec_model(
        {
            "DH_model": {
                "mu_DH_model": {"differential_mode": "autodiff"},
            },
        }
    )

    _assert_array(ad_state.dadx()["terms"]["ion"], base_state.dadx()["terms"]["ion"], rtol=1e-7, atol=1e-9)


def test_mu_born_autodiff_matches_analytic_born_terms():
    base_state, _ = _ionic_state()
    ad_state, _ = _ionic_state_with_elec_model(
        {
            "born_model": {
                "mu_born_model": {"differential_mode": "autodiff"},
            },
        }
    )

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
    mix = ePCSAFTMixture.from_params(
        {
            "m": np.asarray([1.0]),
            "s": np.asarray([3.0]),
            "e": np.asarray([150.0]),
        }
    )

    with pytest.raises(epcsaft.InputError):
        mix.state(T=300.0, x=np.asarray([1.0]), P=0.0)

    with pytest.raises(epcsaft.InputError):
        mix.state(T=300.0, x=np.asarray([1.0]), rho=-1.0)


def test_state_constructor_rejects_composition_length_mismatch_before_native_call():
    mix = ePCSAFTMixture.from_params(
        {
            "m": np.asarray([1.0]),
            "s": np.asarray([3.0]),
            "e": np.asarray([150.0]),
        }
    )

    with pytest.raises(epcsaft.InputError, match="composition length"):
        mix.state(T=300.0, x=np.asarray([0.5, 0.5]), rho=100.0)


def test_pressure_based_state_matches_equivalent_density_state():
    mix = ePCSAFTMixture.from_params(
        {
            "m": np.asarray([2.8149]),
            "s": np.asarray([3.7169]),
            "e": np.asarray([285.69]),
        },
        species=["Toluene"],
    )

    state_tp = mix.state(T=320.0, x=np.asarray([1.0]), P=101325.0, phase="liq")
    state_trho = mix.state(T=320.0, x=np.asarray([1.0]), rho=state_tp.density(), phase="liq")

    assert state_tp.density() == pytest.approx(state_trho.density())
    assert state_tp.pressure() == pytest.approx(state_trho.pressure())
    assert state_tp.compressibility_factor() == pytest.approx(state_trho.compressibility_factor())
    assert state_tp.residual_helmholtz() == pytest.approx(state_trho.residual_helmholtz())


def test_pressure_based_state_accepts_density_guess_without_changing_closure():
    mix = ePCSAFTMixture.from_params(
        {
            "m": np.asarray([2.8149]),
            "s": np.asarray([3.7169]),
            "e": np.asarray([285.69]),
        },
        species=["Toluene"],
    )

    base = mix.state(T=320.0, x=np.asarray([1.0]), P=101325.0, phase="liq")
    seeded = mix.state(T=320.0, x=np.asarray([1.0]), P=101325.0, phase="liq", rho_guess=base.density())

    assert seeded.density() == pytest.approx(base.density())
    assert seeded.pressure() == pytest.approx(base.pressure())
    assert seeded.fugacity_coefficient() == pytest.approx(base.fugacity_coefficient())


def test_pressure_based_state_with_poor_density_guess_falls_back_safely():
    mix = ePCSAFTMixture.from_params(
        {
            "m": np.asarray([2.8149]),
            "s": np.asarray([3.7169]),
            "e": np.asarray([285.69]),
        },
        species=["Toluene"],
    )

    reference = mix.state(T=320.0, x=np.asarray([1.0]), P=101325.0, phase="liq")
    seeded = mix.state(T=320.0, x=np.asarray([1.0]), P=101325.0, phase="liq", rho_guess=1.0e-6)

    assert seeded.density() == pytest.approx(reference.density())
    assert seeded.pressure() == pytest.approx(reference.pressure())


@pytest.mark.parametrize("rho_guess", [0.0, -1.0, np.inf, np.nan])
def test_pressure_based_state_rejects_invalid_density_guess(rho_guess):
    mix = ePCSAFTMixture.from_params(
        {
            "m": np.asarray([1.0]),
            "s": np.asarray([3.0]),
            "e": np.asarray([150.0]),
        }
    )

    with pytest.raises(epcsaft.InputError, match="rho_guess must be finite and positive"):
        mix.state(T=300.0, x=np.asarray([1.0]), P=101325.0, rho_guess=rho_guess)


def test_density_guess_is_only_valid_for_pressure_based_states():
    mix = ePCSAFTMixture.from_params(
        {
            "m": np.asarray([1.0]),
            "s": np.asarray([3.0]),
            "e": np.asarray([150.0]),
        }
    )

    with pytest.raises(epcsaft.InputError, match="rho_guess is only supported"):
        mix.state(T=300.0, x=np.asarray([1.0]), rho=100.0, rho_guess=100.0)


def test_check_density_reports_pressure_consistency_diagnostics():
    mix = ePCSAFTMixture.from_params(
        {
            "m": np.asarray([2.8149]),
            "s": np.asarray([3.7169]),
            "e": np.asarray([285.69]),
        },
        species=["Toluene"],
    )

    state = mix.state(T=320.0, x=np.asarray([1.0]), P=101325.0, phase="liq")
    matching = mix.check_density(T=320.0, x=np.asarray([1.0]), P=101325.0, rho=state.density(), phase="liq")
    perturbed = mix.check_density(T=320.0, x=np.asarray([1.0]), P=101325.0, rho=state.density() * 0.9, phase="liq")

    assert matching["within_tolerance"] is True
    assert matching["pressure_residual"] == pytest.approx(0.0, abs=1.0e-6)
    assert matching["state"].density() == pytest.approx(state.density())
    assert perturbed["within_tolerance"] is False
    assert abs(perturbed["pressure_residual"]) > abs(matching["pressure_residual"])
    assert set(perturbed) == {
        "density",
        "pressure_target",
        "pressure_from_density",
        "pressure_residual",
        "relative_pressure_residual",
        "within_tolerance",
        "state",
    }


def test_pressure_based_state_raises_solver_error_during_construction():
    mix = ePCSAFTMixture.from_params(
        {
            "m": np.asarray([1.0]),
            "s": np.asarray([3.0]),
            "e": np.asarray([150.0]),
        },
        species=["A"],
    )

    with pytest.raises(epcsaft.SolutionError) as excinfo:
        mix.state(T=300.0, x=np.asarray([1.0]), P=1.0e-12, phase="liq")

    message = str(excinfo.value)
    assert "pressure-based state solve failed" in message
    assert "T=300.0" in message
    assert "P=1e-12" in message
    assert "phase=liq" in message
    assert "x=[1.0]" in message
    assert excinfo.value.__cause__ is not None


def test_density_based_native_constructor_failure_raises_public_solution_error(monkeypatch):
    mix = ePCSAFTMixture.from_params(
        {
            "m": np.asarray([1.0]),
            "s": np.asarray([3.0]),
            "e": np.asarray([150.0]),
        },
        species=["A"],
    )

    original_native_state = epcsaft_module._core.NativeState

    def raising_native_state(*args, **kwargs):
        raise RuntimeError("simulated density native failure")

    monkeypatch.setattr(epcsaft_module._core, "NativeState", raising_native_state)
    with pytest.raises(epcsaft.SolutionError) as excinfo:
        mix.state(T=300.0, x=np.asarray([1.0]), rho=100.0, phase="liq")
    monkeypatch.setattr(epcsaft_module._core, "NativeState", original_native_state)

    message = str(excinfo.value)
    assert "density-based state solve failed" in message
    assert "T=300.0" in message
    assert "rho=100.0" in message
    assert "phase=liq" in message
    assert "ncomp=1" in message
    assert "x=[1.0]" in message
    assert "simulated density native failure" in message
    assert excinfo.value.__cause__ is not None
