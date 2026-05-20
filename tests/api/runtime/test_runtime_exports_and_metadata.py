"""Object-oriented integration tests for the pybind11 native ePC-SAFT runtime."""

import inspect

import numpy as np
import pytest

import epcsaft
from epcsaft import _core, ePCSAFTMixture, runtime
from tests.helpers.runtime_cases import (
    _assert_array,
    _ionic_params,
    _ionic_state,
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


def test_mixture_equilibrium_facades_do_not_expose_phase_seed_keyword():
    removed_keyword = "initial" + "_phases"

    assert removed_keyword not in inspect.signature(epcsaft.ePCSAFTMixture.equilibrium).parameters
    assert removed_keyword not in inspect.signature(epcsaft.ePCSAFTMixture.equilibrium_curve).parameters
    assert removed_keyword not in inspect.signature(epcsaft.ePCSAFTMixture.lle_tp).parameters
    assert removed_keyword not in inspect.signature(epcsaft.ePCSAFTMixture.electrolyte_lle_tp).parameters


def test_ipopt_backend_info_missing_smoke_uses_probe_status(monkeypatch: pytest.MonkeyPatch):
    monkeypatch.delattr(_core, "_native_ipopt_smoke")

    info = runtime._native_ipopt_backend_info()

    assert info["backend"] == "ipopt"
    assert info["status"] == "ipopt_probe_missing"
    assert info["available"] is False
    assert info["adapter_kind"] == "native_tnlp_adapter"


def test_runtime_build_info_and_capabilities_are_json_like():
    info = epcsaft.runtime_build_info()
    assert info["package_version"] == epcsaft.__version__
    assert "source_git_commit" in info
    assert info["native_extension_available"] is True
    assert "native_dependencies" in info

    capabilities = epcsaft.capabilities()
    assert capabilities["native_extension"] is True
    ipopt = capabilities["optimizers"]["ipopt"]
    assert ipopt["backend"] == "ipopt"
    mixed_regression = capabilities["regression"]["reactive_electrolyte_batch_context"][
        "mixed_pressure_speciation_residual_context"
    ]
    assert mixed_regression["available"] is True
    assert "status" not in mixed_regression
    assert mixed_regression["production_optimizer"] is False
    assert mixed_regression["optimizer"] is None
    assert mixed_regression["supports_pressure_targets"] is True
    assert mixed_regression["supports_speciation_targets"] is True
    assert mixed_regression["validates_parameter_bounds"] is True
    assert mixed_regression["native_hot_loop"] is False
    assert mixed_regression["ceres"]["production"] is False
    assert ipopt["available"] is info["native_dependencies"]["ipopt"]["available"]
    assert ipopt["formulations"] == ["thermodynamic_constrained_nlp"]
    assert ipopt["adapter_available"] is info["native_dependencies"]["ipopt"]["adapter_available"]
    assert ipopt["production"] is ipopt["available"]
    assert ipopt["public_routes"] == [
        "reactive_speciation:ideal_mole_fraction",
        "reactive_speciation:mole_fraction_activity",
        "reactive_speciation:concentration",
        "neutral_tp_flash",
        "neutral_stability",
        "electrolyte_stability",
        "neutral_lle_flash",
        "neutral_bubble_p",
        "neutral_bubble_t",
        "neutral_dew_p",
        "neutral_dew_t",
        "electrolyte_lle",
        "electrolyte_bubble_pressure",
    ]
    assert info["native_dependencies"]["ipopt"]["available"] is ipopt["available"]
    cppad = info["native_dependencies"]["cppad"]
    assert cppad["backend"] == "cppad"
    assert cppad["status"] == "enabled_available"
    assert cppad["compiled"] is True
    assert cppad["available"] is True
    assert capabilities["derivatives"]["cppad"] == {
        **cppad,
        "scope": "package-wide AD substrate; production derivative routes are listed in coverage_matrix",
    }
    ceres = info["native_dependencies"]["ceres"]
    assert ceres["backend"] == "ceres"
    assert ceres["status"] == "enabled_available"
    assert ceres["compiled"] is True
    assert ceres["available"] is True
    assert capabilities["optimizers"]["ceres"]["available"] is ceres["available"]
    assert capabilities["optimizers"]["ceres"]["production"] is ceres["available"]
    assert capabilities["optimizers"]["ceres"]["native_hot_loop"] is ceres["available"]
    equilibrium = capabilities["equilibrium"]
    assert "reactive_electrolyte_bubble" not in equilibrium
    assert "reactive_phase_equilibrium" not in equilibrium
    unpromoted_electrolyte_vle_routes = {
        "electrolyte_bubble_temperature",
        "electrolyte_dew_pressure",
        "electrolyte_dew_temperature",
        "electrolyte_tp_flash",
        "electrolyte_vle",
    }
    assert unpromoted_electrolyte_vle_routes.isdisjoint(equilibrium)
    assert unpromoted_electrolyte_vle_routes.isdisjoint(ipopt["public_routes"])

    neutral_tp = equilibrium["neutral_tp_flash"]
    assert neutral_tp["available"] is ipopt["available"]
    assert neutral_tp["backend"] == "native_ipopt_equilibrium_nlp"
    assert neutral_tp["methods"] == ["tp_flash", "flash_tp"]
    neutral_stability = equilibrium["neutral_stability"]
    assert neutral_stability["available"] is ipopt["available"]
    assert neutral_stability["backend"] == "native_ipopt_equilibrium_nlp"
    assert neutral_stability["methods"] == ["stability_tp"]
    assert neutral_stability["route"] == "tangent_plane_distance"
    electrolyte_stability = equilibrium["electrolyte_stability"]
    assert electrolyte_stability["available"] is ipopt["available"]
    assert electrolyte_stability["backend"] == "native_ipopt_equilibrium_nlp"
    assert electrolyte_stability["methods"] == ["electrolyte_stability_tp"]
    assert electrolyte_stability["route"] == "charge_constrained_tangent_plane_distance"
    neutral_bubble_dew = equilibrium["neutral_bubble_dew"]
    assert neutral_bubble_dew["available"] is ipopt["available"]
    assert neutral_bubble_dew["backend"] == "native_ipopt_equilibrium_nlp"
    assert neutral_bubble_dew["methods"] == ["bubble_p", "bubble_t", "dew_p", "dew_t"]
    neutral_lle = equilibrium["neutral_lle_flash"]
    assert neutral_lle["available"] is ipopt["available"]
    assert neutral_lle["backend"] == "native_ipopt_equilibrium_nlp"
    assert neutral_lle["methods"] == ["lle_flash", "lle_tp"]
    electrolyte_bubble = equilibrium["electrolyte_bubble_pressure"]
    assert electrolyte_bubble["available"] is ipopt["available"]
    assert electrolyte_bubble["backend"] == "native_ipopt_equilibrium_nlp"
    assert electrolyte_bubble["scope"] == "fixed liquid composition with neutral vapor species; ions remain liquid-only"
    electrolyte_lle = equilibrium["electrolyte_lle"]
    assert electrolyte_lle["available"] is ipopt["available"]
    assert electrolyte_lle["backend"] == "native_ipopt_equilibrium_nlp"
    assert electrolyte_lle["methods"] == ["electrolyte_lle", "electrolyte_lle_tp"]
    reactive_speciation = equilibrium["reactive_speciation"]
    assert reactive_speciation["available"] is ipopt["available"]
    assert reactive_speciation["backend"] == "native_ipopt_equilibrium_nlp"
    assert reactive_speciation["sweep_available"] is ipopt["available"]
    assert reactive_speciation["solver_backends"] == ["auto", "ipopt"]
    assert reactive_speciation["ipopt_routes"] == [
        "reactive_speciation:ideal_mole_fraction",
        "reactive_speciation:mole_fraction_activity",
        "reactive_speciation:concentration",
    ]
    assert reactive_speciation["ideal_speciation_nlp_available"] is ipopt["available"]
    assert reactive_speciation["nonideal_speciation_nlp_available"] is ipopt["available"]
    assert reactive_speciation["implemented_standard_states"] == [
        "ideal_mole_fraction",
        "mole_fraction_activity",
        "concentration",
    ]
    assert (
        reactive_speciation["jacobian_auto_policy"]
        == "ideal_analytic_nonideal_cppad_explicit_density_else_raise"
    )
    assert reactive_speciation["jacobian_auto_supported_standard_states"] == [
        "ideal_mole_fraction",
        "mole_fraction_activity",
        "concentration",
    ]
    assert reactive_speciation["auto_request"] == "implemented_standard_states_route_to_native_ipopt"
    assert reactive_speciation["nonideal_derivative_backend"] == "cppad_explicit_density"
    assert reactive_speciation["mixed_standard_state_policy"] == "raise_until_single_objective_is_specified"
    assert capabilities["regression"]["pure_neutral"]["backend"] == "native_ceres"
    reactive_regression = capabilities["regression"]["reactive_electrolyte_residuals"]
    assert reactive_regression["available"] is True
    assert reactive_regression["backend"] == "structured_residual_evaluation"
    assert "not a production optimizer" in reactive_regression["scope"]
    batch_context = capabilities["regression"]["reactive_electrolyte_batch_context"]
    assert batch_context["available"] is True
    assert batch_context["backend"] == "batch_residual_evaluation_context"
    assert "ReactiveElectrolyteRegressionContext" in batch_context["classes"]
    assert batch_context["methods"] == ["evaluate_objective"]
    assert equilibrium["problem_objects"]["entrypoint"] == "mixture.solve_equilibrium(problem)"
    assert "ReactivePhaseEquilibriumProblem" in equilibrium["problem_objects"]["classes"]
    assert equilibrium["repeated_state_properties"]["density_seed_parameter"] == "rho_guess"

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

def test_fugacity_helper_accepts_canonical_rho_guess() -> None:
    state, _ = _ionic_state()
    mix = state.mixture
    seeded = epcsaft.evaluate_fugacity_coefficients(
        mix,
        T=state.T,
        x=state.x,
        P=state.pressure(),
        phase="liq",
        rho_guess=state.molar_density(),
    )

    assert seeded["density"] == pytest.approx(state.molar_density())

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
