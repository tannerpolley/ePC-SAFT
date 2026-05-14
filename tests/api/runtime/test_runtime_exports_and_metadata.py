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
    assert mixed_regression["ceres"]["production"] is False
    assert ipopt["available"] is ipopt_backend.cyipopt_available()
    assert ipopt["formulations"] == ["bound_constrained_residual_minimization"]
    assert ipopt["full_constrained_nlp_available"] is False
    assert ipopt["default_auto_uses_ipopt"] is False
    assert ipopt["exact_hessian_available"] is False
    assert info["optional_dependencies"]["cyipopt"]["available"] is ipopt["available"]
    cppad = info["optional_dependencies"]["cppad"]
    assert cppad["backend"] == "cppad"
    assert cppad["status"] in {"disabled", "enabled_available", "enabled_missing", "not_configured"}
    assert cppad["compiled"] is (cppad["status"] == "enabled_available")
    assert capabilities["derivatives"]["finite_difference"] == {
        "available": False,
        "production": False,
        "reason": "finite_difference_derivatives_forbidden",
    }
    assert capabilities["derivatives"]["eigen_forward"]["scope"] == "legacy/local forward-mode AD"
    assert capabilities["derivatives"]["eigen_forward"]["available"] is True
    assert capabilities["derivatives"]["cppad"] == {
        **cppad,
        "production": False,
        "reason": "dependency_not_compiled" if not cppad["available"] else "not_validated_for_production",
        "scope": "package-wide AD substrate",
        "production_eos_coverage": False,
    }
    ceres = info["optional_dependencies"]["ceres"]
    assert ceres["backend"] == "ceres"
    assert ceres["status"] in {"disabled", "enabled_available", "not_configured"}
    assert ceres["compiled"] is (ceres["status"] == "enabled_available")
    assert capabilities["optimizers"]["ceres"]["available"] is ceres["available"]
    assert capabilities["optimizers"]["ceres"]["production"] is False
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
        == "analytic_inner_solve_with_activity_fixed_point_for_nonideal_standard_states"
    )
    assert (
        capabilities["equilibrium"]["reactive_speciation"]["derivative_gap_status"]
        == "activity_derivatives_not_used_by_fixed_point_outer_iteration"
    )
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
