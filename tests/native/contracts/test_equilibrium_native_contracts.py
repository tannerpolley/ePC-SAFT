from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

import epcsaft
from epcsaft import _core, ePCSAFTMixture

REPO_ROOT = Path(__file__).resolve().parents[3]


def _hydrocarbon_mixture() -> ePCSAFTMixture:
    params = {
        "m": np.asarray([1.0, 1.6069, 2.0020]),
        "s": np.asarray([3.7039, 3.5206, 3.6184]),
        "e": np.asarray([150.03, 191.42, 208.11]),
        "k_ij": np.asarray(
            [
                [0.0, 3.0e-4, 1.15e-2],
                [3.0e-4, 0.0, 5.10e-3],
                [1.15e-2, 5.10e-3, 0.0],
            ]
        ),
    }
    return ePCSAFTMixture.from_params(params, species=["Methane", "Ethane", "Propane"])


def _electrolyte_mixture() -> ePCSAFTMixture:
    feed = np.asarray([0.55, 0.40, 0.025, 0.025], dtype=float)
    return ePCSAFTMixture.from_dataset("2022_Ascani", ["H2O", "Butanol", "Na+", "Cl-"], feed, 298.15)


def test_native_equilibrium_entrypoint_is_exposed() -> None:
    assert hasattr(_core, "_solve_equilibrium_native")
    assert hasattr(_core, "_evaluate_electrolyte_lle_residual_native")


def test_public_equilibrium_result_comes_from_native_backend() -> None:
    mix = _hydrocarbon_mixture()

    result = mix.equilibrium(kind="tp_flash", T=220.0, P=1.0e5, z=[0.1, 0.3, 0.6], backend="native")

    assert isinstance(result, epcsaft.EquilibriumResult)
    assert result.backend == "neutral_vle"
    assert result.diagnostics["solver_language"] == "c++"
    assert result.diagnostics["native_entrypoint"] == "_solve_equilibrium_native"


def test_native_electrolyte_stability_entrypoint_runs_in_cpp() -> None:
    mix = _electrolyte_mixture()
    request = {
        "kind": "electrolyte_stability",
        "T": 298.15,
        "P": 1.013e5,
        "z": [0.55, 0.40, 0.025, 0.025],
        "species": mix.species,
        "options": {
            "max_iterations": 60,
            "tolerance": 1.0e-8,
            "damping": 0.5,
            "min_composition": 1.0e-12,
            "include_phase_diagnostics": False,
            "stability_precheck": True,
        },
    }

    payload = _core._solve_equilibrium_native(mix._native, request)

    assert payload["result_type"] == "stability"
    assert payload["backend"] == "electrolyte_tpd"
    assert payload["diagnostics"]["solver_language"] == "c++"
    assert payload["diagnostics"]["native_entrypoint"] == "_solve_equilibrium_native"
    assert payload["diagnostics"]["tpd_method"] == "native_tpd_global_search"
    assert payload["diagnostics"]["tpd_trial_count"] == len(payload["trials"])
    assert payload["diagnostics"]["tpd_multistart_count"] > 0
    assert payload["diagnostics"]["tpd_polish_iterations"] > 0
    assert payload["diagnostics"]["tpd_best_seed_name"] == payload["diagnostics"]["seed_name"]
    assert payload["diagnostics"]["phase_charge_balance"]["trial"] == pytest.approx(0.0, abs=1.0e-8)


def test_native_electrolyte_stability_honors_explicit_max_iterations() -> None:
    mix = _electrolyte_mixture()
    request = {
        "kind": "electrolyte_stability",
        "T": 298.15,
        "P": 1.013e5,
        "z": [0.55, 0.40, 0.025, 0.025],
        "species": mix.species,
        "options": {
            "max_iterations": 2,
            "tolerance": 1.0e-8,
            "damping": 0.5,
            "min_composition": 1.0e-12,
            "include_phase_diagnostics": False,
            "stability_precheck": True,
        },
    }

    payload = _core._solve_equilibrium_native(mix._native, request)

    assert payload["diagnostics"]["requested_max_iterations"] == 2
    assert payload["diagnostics"]["effective_max_iterations"] == 2


def test_public_electrolyte_stability_uses_native_backend() -> None:
    mix = _electrolyte_mixture()

    result = mix.equilibrium(
        kind="electrolyte_stability",
        T=298.15,
        P=1.013e5,
        z=[0.55, 0.40, 0.025, 0.025],
        backend="native",
        options=epcsaft.EquilibriumOptions(max_iterations=60, tolerance=1.0e-8),
    )

    assert result.backend == "electrolyte_tpd"
    assert result.diagnostics["solver_language"] == "c++"
    assert result.diagnostics["native_entrypoint"] == "_solve_equilibrium_native"


def test_public_electrolyte_lle_reports_unavailable_solver_derivatives() -> None:
    mix = _electrolyte_mixture()
    aq = np.asarray([0.798324680201737, 0.016320352824141723, 0.09267748348706063, 0.09267748348706063], dtype=float)
    org = np.asarray([0.37006036048879404, 0.6214918588210971, 0.004223890345054407, 0.004223890345054407], dtype=float)
    beta_org = 0.613766575013417
    feed = (1.0 - beta_org) * aq + beta_org * org

    with pytest.raises(epcsaft.InputError, match="not_available"):
        mix.equilibrium(
            kind="electrolyte_lle",
            T=298.15,
            P=1.013e5,
            z=feed,
            backend="native",
            initial_phases={"aq": aq, "org": org, "phase_fraction": beta_org},
            options=epcsaft.EquilibriumOptions(max_iterations=80, tolerance=1.0e-8),
        )


def test_native_electrolyte_lle_residual_evaluator_reports_unavailable_derivatives() -> None:
    mix = _electrolyte_mixture()
    aq = np.asarray([0.798324680201737, 0.016320352824141723, 0.09267748348706063, 0.09267748348706063], dtype=float)
    org = np.asarray([0.37006036048879404, 0.6214918588210971, 0.004223890345054407, 0.004223890345054407], dtype=float)
    beta_org = 0.613766575013417
    feed = ((1.0 - beta_org) * aq + beta_org * org).tolist()
    request = {
        "T": 298.15,
        "P": 1.013e5,
        "z": feed,
        "species": mix.species,
        "initial_phases": {"aq": aq.tolist(), "org": org.tolist(), "phase_fraction": beta_org},
        "options": {
            "max_iterations": 80,
            "tolerance": 1.0e-8,
            "min_composition": 1.0e-12,
            "jacobian_backend": "auto",
        },
    }

    with pytest.raises(_core.NativeValueError, match="not_available"):
        _core._evaluate_electrolyte_lle_residual_native(mix._native, request)


def test_native_electrolyte_lle_residual_evaluator_rejects_auto_without_autodiff() -> None:
    mix = _electrolyte_mixture()
    aq = np.asarray([0.798324680201737, 0.016320352824141723, 0.09267748348706063, 0.09267748348706063], dtype=float)
    org = np.asarray([0.37006036048879404, 0.6214918588210971, 0.004223890345054407, 0.004223890345054407], dtype=float)
    beta_org = 0.613766575013417
    feed = ((1.0 - beta_org) * aq + beta_org * org).tolist()
    request = {
        "T": 298.15,
        "P": 1.013e5,
        "z": feed,
        "species": mix.species,
        "initial_phases": {"aq": aq.tolist(), "org": org.tolist(), "phase_fraction": beta_org},
        "options": {"max_iterations": 80, "tolerance": 1.0e-8, "min_composition": 1.0e-12},
    }

    with pytest.raises(_core.NativeValueError, match="not_available"):
        _core._evaluate_electrolyte_lle_residual_native(mix._native, request)


def test_equilibrium_runtime_does_not_import_external_optimizers() -> None:
    source = (REPO_ROOT / "src" / "epcsaft" / "equilibrium.py").read_text(encoding="utf-8")

    external_optimizer = "sci" + "py.optimize"
    forbidden = (
        external_optimizer,
        "least_squares",
        "differential_evolution",
        "minimize_scalar",
        "_solve_predictive_electrolyte",
        "_electrolyte_gibbs_seed_from_trial",
        "_electrolyte_stability_from_basis",
    )

    for token in forbidden:
        assert token not in source


def test_package_runtime_has_no_external_optimizer_dependency_or_imports() -> None:
    pyproject = (REPO_ROOT / "pyproject.toml").read_text(encoding="utf-8")
    docs_requirements = (REPO_ROOT / "docs" / "requirements.txt").read_text(encoding="utf-8")
    package_sources = "\n".join(
        path.read_text(encoding="utf-8")
        for path in (REPO_ROOT / "src" / "epcsaft").rglob("*.py")
        if path.name != "__pycache__"
    )

    dependency_token = "sci" + "py"
    assert dependency_token not in pyproject.lower()
    assert dependency_token not in docs_requirements.lower()
    assert f"from {dependency_token}" not in package_sources
    assert f"import {dependency_token}" not in package_sources


def test_public_equilibrium_does_not_expose_python_backend_tokens() -> None:
    source = (REPO_ROOT / "src" / "epcsaft" / "epcsaft.py").read_text(encoding="utf-8")
    equilibrium_source = (REPO_ROOT / "src" / "epcsaft" / "equilibrium.py").read_text(encoding="utf-8")

    assert '"python"' not in source
    assert "Python-first" not in equilibrium_source
    assert "np.linalg.lstsq" not in equilibrium_source
