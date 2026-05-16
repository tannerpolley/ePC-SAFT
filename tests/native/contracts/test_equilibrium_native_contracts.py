from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

import epcsaft
from epcsaft import _core, ePCSAFTMixture
from tests.equilibrium.core.test_stability import _assert_stability_route_pending
from tests.equilibrium.core.test_vle import _assert_tp_flash_route_pending
from tests.equilibrium.electrolyte.test_electrolyte_lle_smokes import _assert_electrolyte_lle_route_pending

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
    assert hasattr(_core, "_evaluate_electrolyte_lle_residual_native")


def test_public_tp_flash_requires_native_ipopt_route() -> None:
    mix = _hydrocarbon_mixture()

    with pytest.raises(epcsaft.InputError) as excinfo:
        mix.equilibrium(kind="tp_flash", T=220.0, P=1.0e5, z=[0.1, 0.3, 0.6], backend="native")

    _assert_tp_flash_route_pending(excinfo)


def test_public_electrolyte_stability_requires_native_ipopt_route() -> None:
    mix = _electrolyte_mixture()

    with pytest.raises(epcsaft.InputError) as excinfo:
        mix.equilibrium(
            kind="electrolyte_stability",
            T=298.15,
            P=1.013e5,
            z=[0.55, 0.40, 0.025, 0.025],
            backend="native",
            options=epcsaft.EquilibriumOptions(max_iterations=60, tolerance=1.0e-8),
        )

    _assert_stability_route_pending(excinfo, route="electrolyte_stability")


def test_public_electrolyte_lle_requires_native_ipopt_route() -> None:
    mix = _electrolyte_mixture()
    aq = np.asarray([0.798324680201737, 0.016320352824141723, 0.09267748348706063, 0.09267748348706063], dtype=float)
    org = np.asarray([0.37006036048879404, 0.6214918588210971, 0.004223890345054407, 0.004223890345054407], dtype=float)
    beta_org = 0.613766575013417
    feed = (1.0 - beta_org) * aq + beta_org * org

    with pytest.raises(epcsaft.InputError) as excinfo:
        mix.equilibrium(
            kind="electrolyte_lle",
            T=298.15,
            P=1.013e5,
            z=feed,
            backend="native",
            initial_phases={"aq": aq, "org": org, "phase_fraction": beta_org},
            options=epcsaft.EquilibriumOptions(max_iterations=80, tolerance=1.0e-8),
        )

    _assert_electrolyte_lle_route_pending(excinfo)


def test_native_electrolyte_lle_residual_evaluator_reports_cppad_implicit_derivatives() -> None:
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

    payload = _core._evaluate_electrolyte_lle_residual_native(mix._native, request)

    assert payload["variable_model"] == "ascani_transformed_salt_pairs"
    assert payload["jacobian_backend"] == "cppad_implicit"
    assert payload["jacobian_row_major"]
    assert payload["gradient"]
    assert payload["diagnostics"]["jacobian_available"] is True
    assert payload["diagnostics"]["derivative_available"] is True
    jacobian_rows, jacobian_cols = payload["jacobian_shape"]
    assert len(payload["jacobian_row_major"]) == jacobian_rows * jacobian_cols
    assert len(payload["gradient"]) == jacobian_cols
    assert payload["diagnostics"]["residual_surface"] == "native_electrolyte_lle_transformed_variables"
    assert payload["material_balance_error"] <= 1.0e-10
    assert payload["charge_balance_error"] <= 1.0e-8


def test_native_electrolyte_lle_residual_evaluator_defaults_to_cppad_implicit_derivatives() -> None:
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

    payload = _core._evaluate_electrolyte_lle_residual_native(mix._native, request)

    removed_status = "not" + "_available"
    payload_text = str(payload).lower()
    assert payload["jacobian_backend"] == "cppad_implicit"
    assert payload["diagnostics"]["derivative_available"] is True
    assert removed_status not in payload_text
    assert "numerical_derivative" not in payload_text


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
