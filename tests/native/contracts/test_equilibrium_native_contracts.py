from __future__ import annotations

from pathlib import Path

import numpy as np

import epcsaft
from epcsaft import _core, ePCSAFTMixture

REPO_ROOT = Path(__file__).resolve().parents[3]


def _electrolyte_mixture() -> ePCSAFTMixture:
    feed = np.asarray([0.55, 0.40, 0.025, 0.025], dtype=float)
    return ePCSAFTMixture.from_dataset("2022_Ascani", ["H2O", "Butanol", "Na+", "Cl-"], feed, 298.15)


def test_native_equilibrium_entrypoint_is_exposed() -> None:
    assert hasattr(_core, "_evaluate_electrolyte_lle_residual_native")


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

    assert payload["jacobian_backend"] == "cppad_implicit"
    assert payload["diagnostics"]["derivative_available"] is True


def test_equilibrium_runtime_does_not_import_external_optimizers() -> None:
    source = (REPO_ROOT / "src" / "epcsaft" / "equilibrium.py").read_text(encoding="utf-8")

    external_optimizer = "sci" + "py.optimize"
    forbidden = (
        external_optimizer,
        "least" + "_squares",
        "differential" + "_evolution",
        "minimize" + "_scalar",
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
    assert "np.linalg." + "lstsq" not in equilibrium_source
