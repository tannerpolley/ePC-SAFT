from __future__ import annotations

import numpy as np
import pytest

import epcsaft
from tests.helpers.native_cases import _neutral_state


def _neutral_mixture() -> tuple[epcsaft.ePCSAFTMixture, np.ndarray]:
    mix, _species, _pressure, _density, _temperature, composition = _neutral_state()
    return mix, composition


def test_equilibrium_options_reject_removed_numerical_derivative_backend() -> None:
    mix, feed = _neutral_mixture()
    removed_backend = "finite" + "_difference"

    with pytest.raises(epcsaft.InputError, match="jacobian_backend"):
        mix.flash_tp(
            T=220.0,
            P=1.0e5,
            z=feed,
            options=epcsaft.EquilibriumOptions(jacobian_backend=removed_backend),
        )


def test_equilibrium_options_reject_legacy_autodiff_backend() -> None:
    mix, feed = _neutral_mixture()

    with pytest.raises(epcsaft.InputError, match="jacobian_backend"):
        mix.flash_tp(
            T=220.0,
            P=1.0e5,
            z=feed,
            options=epcsaft.EquilibriumOptions(jacobian_backend="autodiff"),
        )


def test_auto_equilibrium_diagnostics_do_not_fallback_to_numerical_derivative() -> None:
    mix, feed = _neutral_mixture()

    result = mix.flash_tp(T=220.0, P=1.0e5, z=feed, options=epcsaft.EquilibriumOptions(jacobian_backend="auto"))

    diagnostics = result.diagnostics
    assert diagnostics["requested_jacobian_backend"] == "auto"
    assert diagnostics["jacobian_fallback_used"] is False
    assert diagnostics["derivative_backend"] == "not_applicable"
    assert diagnostics["derivative_status"] == "not_applicable"
    removed_reason_key = "not" + "_available_reason"
    assert removed_reason_key not in diagnostics
    assert diagnostics["thermodynamic_backend"] == "epcsaft_state_fugacity_activity_property_api"
    assert "density_roots" in diagnostics["solved_internal_states"]
    assert diagnostics["derivative_backend_by_block"]["density_root"] == "not_applicable"
    assert "numerical_derivative" not in str(diagnostics).lower()


def test_reactive_ideal_speciation_uses_analytic_residual_derivatives() -> None:
    mix = epcsaft.ePCSAFTMixture.from_params(
        {
            "m": np.asarray([1.0, 1.0]),
            "s": np.asarray([3.0, 3.0]),
            "e": np.asarray([200.0, 200.0]),
        },
        species=["A", "B"],
    )

    result = mix.chemical_equilibrium(
        T=298.15,
        P=1.0e5,
        balances={"total": {"A": 1.0, "B": 1.0}},
        totals={"total": 1.0},
        reactions=[
            epcsaft.ReactionDefinition(
                {"A": -1.0, "B": 1.0},
                np.log(3.0),
                standard_state="ideal_mole_fraction",
            )
        ],
        initial_x=[0.5, 0.5],
        options=epcsaft.ReactiveSpeciationOptions(jacobian_backend="auto", tolerance=1.0e-10),
    )

    diagnostics = result.diagnostics
    assert diagnostics["derivative_backend"] == "analytic"
    assert diagnostics["derivative_status"] == "analytic"
    assert diagnostics["derivative_backend_by_block"]["reaction_residual_jacobian"] == "analytic"
    assert diagnostics["jacobian_fallback_used"] is False
    backend_values = {str(value).lower() for value in diagnostics["derivative_backend_by_block"].values()}
    assert "numerical_derivative" not in backend_values


def test_activity_coupled_reactive_speciation_uses_native_activity_residual() -> None:
    mix = epcsaft.ePCSAFTMixture.from_params(
        {
            "m": np.asarray([1.0, 1.0]),
            "s": np.asarray([3.0, 3.0]),
            "e": np.asarray([200.0, 200.0]),
        },
        species=["A", "B"],
    )

    result = mix.chemical_equilibrium(
        T=298.15,
        P=1.0e5,
        balances={"total": {"A": 1.0, "B": 1.0}},
        totals={"total": 1.0},
        reactions=[epcsaft.ReactionDefinition({"A": -1.0, "B": 1.0}, np.log(3.0))],
        initial_x=[0.5, 0.5],
        options=epcsaft.ReactiveSpeciationOptions(jacobian_backend="auto"),
    )

    assert result.success is True
    assert result.diagnostics["solver_language"] == "c++"
    assert result.diagnostics["activity_fixed_point"] is False
    assert result.diagnostics["activity_or_fugacity_terms_in_residual"] is True
    assert result.diagnostics["derivative_backend_by_block"]["reactive_speciation_variables"] == "analytic_implicit"


def test_explicit_cppad_reactive_speciation_fails_until_supported() -> None:
    mix = epcsaft.ePCSAFTMixture.from_params(
        {
            "m": np.asarray([1.0, 1.0]),
            "s": np.asarray([3.0, 3.0]),
            "e": np.asarray([200.0, 200.0]),
        },
        species=["A", "B"],
    )

    with pytest.raises(epcsaft.InputError, match="CppAD chemical-equilibrium residual jacobian"):
        mix.chemical_equilibrium(
            T=298.15,
            P=1.0e5,
            balances={"total": {"A": 1.0, "B": 1.0}},
            totals={"total": 1.0},
            reactions=[
                epcsaft.ReactionDefinition(
                    {"A": -1.0, "B": 1.0},
                    np.log(3.0),
                    standard_state="ideal_mole_fraction",
                )
            ],
            initial_x=[0.5, 0.5],
            options=epcsaft.ReactiveSpeciationOptions(jacobian_backend="cppad"),
        )
