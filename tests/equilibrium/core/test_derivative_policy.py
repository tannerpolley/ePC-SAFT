from __future__ import annotations

import numpy as np
import pytest

import epcsaft
from tests.api.reactive.test_reactive_speciation_options import _assert_reactive_speciation_route_pending
from tests.equilibrium.core.test_vle import _assert_tp_flash_route_pending
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


def test_auto_equilibrium_diagnostics_reject_numerical_derivative_route() -> None:
    mix, feed = _neutral_mixture()

    with pytest.raises(epcsaft.InputError) as excinfo:
        mix.flash_tp(T=220.0, P=1.0e5, z=feed, options=epcsaft.EquilibriumOptions(jacobian_backend="auto"))

    _assert_tp_flash_route_pending(excinfo)
    assert "numerical_derivative" not in str(excinfo.value).lower()


def test_reactive_ideal_speciation_auto_requires_native_ipopt_route() -> None:
    mix = epcsaft.ePCSAFTMixture.from_params(
        {
            "m": np.asarray([1.0, 1.0]),
            "s": np.asarray([3.0, 3.0]),
            "e": np.asarray([200.0, 200.0]),
        },
        species=["A", "B"],
    )

    with pytest.raises(epcsaft.InputError) as excinfo:
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
            options=epcsaft.ReactiveSpeciationOptions(jacobian_backend="auto", tolerance=1.0e-10),
        )

    _assert_reactive_speciation_route_pending(excinfo)


def test_activity_coupled_reactive_speciation_auto_requires_native_ipopt_route() -> None:
    mix = epcsaft.ePCSAFTMixture.from_params(
        {
            "m": np.asarray([1.0, 1.0]),
            "s": np.asarray([3.0, 3.0]),
            "e": np.asarray([200.0, 200.0]),
        },
        species=["A", "B"],
    )

    with pytest.raises(epcsaft.InputError) as excinfo:
        mix.chemical_equilibrium(
            T=298.15,
            P=1.0e5,
            balances={"total": {"A": 1.0, "B": 1.0}},
            totals={"total": 1.0},
            reactions=[epcsaft.ReactionDefinition({"A": -1.0, "B": 1.0}, np.log(3.0))],
            initial_x=[0.5, 0.5],
            options=epcsaft.ReactiveSpeciationOptions(jacobian_backend="auto"),
        )

    _assert_reactive_speciation_route_pending(excinfo)


def test_explicit_cppad_reactive_speciation_fails_until_supported() -> None:
    mix = epcsaft.ePCSAFTMixture.from_params(
        {
            "m": np.asarray([1.0, 1.0]),
            "s": np.asarray([3.0, 3.0]),
            "e": np.asarray([200.0, 200.0]),
        },
        species=["A", "B"],
    )

    with pytest.raises(epcsaft.InputError) as excinfo:
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

    _assert_reactive_speciation_route_pending(excinfo)
