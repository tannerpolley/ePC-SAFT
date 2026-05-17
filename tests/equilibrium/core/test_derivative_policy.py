from __future__ import annotations

import numpy as np
import pytest

import epcsaft
from tests.equilibrium.core.test_vle import _assert_tp_flash_route_pending
from tests.helpers.native_cases import _neutral_state


def _neutral_mixture() -> tuple[epcsaft.ePCSAFTMixture, np.ndarray]:
    mix, _species, _pressure, _density, _temperature, composition = _neutral_state()
    return mix, composition


def test_equilibrium_options_reject_removed_nonexact_derivative_backend() -> None:
    mix, feed = _neutral_mixture()
    removed_backend = "finite" + "_difference"

    with pytest.raises(epcsaft.InputError, match="jacobian_backend"):
        mix.flash_tp(
            T=220.0,
            P=1.0e5,
            z=feed,
            options=epcsaft.EquilibriumOptions(jacobian_backend=removed_backend),
        )


def test_equilibrium_options_reject_removed_autodiff_backend() -> None:
    mix, feed = _neutral_mixture()

    with pytest.raises(epcsaft.InputError, match="jacobian_backend"):
        mix.flash_tp(
            T=220.0,
            P=1.0e5,
            z=feed,
            options=epcsaft.EquilibriumOptions(jacobian_backend="autodiff"),
        )


def test_auto_equilibrium_diagnostics_reject_nonexact_derivative_route() -> None:
    mix, feed = _neutral_mixture()

    with pytest.raises(epcsaft.InputError) as excinfo:
        mix.flash_tp(T=220.0, P=1.0e5, z=feed, options=epcsaft.EquilibriumOptions(jacobian_backend="auto"))

    _assert_tp_flash_route_pending(excinfo)
    assert "numerical" + "_derivative" not in str(excinfo.value).lower()


def test_reactive_ideal_speciation_auto_uses_native_ipopt_route() -> None:
    from epcsaft import _core

    mix = epcsaft.ePCSAFTMixture.from_params(
        {
            "m": np.asarray([1.0, 1.0]),
            "s": np.asarray([3.0, 3.0]),
            "e": np.asarray([200.0, 200.0]),
        },
        species=["A", "B"],
    )
    kwargs = {
        "T": 298.15,
        "P": 1.0e5,
        "balances": {"total": {"A": 1.0, "B": 1.0}},
        "totals": {"total": 1.0},
        "reactions": [
            epcsaft.ReactionDefinition(
                {"A": -1.0, "B": 1.0},
                np.log(3.0),
                standard_state="ideal_mole_fraction",
            )
        ],
        "initial_x": [0.5, 0.5],
        "options": epcsaft.ReactiveSpeciationOptions(jacobian_backend="auto", tolerance=1.0e-10),
    }

    if not _core._native_ipopt_smoke()["compiled"]:
        with pytest.raises(epcsaft.SolutionError, match=r"EPCSAFT_ENABLE_IPOPT=ON"):
            mix.chemical_equilibrium(**kwargs)
        return

    result = mix.chemical_equilibrium(**kwargs)

    assert result.success is True
    assert result.diagnostics["requested_solver_backend"] == "auto"
    assert result.diagnostics["selected_solver_backend"] == "native_ipopt"


def test_explicit_cppad_ideal_reactive_speciation_routes_to_native_ipopt() -> None:
    from epcsaft import _core

    mix = epcsaft.ePCSAFTMixture.from_params(
        {
            "m": np.asarray([1.0, 1.0]),
            "s": np.asarray([3.0, 3.0]),
            "e": np.asarray([200.0, 200.0]),
        },
        species=["A", "B"],
    )
    kwargs = {
        "T": 298.15,
        "P": 1.0e5,
        "balances": {"total": {"A": 1.0, "B": 1.0}},
        "totals": {"total": 1.0},
        "reactions": [
            epcsaft.ReactionDefinition(
                {"A": -1.0, "B": 1.0},
                np.log(3.0),
                standard_state="ideal_mole_fraction",
            )
        ],
        "initial_x": [0.5, 0.5],
        "options": epcsaft.ReactiveSpeciationOptions(jacobian_backend="cppad"),
    }

    if not _core._native_ipopt_smoke()["compiled"]:
        with pytest.raises(epcsaft.SolutionError, match=r"EPCSAFT_ENABLE_IPOPT=ON"):
            mix.chemical_equilibrium(**kwargs)
        return

    result = mix.chemical_equilibrium(**kwargs)

    assert result.success is True
    assert result.diagnostics["requested_jacobian_backend"] == "cppad"
    assert result.diagnostics["derivative_backend"] == "cppad"
    assert result.diagnostics["implicit_solve_results"]["reactive_speciation_variables"]["backend"] == "cppad_implicit"
