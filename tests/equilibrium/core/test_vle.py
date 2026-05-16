from __future__ import annotations

import numpy as np
import pytest

import epcsaft
from epcsaft import ePCSAFTMixture


def _hydrocarbon_basis_mixture() -> ePCSAFTMixture:
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


def _assert_tp_flash_route_pending(excinfo: pytest.ExceptionInfo[epcsaft.InputError]) -> None:
    message = str(excinfo.value)
    assert "tp_flash requires a native Ipopt equilibrium NLP route" in message
    assert "previous native TP flash route is disabled" in message


def test_ternary_hydrocarbon_basis_tp_flash_requires_native_ipopt_route() -> None:
    mix = _hydrocarbon_basis_mixture()
    feed = np.asarray([0.1, 0.3, 0.6])

    with pytest.raises(epcsaft.InputError) as excinfo:
        mix.equilibrium(kind="tp_flash", T=220.0, P=1.0e5, z=feed)

    _assert_tp_flash_route_pending(excinfo)


def test_tp_flash_accepts_stability_precheck_option_before_route_gate() -> None:
    mix = _hydrocarbon_basis_mixture()

    with pytest.raises(epcsaft.InputError) as excinfo:
        mix.equilibrium(
            kind="tp_flash",
            T=300.0,
            P=1.0e5,
            z=[0.1, 0.3, 0.6],
            options=epcsaft.EquilibriumOptions(stability_precheck=False),
        )

    _assert_tp_flash_route_pending(excinfo)


def test_tp_flash_phase_diagnostics_option_requires_native_ipopt_route() -> None:
    mix = _hydrocarbon_basis_mixture()

    with pytest.raises(epcsaft.InputError) as excinfo:
        mix.equilibrium(
            kind="tp_flash",
            T=220.0,
            P=1.0e5,
            z=[0.1, 0.3, 0.6],
            options=epcsaft.EquilibriumOptions(include_phase_diagnostics=True),
        )

    _assert_tp_flash_route_pending(excinfo)
