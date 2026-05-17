from __future__ import annotations

import numpy as np
import pytest

import epcsaft
from epcsaft import ePCSAFTMixture


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


def _assert_stability_route_pending(excinfo: pytest.ExceptionInfo[epcsaft.InputError], route: str = "stability") -> None:
    message = str(excinfo.value)
    assert f"{route} requires a native Ipopt equilibrium stability NLP route" in message


def test_stability_public_exports_are_available() -> None:
    assert hasattr(epcsaft, "StabilityTrial")
    assert hasattr(epcsaft, "StabilityResult")


def test_stability_requires_native_ipopt_route_after_validation() -> None:
    mix = _hydrocarbon_mixture()

    with pytest.raises(epcsaft.InputError) as excinfo:
        mix.equilibrium(
            kind="stability",
            T=300.0,
            P=1.0e5,
            z=[0.1, 0.3, 0.6],
            parent_phase="vap",
            trial_phases=("vap",),
        )

    _assert_stability_route_pending(excinfo)


@pytest.mark.parametrize(
    ("kwargs", "match"),
    [
        ({"kind": "stability", "P": 1.0e5, "z": [0.1, 0.3, 0.6]}, "T"),
        ({"kind": "stability", "T": 300.0, "z": [0.1, 0.3, 0.6]}, "P"),
        ({"kind": "stability", "T": 300.0, "P": 1.0e5}, "z"),
        ({"kind": "stability", "T": 300.0, "P": 1.0e5, "z": [1.0]}, "length"),
        ({"kind": "stability", "T": 300.0, "P": 1.0e5, "z": [0.1, -0.3, 0.6]}, "non-negative"),
        (
            {"kind": "stability", "T": 300.0, "P": 1.0e5, "z": [0.1, 0.3, 0.6], "parent_phase": "solid"},
            "parent_phase",
        ),
        (
            {"kind": "stability", "T": 300.0, "P": 1.0e5, "z": [0.1, 0.3, 0.6], "trial_phases": ("liq", "solid")},
            "trial_phases",
        ),
    ],
)
def test_stability_rejects_invalid_public_inputs(kwargs, match) -> None:
    mix = _hydrocarbon_mixture()

    with pytest.raises(epcsaft.InputError, match=match):
        mix.equilibrium(**kwargs)


def test_stability_rejects_ionic_mixtures_for_v3() -> None:
    params = {
        "m": np.asarray([1.2047, 1.0, 1.0]),
        "s": np.asarray([2.7927, 2.8232, 2.7560]),
        "e": np.asarray([353.95, 230.0, 170.0]),
        "z": np.asarray([0.0, 1.0, -1.0]),
        "dielc": np.asarray([78.09, 8.0, 8.0]),
    }
    mix = ePCSAFTMixture.from_params(params, species=["water", "Na+", "Cl-"])

    with pytest.raises(epcsaft.InputError, match="ion-containing"):
        mix.equilibrium(kind="stability", T=298.15, P=1.0e5, z=[0.9998, 1.0e-4, 1.0e-4])
