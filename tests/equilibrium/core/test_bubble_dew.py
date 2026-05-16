from __future__ import annotations

import numpy as np
import pytest

import epcsaft


def _hydrocarbon_mixture() -> epcsaft.ePCSAFTMixture:
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
    return epcsaft.ePCSAFTMixture.from_params(params, species=["Methane", "Ethane", "Propane"])


@pytest.mark.parametrize(
    ("method", "kwargs", "route"),
    [
        ("bubble_p", {"T": 220.0, "x": [0.2, 0.3, 0.5]}, "bubble_p"),
        ("dew_p", {"T": 220.0, "y": [0.1, 0.3, 0.6]}, "dew_p"),
        ("bubble_t", {"P": 1.0e5, "x": [0.2, 0.3, 0.5]}, "bubble_t"),
        ("dew_t", {"P": 1.0e5, "y": [0.1, 0.3, 0.6]}, "dew_t"),
    ],
)
def test_neutral_bubble_dew_requires_native_ipopt_route(method: str, kwargs: dict[str, object], route: str) -> None:
    mix = _hydrocarbon_mixture()

    with pytest.raises(epcsaft.InputError, match=rf"{route} requires a native Ipopt equilibrium NLP route"):
        getattr(mix, method)(**kwargs)
