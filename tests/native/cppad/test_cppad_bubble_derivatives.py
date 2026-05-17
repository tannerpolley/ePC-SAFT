from __future__ import annotations

import numpy as np
import pytest

import epcsaft
from tests.helpers.native_cases import _neutral_state


def test_cppad_pressure_derivative_api_underpins_bubble_policy() -> None:
    mix, _species, _pressure, density, temperature, composition = _neutral_state()
    state = mix.state(T=temperature, x=composition, rho=density)

    pressure_density = state.pressure_density_derivative_result()

    assert pressure_density["backend"] == "cppad"
    assert pressure_density["shape"] == [1, 1]


def test_neutral_bubble_uses_native_ipopt_route_gate() -> None:
    mix, _species, _pressure, _density, _temperature, composition = _neutral_state()

    with pytest.raises(epcsaft.InputError, match=r"bubble_p requires a native Ipopt equilibrium NLP route"):
        mix.bubble_p(T=220.0, x=composition)
