from __future__ import annotations

import pytest

from epcsaft._types import InputError
from tests.helpers.native_cases import _ionic_state


def _state():
    mix, _species, _pressure, density, temperature, composition = _ionic_state()
    return mix.state(T=temperature, x=composition, rho=density)


def test_activity_composition_derivative_route_raises_until_supported() -> None:
    state = _state()

    with pytest.raises(InputError, match="unsupported"):
        state.activity_composition_derivative_result()


def test_activity_parameter_derivative_route_raises_without_ssmds_path() -> None:
    state = _state()

    with pytest.raises(InputError, match="unsupported"):
        state.activity_parameter_derivative_result()
