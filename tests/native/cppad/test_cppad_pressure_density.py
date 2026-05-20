from __future__ import annotations

import numpy as np
import pytest

from epcsaft import _core
from epcsaft import ePCSAFTMixture
from epcsaft.epcsaft import create_struct


def _single_component_args():
    return create_struct({"m": np.asarray([1.0]), "s": np.asarray([3.7039]), "e": np.asarray([150.03])})


def test_cppad_pressure_density_closure_records_density_dependence() -> None:
    params = {"m": np.asarray([1.0]), "s": np.asarray([3.7039]), "e": np.asarray([150.03])}
    args = create_struct(params)
    t = 300.0
    rho = 1000.0

    result = _core._native_cppad_pressure_density(t, rho, [1.0], args)

    assert result["cppad_used"] is True
    assert result["derivative_backend"] == "cppad"
    assert result["outputs"] == ["pressure"]
    assert result["variables"] == ["rho"]
    assert result["shape"] == (1, 1)
    mix = ePCSAFTMixture.from_params(params, species=["A"])
    state = mix.state(T=t, x=np.asarray([1.0]), rho=rho)
    step = rho * 1.0e-4
    slope = (
        mix.state(T=t, x=np.asarray([1.0]), rho=rho + step).pressure()
        - mix.state(T=t, x=np.asarray([1.0]), rho=rho - step).pressure()
    ) / (2.0 * step)
    assert result["value"][0] == pytest.approx(state.pressure())
    assert result["jacobian_row_major"][0] == pytest.approx(slope, rel=5.0e-4)
