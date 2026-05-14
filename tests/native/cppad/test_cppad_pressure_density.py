from __future__ import annotations

import numpy as np
import pytest

from epcsaft import _core
from epcsaft.epcsaft import create_struct


def _single_component_args():
    return create_struct({"m": np.asarray([1.0]), "s": np.asarray([3.7039]), "e": np.asarray([150.03])})


def test_cppad_pressure_density_closure_records_density_dependence() -> None:
    args = _single_component_args()
    t = 300.0
    rho = 12.5

    result = _core._native_cppad_pressure_density(t, rho, [1.0], args)

    if not result["cppad_compiled"]:
        assert result["derivative_backend"] == "backend_unavailable"
        return

    assert result["cppad_used"] is True
    assert result["derivative_backend"] == "cppad"
    assert result["shape"] == (1, 1)
    assert result["value"][0] == pytest.approx(1.380649e-23 * t * rho * 6.02214076e23)
    assert result["jacobian_row_major"][0] == pytest.approx(1.380649e-23 * t * 6.02214076e23)
