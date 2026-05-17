from __future__ import annotations

import pytest

from epcsaft import _core


def test_native_cppad_smoke_reports_exact_cppad_derivative() -> None:
    smoke = _core._native_cppad_smoke()

    assert smoke["status"] == "enabled_available"
    assert smoke["cppad_compiled"] is True
    assert smoke["cppad_used"] is True
    assert smoke["derivative_backend"] == "cppad"
    assert smoke["outputs"] == ["x_squared"]
    assert smoke["variables"] == ["x"]
    assert smoke["value"] == pytest.approx([9.0])
    assert smoke["jacobian_row_major"] == pytest.approx([6.0])
    assert smoke["shape"] == (1, 1)
