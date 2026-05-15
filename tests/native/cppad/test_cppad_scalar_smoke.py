from __future__ import annotations

import pytest

from epcsaft import _core


def test_native_cppad_smoke_reports_disabled_or_exact_cppad_derivative() -> None:
    smoke = _core._native_cppad_smoke()

    assert smoke["status"] in {"disabled", "enabled_available", "enabled_missing", "not_configured"}
    assert smoke["cppad_compiled"] is (smoke["status"] == "enabled_available")

    if not smoke["cppad_compiled"]:
        assert smoke["cppad_used"] is False
        assert smoke["derivative_backend"] == "not_available"
        return

    assert smoke["cppad_used"] is True
    assert smoke["derivative_backend"] == "cppad"
    assert smoke["outputs"] == ["x_squared"]
    assert smoke["variables"] == ["x"]
    assert smoke["value"] == pytest.approx([9.0])
    assert smoke["jacobian_row_major"] == pytest.approx([6.0])
    assert smoke["shape"] == (1, 1)
