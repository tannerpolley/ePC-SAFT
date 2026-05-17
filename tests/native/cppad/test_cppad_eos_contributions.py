from __future__ import annotations

import numpy as np
import pytest

from epcsaft import _core
from epcsaft.epcsaft import create_struct


def _neutral_args():
    return create_struct(
        {
            "m": np.asarray([1.0, 1.6069]),
            "s": np.asarray([3.7039, 3.5206]),
            "e": np.asarray([150.03, 191.42]),
            "k_ij": np.zeros((2, 2)),
        }
    )


def test_cppad_eos_contribution_recording_matches_double_value_path() -> None:
    args = _neutral_args()
    x = [0.35, 0.65]
    t = 310.0
    rho = 8200.0

    result = _core._native_cppad_eos_contributions(t, rho, x, args)

    state = _core.NativeState(_core.NativeMixture(args), t, x, 0, False, 0.0, True, rho, False, 0.0)
    expected = state.residual_helmholtz_result()

    assert result["derivative_backend"] == "cppad"
    assert result["cppad_used"] is True
    assert result["outputs"] == ["hc", "disp", "assoc", "ion", "born", "total"]
    assert result["variables"] == ["x_0", "x_1"]
    assert result["shape"] == (6, 2)
    assert result["value"] == pytest.approx(
        [expected.hc, expected.disp, expected.assoc, expected.ion, expected.born, expected.total]
    )

    jacobian = np.asarray(result["jacobian_row_major"], dtype=float).reshape(result["shape"])
    assert np.all(np.isfinite(jacobian))
    assert jacobian[-1] == pytest.approx(jacobian[:-1].sum(axis=0))


def test_cppad_eos_contribution_recording_rejects_active_association() -> None:
    args = _neutral_args()
    args.assoc_num = [1, 1]

    with pytest.raises(_core.NativeValueError, match="unsupported"):
        _core._native_cppad_eos_contributions(310.0, 8200.0, [0.35, 0.65], args)
