from __future__ import annotations

import numpy as np

from epcsaft import _core
from tests.helpers.numeric import assert_allclose
from tests.native.equilibrium.test_electrolyte_lle_residual_surface import _electrolyte_mixture, _initial_request


def test_electrolyte_lle_residual_jacobian_reports_real_transformed_surface() -> None:
    mix = _electrolyte_mixture()
    payload = _core._evaluate_electrolyte_lle_residual_native(mix._native, _initial_request(mix))
    diagnostics = payload["diagnostics"]

    rows, cols = payload["jacobian_shape"]
    jacobian = np.asarray(payload["jacobian_row_major"], dtype=float).reshape((rows, cols))
    residual = np.asarray(payload["residual"], dtype=float)
    gradient = np.asarray(payload["gradient"], dtype=float)

    assert rows == residual.size
    assert cols == len(payload["variables"])
    assert diagnostics["jacobian_available"] is True
    assert diagnostics["derivative_available"] is True
    assert payload["jacobian_backend"] == "cppad_implicit"
    assert diagnostics["jacobian_scope"] == "transformed_variables_phase_state_implicit_density"
    assert np.all(np.isfinite(jacobian))
    assert np.any(np.abs(jacobian[: diagnostics["phase_equilibrium_residual_size"]]) > 1.0e-12)
    assert_allclose(gradient, jacobian.T @ residual, rtol=1.0e-12, atol=1.0e-12)
    assert "numerical_derivative" not in str(payload).lower()
