from __future__ import annotations

import numpy as np
import pytest

from epcsaft import _core
from epcsaft.epcsaft import create_struct
from tests.helpers.native_cases import _ionic_state


def _state():
    mix, _species, _pressure, density, temperature, composition = _ionic_state()
    return mix.state(T=temperature, x=composition, rho=density)


def test_active_association_reports_implicit_backend_not_direct_cppad() -> None:
    state = _state()

    association_rows = [row for row in state.derivative_coverage_matrix() if row["quantity"] == "association"]

    assert association_rows
    assert association_rows[0]["backend"] == "analytic_implicit"
    assert association_rows[0]["backend"] != "cppad"


def test_association_composition_derivative_includes_solved_site_fraction_response() -> None:
    mix, _species, _pressure, density, temperature, composition = _ionic_state()
    native_state = _core.NativeState(
        mix._native,
        temperature,
        composition.tolist(),
        0,
        False,
        0.0,
        True,
        density,
        False,
        0.0,
    )

    association_jacobian = np.asarray(
        native_state.composition_derivative_residual_helmholtz_result().dadx.assoc,
        dtype=float,
    )

    step = 1.0e-6
    central_perturbation = []
    for index in range(composition.size):
        plus = composition.copy()
        minus = composition.copy()
        plus[index] += step
        minus[index] -= step
        plus_value = _core.NativeState(
            mix._native,
            temperature,
            plus.tolist(),
            0,
            False,
            0.0,
            True,
            density,
            False,
            0.0,
        ).residual_helmholtz_result().assoc
        minus_value = _core.NativeState(
            mix._native,
            temperature,
            minus.tolist(),
            0,
            False,
            0.0,
            True,
            density,
            False,
            0.0,
        ).residual_helmholtz_result().assoc
        central_perturbation.append((plus_value - minus_value) / (2.0 * step))

    assert np.any(np.abs(association_jacobian) > 1.0e-8)
    np.testing.assert_allclose(association_jacobian, central_perturbation, rtol=2.0e-7, atol=2.0e-7)

    public_result = mix.state(T=temperature, x=composition, rho=density).ares_composition_derivative_result()
    assert public_result["backend"] == "analytic_implicit"
    assert public_result["backend_details"]["assoc"] == "analytic_implicit"


def test_direct_cppad_eos_contribution_recording_rejects_active_association() -> None:
    mix, _species, _pressure, density, temperature, composition = _ionic_state()
    args = create_struct(mix.parameters)

    if not _core._native_cppad_smoke()["cppad_compiled"]:
        return

    with pytest.raises(_core.NativeValueError, match="not_available"):
        _core._native_cppad_eos_contributions(temperature, density, composition.tolist(), args)
