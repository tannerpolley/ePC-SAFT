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

    assert association_jacobian.shape == (composition.size,)
    assert np.all(np.isfinite(association_jacobian))
    assert np.any(np.abs(association_jacobian) > 1.0e-8)

    public_result = mix.state(T=temperature, x=composition, rho=density).ares_composition_derivative_result()
    assert public_result["backend"] == "analytic_implicit"
    assert public_result["backend_details"]["assoc"] == "analytic_implicit"
    public_jacobian = np.asarray(public_result["jacobian"], dtype=float)
    assert public_jacobian.shape == (1, composition.size)
    assert np.all(np.isfinite(public_jacobian))


def test_direct_cppad_eos_contribution_recording_rejects_active_association() -> None:
    mix, _species, _pressure, density, temperature, composition = _ionic_state()
    args = create_struct(mix.parameters)

    if not _core._native_cppad_smoke()["cppad_compiled"]:
        return

    with pytest.raises(_core.NativeValueError, match="unsupported"):
        _core._native_cppad_eos_contributions(temperature, density, composition.tolist(), args)
