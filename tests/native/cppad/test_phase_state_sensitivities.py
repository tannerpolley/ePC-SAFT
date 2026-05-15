from __future__ import annotations

import numpy as np

from epcsaft import _core
from epcsaft.epcsaft import create_struct
from tests.helpers.native_cases import _ionic_state, _neutral_state


def _neutral_pressure_state():
    mix, _species, pressure, _density, temperature, composition = _neutral_state()
    return mix, pressure, temperature, composition


def _phase_state_sensitivity(mix, temperature, pressure, composition):
    return _core._native_phase_state_ln_fugacity_composition_sensitivity(
        temperature,
        pressure,
        composition.tolist(),
        0,
        create_struct(mix.parameters),
    )


def _assert_projected_phase_state_sensitivity(raw, mix, temperature, pressure, composition, *, atol):
    shape = raw["shape"]
    jacobian = np.asarray(raw["jacobian_row_major"], dtype=float).reshape(shape)
    weighted_columns = jacobian @ composition
    step = 1.0e-6

    for column in range(composition.size):
        direction = -composition.copy()
        direction[column] += 1.0
        plus = composition + step * direction
        minus = composition - step * direction
        assert np.all(plus > 0.0)
        assert np.all(minus > 0.0)
        np.testing.assert_allclose(plus.sum(), 1.0, rtol=0.0, atol=1.0e-12)
        np.testing.assert_allclose(minus.sum(), 1.0, rtol=0.0, atol=1.0e-12)

        plus_raw = _phase_state_sensitivity(mix, temperature, pressure, plus)
        minus_raw = _phase_state_sensitivity(mix, temperature, pressure, minus)
        perturbation_derivative = (
            np.asarray(plus_raw["ln_fugacity"], dtype=float)
            - np.asarray(minus_raw["ln_fugacity"], dtype=float)
        ) / (2.0 * step)
        projected_column = jacobian[:, column] - weighted_columns
        np.testing.assert_allclose(projected_column, perturbation_derivative, rtol=2.0e-5, atol=atol)


def test_native_phase_state_sensitivity_uses_implicit_density_chain_rule() -> None:
    mix, pressure, temperature, composition = _neutral_pressure_state()
    raw = _phase_state_sensitivity(mix, temperature, pressure, composition)

    if not _core._native_cppad_smoke()["cppad_compiled"]:
        assert raw["supported"] is False
        return

    assert raw["supported"] is True
    assert raw["backend"] == "cppad_implicit"
    assert raw["density_backend"] == "implicit_density_root"
    assert raw["shape"] == (composition.size, composition.size)

    state = _core.NativeState(
        mix._native,
        temperature,
        composition.tolist(),
        0,
        True,
        pressure,
        False,
        0.0,
        False,
        0.0,
    )
    np.testing.assert_allclose(raw["ln_fugacity"], state.ln_fugacity_coefficient(), rtol=0.0, atol=1.0e-12)

    dpdrho = float(raw["pressure_density_derivative"])
    drhodx = np.asarray(raw["density_composition_derivative"], dtype=float)
    dpdx_fixed = np.asarray(raw["pressure_composition_fixed_density_derivative"], dtype=float)
    np.testing.assert_allclose(dpdx_fixed + dpdrho * drhodx, 0.0, rtol=1.0e-12, atol=1.0e-6)

    shape = raw["shape"]
    fixed_jacobian = np.asarray(raw["fixed_density_jacobian_row_major"], dtype=float).reshape(shape)
    total_jacobian = np.asarray(raw["jacobian_row_major"], dtype=float).reshape(shape)
    dlnphi_drho = np.asarray(raw["ln_fugacity_density_derivative"], dtype=float)
    np.testing.assert_allclose(
        total_jacobian,
        fixed_jacobian + np.outer(dlnphi_drho, drhodx),
        rtol=1.0e-12,
        atol=1.0e-12,
    )
    _assert_projected_phase_state_sensitivity(raw, mix, temperature, pressure, composition, atol=5.0e-6)


def test_public_pressure_state_ln_fugacity_composition_derivative_is_supported() -> None:
    mix, pressure, temperature, composition = _neutral_pressure_state()
    state = mix.state(T=temperature, P=pressure, x=composition)

    result = state.ln_fugacity_composition_derivative_result()

    if not _core._native_cppad_smoke()["cppad_compiled"]:
        assert result["supported"] is False
        return

    assert result["supported"] is True
    assert result["backend"] == "cppad_implicit"
    assert result["density_backend"] == "implicit_density_root"
    assert result["shape"] == [composition.size, composition.size]
    assert np.asarray(result["jacobian"], dtype=float).shape == (composition.size, composition.size)
    assert "numerical_derivative" not in str(result).lower()


def test_phase_state_sensitivity_supports_active_association_implicit_response() -> None:
    mix, _species, pressure, _density, temperature, composition = _ionic_state()
    raw = _phase_state_sensitivity(mix, temperature, pressure, composition)

    if not _core._native_cppad_smoke()["cppad_compiled"]:
        assert raw["supported"] is False
        return

    assert raw["supported"] is True
    assert raw["backend"] == "cppad_implicit"
    assert raw["density_backend"] == "implicit_density_root"
    assert raw["shape"] == (composition.size, composition.size)

    dpdrho = float(raw["pressure_density_derivative"])
    drhodx = np.asarray(raw["density_composition_derivative"], dtype=float)
    dpdx_fixed = np.asarray(raw["pressure_composition_fixed_density_derivative"], dtype=float)
    np.testing.assert_allclose(dpdx_fixed + dpdrho * drhodx, 0.0, rtol=1.0e-12, atol=1.0e-6)

    shape = raw["shape"]
    fixed_jacobian = np.asarray(raw["fixed_density_jacobian_row_major"], dtype=float).reshape(shape)
    total_jacobian = np.asarray(raw["jacobian_row_major"], dtype=float).reshape(shape)
    dlnphi_drho = np.asarray(raw["ln_fugacity_density_derivative"], dtype=float)
    np.testing.assert_allclose(
        total_jacobian,
        fixed_jacobian + np.outer(dlnphi_drho, drhodx),
        rtol=1.0e-12,
        atol=1.0e-12,
    )
    _assert_projected_phase_state_sensitivity(raw, mix, temperature, pressure, composition, atol=2.0e-4)
