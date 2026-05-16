from __future__ import annotations

import numpy as np

from tests.helpers.native_cases import _neutral_state


def test_cppad_pressure_derivative_api_underpins_bubble_policy() -> None:
    mix, _species, _pressure, density, temperature, composition = _neutral_state()
    state = mix.state(T=temperature, x=composition, rho=density)

    pressure_density = state.pressure_density_derivative_result()

    assert pressure_density["backend"] in {"cppad", "not_available"}
    assert pressure_density["shape"] == [1, 1]
    assert "numerical_derivative" not in str(pressure_density).lower()


def test_neutral_bubble_reports_no_numerical_derivative_root_derivative() -> None:
    mix, _species, _pressure, _density, _temperature, composition = _neutral_state()
    flash = mix.flash_tp(T=220.0, P=1.0e5, z=composition)
    liquid = flash.phases[0]

    bubble = mix.bubble_p(T=220.0, x=liquid.composition)

    diagnostics = bubble.diagnostics
    assert bubble.problem_kind == "bubble_p"
    assert diagnostics["thermodynamic_backend"] == "epcsaft_state_fugacity_activity_property_api"
    assert diagnostics["derivative_backend"] == "not_available"
    assert diagnostics["derivative_backend_by_block"]["bubble_pressure_root"] == "not_available"
    assert diagnostics["jacobian_fallback_used"] is False
    assert np.isfinite(diagnostics["residual_norm"])
    assert "numerical_derivative" not in str(diagnostics).lower()
