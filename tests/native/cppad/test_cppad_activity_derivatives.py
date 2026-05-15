from __future__ import annotations

from tests.helpers.native_cases import _ionic_state


def _state():
    mix, _species, _pressure, density, temperature, composition = _ionic_state()
    return mix.state(T=temperature, x=composition, rho=density)


def test_activity_composition_derivative_result_truthfully_unavailable() -> None:
    state = _state()

    result = state.activity_composition_derivative_result()

    assert result["supported"] is False
    assert result["backend"] == "not_available"
    assert result["shape"] == [state.x.size, state.x.size]
    assert "numerical_derivative" not in str(result).lower()


def test_activity_parameter_derivative_result_reports_not_available_without_ssmds_path() -> None:
    state = _state()

    result = state.activity_parameter_derivative_result()

    assert result["supported"] is False
    assert result["backend"] == "not_available"
    assert result["shape"] == [state.x.size, 0]
    assert "numerical_derivative" not in str(result).lower()
