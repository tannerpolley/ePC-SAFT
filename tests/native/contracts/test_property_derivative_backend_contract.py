from __future__ import annotations

import epcsaft
import pytest
from epcsaft._types import InputError
from tests.helpers.native_cases import _neutral_state


def _state():
    mix, _species, _pressure, density, temperature, composition = _neutral_state()
    return mix.state(T=temperature, x=composition, rho=density)


def test_capabilities_expose_property_derivative_result_apis() -> None:
    capabilities = epcsaft.capabilities()

    payload = capabilities["derivatives"]["property_derivative_result_apis"]

    assert payload["available"] is True
    assert payload["numerical_derivative_backend_available"] is False
    assert "outputs" in payload["result_shape"]
    assert "variables" in payload["result_shape"]
    assert "derivative_backend" in payload["result_shape"]
    assert "pressure_density_derivative_result" in payload["state_methods"]
    assert "relative_permittivity_parameter_derivative_result" in payload["state_methods"]
    assert "numerical_derivative" not in str(payload["backend_labels"])
    assert "not_available" not in str(payload["backend_labels"])


def test_public_derivative_result_methods_share_required_shape() -> None:
    state = _state()

    for method_name in (
        "pressure_density_derivative_result",
        "ares_composition_derivative_result",
    ):
        result = getattr(state, method_name)()
        assert set(
            ("supported", "backend", "derivative_backend", "message", "value", "jacobian", "outputs", "variables", "shape")
        ).issubset(result)
        assert isinstance(result["shape"], list)
        assert len(result["shape"]) == 2
        assert isinstance(result["outputs"], list)
        assert isinstance(result["variables"], list)
        assert "numerical_derivative" not in str(result).lower()


def test_unsupported_property_derivative_methods_raise() -> None:
    state = _state()

    for method_name in (
        "pressure_composition_derivative_result",
        "pressure_parameter_derivative_result",
        "density_pressure_derivative_result",
        "chemical_potential_composition_derivative_result",
        "chemical_potential_parameter_derivative_result",
        "ln_fugacity_composition_derivative_result",
        "ln_fugacity_parameter_derivative_result",
    ):
        with pytest.raises(InputError, match="unsupported"):
            getattr(state, method_name)()
