from __future__ import annotations

import json

import numpy as np
import pytest

import epcsaft


def _co2_water_salt_mixture() -> epcsaft.ePCSAFTMixture:
    params = {
        "m": np.asarray([2.0729, 1.2047, 1.0, 1.0]),
        "s": np.asarray([2.7852, 2.7927, 2.8232, 2.7560]),
        "e": np.asarray([169.21, 353.95, 230.0, 170.0]),
        "z": np.asarray([0.0, 0.0, 1.0, -1.0]),
        "dielc": np.asarray([1.6, 78.09, 8.0, 8.0]),
        "d_born": np.asarray([0.0, 0.0, 3.445, 4.1]),
        "MW": np.asarray([44.0095e-3, 18.01528e-3, 22.989e-3, 35.45e-3]),
    }
    return epcsaft.ePCSAFTMixture.from_params(params, species=["CO2", "H2O", "Na+", "Cl-"])


def test_electrolyte_bubble_pressure_solves_neutral_vapor_over_ionic_liquid() -> None:
    mix = _co2_water_salt_mixture()
    x_liq = np.asarray([0.02, 0.979, 0.0005, 0.0005], dtype=float)

    result = mix.equilibrium(
        kind="electrolyte_bubble_pressure",
        T=313.15,
        x_liq=x_liq,
        volatile_species=["CO2", "H2O"],
        vapor_species=["CO2", "H2O"],
        nonvolatile_species=["Na+", "Cl-"],
        options=epcsaft.ElectrolyteBubbleOptions(initial_pressure=1.0e5, max_iterations=80),
    )

    assert isinstance(result, epcsaft.ElectrolyteBubbleResult)
    assert result.success is True
    assert result.P > 0.0
    np.testing.assert_allclose(result.x_liq, x_liq)
    assert set(result.y_vap) == {"CO2", "H2O"}
    assert sum(result.y_vap.values()) == pytest.approx(1.0, abs=1.0e-8)
    assert "Na+" not in result.y_vap
    assert "Cl-" not in result.y_vap
    assert result.charge_residual == pytest.approx(0.0, abs=1.0e-10)
    assert result.fugacity_residual_norm <= result.diagnostics["tolerance"]
    json.dumps(result.to_dict(), allow_nan=False)


@pytest.mark.parametrize(
    ("kwargs", "match"),
    [
        ({}, "x_liq"),
        ({"x_liq": [0.02, 0.979, 0.0004, 0.0006]}, "charge neutral"),
        ({"x_liq": [0.02, 0.979, 0.0005, 0.0005], "vapor_species": ["Na+"]}, "neutral vapor"),
        ({"x_liq": [0.02, 0.979, 0.0005, 0.0005], "volatile_species": ["missing"]}, "Unknown species"),
    ],
)
def test_electrolyte_bubble_pressure_rejects_invalid_public_inputs(kwargs, match) -> None:
    mix = _co2_water_salt_mixture()
    payload = {
        "kind": "electrolyte_bubble_pressure",
        "T": 313.15,
        "volatile_species": ["CO2", "H2O"],
        "vapor_species": ["CO2", "H2O"],
        "nonvolatile_species": ["Na+", "Cl-"],
    }
    payload.update(kwargs)

    with pytest.raises(epcsaft.InputError, match=match):
        mix.equilibrium(**payload)
