from __future__ import annotations

import json

import numpy as np
import pytest

import epcsaft


def _hydrocarbon_mixture() -> epcsaft.ePCSAFTMixture:
    params = {
        "m": np.asarray([1.0, 1.6069, 2.0020]),
        "s": np.asarray([3.7039, 3.5206, 3.6184]),
        "e": np.asarray([150.03, 191.42, 208.11]),
        "k_ij": np.asarray(
            [
                [0.0, 3.0e-4, 1.15e-2],
                [3.0e-4, 0.0, 5.10e-3],
                [1.15e-2, 5.10e-3, 0.0],
            ]
        ),
    }
    return epcsaft.ePCSAFTMixture.from_params(params, species=["Methane", "Ethane", "Propane"])


def test_neutral_bubble_and_dew_pressure_match_tp_flash_endpoint() -> None:
    mix = _hydrocarbon_mixture()
    flash = mix.flash_tp(T=220.0, P=1.0e5, z=[0.1, 0.3, 0.6])
    liquid, vapor = flash.phases

    bubble = mix.bubble_p(T=220.0, x=liquid.composition)
    dew = mix.dew_p(T=220.0, y=vapor.composition)

    assert bubble.problem_kind == "bubble_p"
    assert dew.problem_kind == "dew_p"
    assert bubble.phases[0].label == "liq"
    assert bubble.phases[1].label == "vap"
    assert dew.phases[0].label == "liq"
    assert dew.phases[1].label == "vap"
    assert bubble.phases[0].phase_fraction == pytest.approx(1.0)
    assert bubble.phases[1].phase_fraction == pytest.approx(0.0)
    assert dew.phases[0].phase_fraction == pytest.approx(0.0)
    assert dew.phases[1].phase_fraction == pytest.approx(1.0)
    assert bubble.phases[1].composition == pytest.approx(vapor.composition, abs=2.0e-5)
    assert dew.phases[0].composition == pytest.approx(liquid.composition, abs=2.0e-5)
    assert bubble.phases[0].pressure == pytest.approx(1.0e5, rel=2.0e-5)
    assert dew.phases[1].pressure == pytest.approx(1.0e5, rel=2.0e-5)
    assert bubble.diagnostics["fugacity_residual_norm"] < 2.0e-5
    assert dew.diagnostics["fugacity_residual_norm"] < 2.0e-5
    assert bubble.diagnostics["equilibrium_route"] == "neutral_vle"
    assert dew.diagnostics["equilibrium_route"] == "neutral_vle"
    assert bubble.diagnostics["neutral_fast_path"] is True
    assert dew.diagnostics["neutral_fast_path"] is True
    assert bubble.diagnostics["neutral_fallback_used"] is False
    assert dew.diagnostics["neutral_fallback_used"] is False
    assert bubble.diagnostics["neutral_fallback_reason"] == ""
    assert dew.diagnostics["neutral_fallback_reason"] == ""
    assert bubble.diagnostics["partial_pressures"]["Methane"] == pytest.approx(
        bubble.phases[1].composition[0] * bubble.phases[1].pressure
    )
    assert dew.diagnostics["partial_pressures"]["Methane"] == pytest.approx(
        dew.phases[1].composition[0] * dew.phases[1].pressure
    )
    assert bubble.diagnostics["volatile_partial_pressure_basis"] == "liquid_fugacity_equilibrium"
    assert dew.diagnostics["partial_pressure_route"] == "vapor_composition_times_total_pressure"
    json.dumps(bubble.to_dict(), allow_nan=False)
    json.dumps(dew.to_dict(), allow_nan=False)


def test_neutral_bubble_and_dew_temperature_match_tp_flash_endpoint() -> None:
    mix = _hydrocarbon_mixture()
    flash = mix.flash_tp(T=220.0, P=1.0e5, z=[0.1, 0.3, 0.6])
    liquid, vapor = flash.phases

    bubble = mix.bubble_t(P=1.0e5, x=liquid.composition)
    dew = mix.dew_t(P=1.0e5, y=vapor.composition)

    assert bubble.problem_kind == "bubble_t"
    assert dew.problem_kind == "dew_t"
    assert bubble.phases[0].temperature == pytest.approx(220.0, abs=2.0e-4)
    assert dew.phases[1].temperature == pytest.approx(220.0, abs=2.0e-4)
    assert bubble.phases[1].composition == pytest.approx(vapor.composition, abs=2.0e-5)
    assert dew.phases[0].composition == pytest.approx(liquid.composition, abs=2.0e-5)
    assert bubble.diagnostics["scalar_residual"] == pytest.approx(0.0, abs=2.0e-5)
    assert dew.diagnostics["scalar_residual"] == pytest.approx(0.0, abs=2.0e-5)
    assert bubble.diagnostics["equilibrium_route"] == "neutral_vle"
    assert dew.diagnostics["equilibrium_route"] == "neutral_vle"
