from __future__ import annotations

from dataclasses import fields

import numpy as np
import pytest

import epcsaft
from epcsaft import InputError, _core


def _salt_mixture() -> epcsaft.ePCSAFTMixture:
    x = np.asarray([0.98, 0.01, 0.01], dtype=float)
    return epcsaft.ePCSAFTMixture.from_dataset("2026_Khudaida", ["H2O", "Na+", "Cl-"], x, 298.15)


def test_electrolyte_bubble_pressure_builds_native_route_before_ipopt_gate(monkeypatch) -> None:
    assert {field.name for field in fields(epcsaft.ElectrolyteBubbleOptions)} == {
        "max_iterations",
        "tolerance",
        "min_composition",
        "charge_tolerance",
        "hessian_mode",
        "ipopt_iteration_history_limit",
        "continuation_state",
    }
    mix = _salt_mixture()
    calls: list[dict[str, object]] = []

    def fake_route(
        _native,
        temperature,
        liquid_composition,
        max_iterations,
        tolerance,
        hessian_mode,
        iteration_history_limit,
        phase_total_tolerance,
        pressure_tolerance,
        charge_tolerance,
        chemical_potential_tolerance,
        phase_distance_tolerance,
        continuation_state,
    ):
        calls.append(
            {
                "temperature": temperature,
                "liquid_composition": liquid_composition,
                "max_iterations": max_iterations,
                "tolerance": tolerance,
                "continuation_state": continuation_state,
                "phase_total_tolerance": phase_total_tolerance,
                "pressure_tolerance": pressure_tolerance,
                "charge_tolerance": charge_tolerance,
                "chemical_potential_tolerance": chemical_potential_tolerance,
                "phase_distance_tolerance": phase_distance_tolerance,
            }
        )
        return {
            "backend": "ipopt",
            "compiled": False,
            "ran": False,
            "accepted": False,
            "status": "ipopt_dependency_required",
            "postsolve": {"accepted": False},
        }

    monkeypatch.setattr(_core, "_native_electrolyte_bubble_p_eos_route_result", fake_route)

    with pytest.raises(InputError, match="native Ipopt equilibrium route builder"):
        mix.equilibrium(
            kind="electrolyte_bubble_pressure",
            T=298.15,
            x_liq=[0.98, 0.01, 0.01],
            vapor_species=["H2O"],
            volatile_species=["H2O"],
            nonvolatile_species=["Na+", "Cl-"],
            backend="native",
        )

    assert len(calls) == 1
    call = calls[0]
    assert call["temperature"] == pytest.approx(298.15)
    assert call["liquid_composition"] == pytest.approx([0.98, 0.01, 0.01])
    assert call["max_iterations"] == 80
    assert call["tolerance"] == pytest.approx(1.0e-6)
    assert call["phase_total_tolerance"] == pytest.approx(1.0e-6)
    assert call["pressure_tolerance"] == pytest.approx(0.1)
    assert call["charge_tolerance"] == pytest.approx(1.0e-8)
    assert call["chemical_potential_tolerance"] == pytest.approx(1.0e-6)
    assert call["phase_distance_tolerance"] > 0.0


def test_electrolyte_bubble_pressure_rejects_python_backend_alias() -> None:
    mix = _salt_mixture()

    with pytest.raises(InputError, match="backend must be None or 'native'"):
        mix.equilibrium(
            kind="electrolyte_bubble_pressure",
            T=298.15,
            x_liq=[0.98, 0.01, 0.01],
            vapor_species=["H2O"],
            backend="python",
        )
