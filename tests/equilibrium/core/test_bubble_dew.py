from __future__ import annotations

import numpy as np
import pytest

import epcsaft
from epcsaft import _core


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


def test_bubble_p_builds_one_native_route_request_before_ipopt_gate(monkeypatch: pytest.MonkeyPatch) -> None:
    mix = _hydrocarbon_mixture()
    calls: list[dict[str, object]] = []

    def fake_route(
        _native,
        temperature,
        liquid_composition,
        max_iterations,
        tolerance,
        phase_total_tolerance,
        pressure_tolerance,
        chemical_potential_tolerance,
        phase_distance_tolerance,
    ):
        calls.append(
            {
                "temperature": temperature,
                "liquid_composition": liquid_composition,
                "max_iterations": max_iterations,
                "tolerance": tolerance,
                "phase_total_tolerance": phase_total_tolerance,
                "pressure_tolerance": pressure_tolerance,
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

    monkeypatch.setattr(_core, "_native_neutral_bubble_p_eos_route_result", fake_route)

    with pytest.raises(epcsaft.InputError, match=r"bubble_p requires a native Ipopt equilibrium NLP route"):
        mix.bubble_p(
            T=220.0,
            x=[0.2, 0.3, 0.5],
            options=epcsaft.EquilibriumOptions(max_iterations=17, tolerance=4.0e-8),
        )

    assert len(calls) == 1
    call = calls[0]
    assert call["temperature"] == pytest.approx(220.0)
    assert call["liquid_composition"] == pytest.approx([0.2, 0.3, 0.5])
    assert call["max_iterations"] == 17
    assert call["tolerance"] == pytest.approx(4.0e-8)
    assert call["phase_total_tolerance"] > 0.0
    assert call["pressure_tolerance"] == pytest.approx(4.0e-3)
    assert call["chemical_potential_tolerance"] > 0.0
    assert call["phase_distance_tolerance"] > 0.0


def test_dew_p_converts_accepted_native_route_payload(monkeypatch: pytest.MonkeyPatch) -> None:
    mix = _hydrocarbon_mixture()
    vapor_composition = [0.1, 0.3, 0.6]
    route_amounts = [[0.5, 0.35, 0.15], [0.1, 0.3, 0.6]]
    route_volumes = [0.001, 0.02]
    solved_pressure = 2.5e5

    def fake_route(
        _native,
        temperature,
        received_vapor_composition,
        max_iterations,
        tolerance,
        phase_total_tolerance,
        pressure_tolerance,
        chemical_potential_tolerance,
        phase_distance_tolerance,
    ):
        assert temperature == pytest.approx(240.0)
        assert received_vapor_composition == pytest.approx(vapor_composition)
        assert max_iterations > 0
        assert tolerance > 0.0
        assert phase_total_tolerance > 0.0
        assert pressure_tolerance > 0.0
        assert chemical_potential_tolerance > 0.0
        assert phase_distance_tolerance > 0.0
        return {
            "backend": "ipopt",
            "compiled": True,
            "ran": True,
            "accepted": True,
            "status": "accepted",
            "variables": [*route_amounts[0], route_volumes[0], *route_amounts[1], route_volumes[1], solved_pressure],
            "phase_amounts": route_amounts,
            "phase_volumes": route_volumes,
            "postsolve": {"accepted": True},
        }

    def fake_result(
        _native,
        temperature,
        pressure,
        phase_amounts,
        phase_volumes,
        feed_amounts,
        material_tolerance,
        pressure_tolerance,
        chemical_potential_tolerance,
        phase_distance_tolerance,
    ):
        assert temperature == pytest.approx(240.0)
        assert pressure == pytest.approx(solved_pressure)
        assert phase_amounts == route_amounts
        assert phase_volumes == route_volumes
        assert feed_amounts == pytest.approx([0.6, 0.65, 0.75])
        assert material_tolerance > 0.0
        assert pressure_tolerance > 0.0
        assert chemical_potential_tolerance > 0.0
        assert phase_distance_tolerance > 0.0
        return {
            "accepted": True,
            "backend": "native_equilibrium_nlp",
            "problem_kind": "neutral_two_phase_eos",
            "stable": False,
            "split_detected": True,
            "phases": [
                {
                    "label": "phase_0",
                    "composition": [0.5, 0.35, 0.15],
                    "density": 750.0,
                    "temperature": 240.0,
                    "pressure": solved_pressure,
                    "phase_fraction": 0.5,
                    "ln_fugacity_coefficient": [0.2, 0.1, 0.0],
                    "fugacity_coefficient": [float(np.exp(0.2)), float(np.exp(0.1)), 1.0],
                },
                {
                    "label": "phase_1",
                    "composition": vapor_composition,
                    "density": 20.0,
                    "temperature": 240.0,
                    "pressure": solved_pressure,
                    "phase_fraction": 0.5,
                    "ln_fugacity_coefficient": [0.0, 0.1, 0.2],
                    "fugacity_coefficient": [1.0, float(np.exp(0.1)), float(np.exp(0.2))],
                },
            ],
        }

    monkeypatch.setattr(_core, "_native_neutral_dew_p_eos_route_result", fake_route)
    monkeypatch.setattr(_core, "_native_neutral_two_phase_eos_result", fake_result)

    result = mix.dew_p(T=240.0, y=vapor_composition)

    assert result.backend == "native_equilibrium_nlp"
    assert result.problem_kind == "neutral_dew_p"
    assert [phase.label for phase in result.phases] == ["liq", "vap"]
    assert result.phases[1].composition == pytest.approx(vapor_composition)
    assert result.phases[0].pressure == pytest.approx(solved_pressure)
