from __future__ import annotations

import numpy as np
import pytest

import epcsaft
from epcsaft import _core, ePCSAFTMixture


def _hydrocarbon_basis_mixture() -> ePCSAFTMixture:
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
    return ePCSAFTMixture.from_params(params, species=["Methane", "Ethane", "Propane"])


def _assert_tp_flash_route_pending(excinfo: pytest.ExceptionInfo[epcsaft.InputError]) -> None:
    message = str(excinfo.value)
    assert "tp_flash requires a native Ipopt equilibrium NLP route" in message
    assert "No package-owned alternate TP flash solver is available" in message


def test_ternary_hydrocarbon_basis_tp_flash_requires_native_ipopt_route() -> None:
    mix = _hydrocarbon_basis_mixture()
    feed = np.asarray([0.1, 0.3, 0.6])

    with pytest.raises(epcsaft.InputError) as excinfo:
        mix.equilibrium(kind="tp_flash", T=220.0, P=1.0e5, z=feed)

    _assert_tp_flash_route_pending(excinfo)


def test_tp_flash_builds_one_native_route_request_before_ipopt_gate(monkeypatch: pytest.MonkeyPatch) -> None:
    mix = _hydrocarbon_basis_mixture()
    feed = np.asarray([0.1, 0.3, 0.6])
    calls: list[dict[str, object]] = []

    def fake_route(
        _native,
        temperature,
        pressure,
        feed_amounts,
        max_iterations,
        tolerance,
        material_tolerance,
        pressure_tolerance,
        chemical_potential_tolerance,
        phase_distance_tolerance,
    ):
        calls.append(
            {
                "temperature": temperature,
                "pressure": pressure,
                "feed_amounts": feed_amounts,
                "max_iterations": max_iterations,
                "tolerance": tolerance,
                "material_tolerance": material_tolerance,
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
            "status": "requires_ipopt_build",
            "postsolve": {"accepted": False},
        }

    monkeypatch.setattr(_core, "_native_neutral_tp_flash_eos_route_result", fake_route)

    with pytest.raises(epcsaft.InputError) as excinfo:
        mix.equilibrium(
            kind="tp_flash",
            T=220.0,
            P=1.0e5,
            z=feed,
            options=epcsaft.EquilibriumOptions(max_iterations=17, tolerance=2.0e-7),
        )

    _assert_tp_flash_route_pending(excinfo)
    assert len(calls) == 1
    call = calls[0]
    assert call["feed_amounts"] == pytest.approx(feed.tolist())
    assert call["max_iterations"] == 17
    assert call["tolerance"] == pytest.approx(2.0e-7)


def test_tp_flash_converts_accepted_native_route_payload(monkeypatch: pytest.MonkeyPatch) -> None:
    mix = _hydrocarbon_basis_mixture()
    feed = np.asarray([0.1, 0.3, 0.6])
    route_amounts = [[0.07, 0.12, 0.21], [0.03, 0.18, 0.39]]
    route_volumes = [0.01, 0.02]

    def fake_route(
        _native,
        temperature,
        pressure,
        feed_amounts,
        max_iterations,
        tolerance,
        material_tolerance,
        pressure_tolerance,
        chemical_potential_tolerance,
        phase_distance_tolerance,
    ):
        assert temperature == pytest.approx(220.0)
        assert pressure == pytest.approx(1.0e5)
        assert feed_amounts == pytest.approx(feed.tolist())
        assert max_iterations > 0
        assert tolerance > 0.0
        assert material_tolerance > 0.0
        assert pressure_tolerance > 0.0
        assert chemical_potential_tolerance > 0.0
        assert phase_distance_tolerance > 0.0
        return {
            "backend": "ipopt",
            "compiled": True,
            "ran": True,
            "accepted": True,
            "status": "accepted",
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
        assert temperature == pytest.approx(220.0)
        assert pressure == pytest.approx(1.0e5)
        assert phase_amounts == route_amounts
        assert phase_volumes == route_volumes
        assert feed_amounts == pytest.approx(feed.tolist())
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
                    "composition": [0.175, 0.3, 0.525],
                    "density": 40.0,
                    "temperature": 220.0,
                    "pressure": 1.0e5,
                    "phase_fraction": 0.4,
                    "ln_fugacity_coefficient": [0.0, 0.1, 0.2],
                    "fugacity_coefficient": [1.0, float(np.exp(0.1)), float(np.exp(0.2))],
                },
                {
                    "label": "phase_1",
                    "composition": [0.05, 0.3, 0.65],
                    "density": 30.0,
                    "temperature": 220.0,
                    "pressure": 1.0e5,
                    "phase_fraction": 0.6,
                    "ln_fugacity_coefficient": [0.2, 0.1, 0.0],
                    "fugacity_coefficient": [float(np.exp(0.2)), float(np.exp(0.1)), 1.0],
                },
            ],
        }

    monkeypatch.setattr(_core, "_native_neutral_tp_flash_eos_route_result", fake_route)
    monkeypatch.setattr(_core, "_native_neutral_two_phase_eos_result", fake_result)

    result = mix.equilibrium(kind="tp_flash", T=220.0, P=1.0e5, z=feed)

    assert result.backend == "native_equilibrium_nlp"
    assert result.problem_kind == "neutral_two_phase_eos"
    assert result.split_detected is True
    assert [phase.label for phase in result.phases] == ["phase_0", "phase_1"]
    assert result.phases[0].fugacity_coefficient == pytest.approx(np.exp(result.phases[0].ln_fugacity_coefficient))


def test_tp_flash_accepts_stability_precheck_option_before_route_gate() -> None:
    mix = _hydrocarbon_basis_mixture()

    with pytest.raises(epcsaft.InputError) as excinfo:
        mix.equilibrium(
            kind="tp_flash",
            T=300.0,
            P=1.0e5,
            z=[0.1, 0.3, 0.6],
            options=epcsaft.EquilibriumOptions(stability_precheck=False),
        )

    _assert_tp_flash_route_pending(excinfo)


def test_tp_flash_phase_diagnostics_option_requires_native_ipopt_route() -> None:
    mix = _hydrocarbon_basis_mixture()

    with pytest.raises(epcsaft.InputError) as excinfo:
        mix.equilibrium(
            kind="tp_flash",
            T=220.0,
            P=1.0e5,
            z=[0.1, 0.3, 0.6],
            options=epcsaft.EquilibriumOptions(include_phase_diagnostics=True),
        )

    _assert_tp_flash_route_pending(excinfo)
