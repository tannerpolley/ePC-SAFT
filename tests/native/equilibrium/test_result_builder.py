from __future__ import annotations

import numpy as np
import pytest

import epcsaft
from epcsaft import _core


def _neutral_binary_mixture() -> epcsaft.ePCSAFTMixture:
    params = {
        "m": np.asarray([1.0, 1.6069]),
        "s": np.asarray([3.7039, 3.5206]),
        "e": np.asarray([150.03, 191.42]),
        "k_ij": np.asarray([[0.0, 3.0e-4], [3.0e-4, 0.0]]),
    }
    return epcsaft.ePCSAFTMixture.from_params(params, species=["Methane", "Ethane"])


def test_neutral_two_phase_eos_result_builder_translates_accepted_candidate() -> None:
    mix = _neutral_binary_mixture()
    temperature = 300.0
    phase_amounts = [
        np.asarray([0.7, 0.3], dtype=float),
        np.asarray([0.1, 0.9], dtype=float),
    ]
    density = 120.0
    volumes = [float(phase.sum() / density) for phase in phase_amounts]
    feed_amounts = phase_amounts[0] + phase_amounts[1]
    target_pressure = mix.state(
        T=temperature,
        rho=density,
        x=phase_amounts[0] / phase_amounts[0].sum(),
        phase="liquid",
    ).pressure()

    payload = _core._native_neutral_two_phase_eos_result(
        mix._native,
        temperature,
        target_pressure,
        [phase.tolist() for phase in phase_amounts],
        volumes,
        feed_amounts.tolist(),
        1.0e-8,
        1.0e12,
        1.0e12,
        1.0e-3,
    )

    assert payload["accepted"] is True
    assert payload["backend"] == "native_equilibrium_nlp"
    assert payload["problem_kind"] == "neutral_two_phase_eos"
    assert payload["derivative_backend"] == "analytic_cppad"
    assert payload["split_detected"] is True
    assert payload["stable"] is False
    assert payload["rejection_reason"] == "accepted"
    assert payload["phase_labels"] == ["phase_0", "phase_1"]
    assert payload["phase_distance"] > 1.0e-3
    assert payload["material_balance_norm"] <= 1.0e-12
    assert payload["diagnostics"]["material_balance_norm"] <= 1.0e-12
    assert len(payload["phases"]) == 2

    total_amount = sum(float(phase.sum()) for phase in phase_amounts)
    for index, phase in enumerate(payload["phases"]):
        assert phase["label"] == f"phase_{index}"
        assert phase["temperature"] == pytest.approx(temperature)
        assert phase["pressure"] == pytest.approx(target_pressure)
        assert phase["density"] == pytest.approx(density)
        assert phase["phase_fraction"] == pytest.approx(float(phase_amounts[index].sum()) / total_amount)
        assert phase["composition"] == pytest.approx((phase_amounts[index] / phase_amounts[index].sum()).tolist())
        assert phase["amount_total"] == pytest.approx(float(phase_amounts[index].sum()))
        assert phase["volume"] == pytest.approx(volumes[index])
        assert len(phase["ln_fugacity_coefficient"]) == 2
        assert phase["fugacity_coefficient"] == pytest.approx(phase["ln_fugacity_coefficient"])
        assert np.all(np.isfinite(np.asarray(phase["ln_fugacity_coefficient"], dtype=float)))
        assert phase["diagnostics"]["volume"] == pytest.approx(volumes[index])


def test_neutral_two_phase_eos_result_builder_returns_rejection_without_phases() -> None:
    mix = _neutral_binary_mixture()
    temperature = 300.0
    phase_amounts = [
        np.asarray([0.4, 0.6], dtype=float),
        np.asarray([0.4, 0.6], dtype=float),
    ]
    density = 120.0
    volumes = [float(phase.sum() / density) for phase in phase_amounts]
    feed_amounts = phase_amounts[0] + phase_amounts[1]
    target_pressure = mix.state(
        T=temperature,
        rho=density,
        x=phase_amounts[0] / phase_amounts[0].sum(),
        phase="liquid",
    ).pressure()

    payload = _core._native_neutral_two_phase_eos_result(
        mix._native,
        temperature,
        target_pressure,
        [phase.tolist() for phase in phase_amounts],
        volumes,
        feed_amounts.tolist(),
        1.0e-8,
        1.0e-6,
        1.0e-6,
        1.0e-3,
    )

    assert payload["accepted"] is False
    assert payload["rejection_reason"] == "phase_distance"
    assert payload["split_detected"] is False
    assert payload["stable"] is False
    assert payload["phase_labels"] == []
    assert payload["phases"] == []
    assert payload["phase_distance"] == pytest.approx(0.0, abs=1.0e-14)
    assert payload["diagnostics"]["rejection_reason"] == "phase_distance"
