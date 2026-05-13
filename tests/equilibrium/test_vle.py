from __future__ import annotations

import json

import numpy as np
import pytest

import epcsaft
from epcsaft import ePCSAFTMixture


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


def test_ternary_hydrocarbon_basis_tp_flash_closes_material_and_fugacity_balance() -> None:
    mix = _hydrocarbon_basis_mixture()
    feed = np.asarray([0.1, 0.3, 0.6])

    result = mix.equilibrium(kind="tp_flash", T=220.0, P=1.0e5, z=feed)

    assert result.split_detected is True
    assert result.stable is False
    assert result.phase_labels == ["liq", "vap"]
    liquid, vapor = result.phases
    assert liquid.phase_fraction == pytest.approx(0.0717673735624358, abs=2.0e-6)
    assert vapor.phase_fraction == pytest.approx(0.9282326264375642, abs=2.0e-6)
    np.testing.assert_allclose(
        liquid.composition,
        np.asarray([0.0012963789214619132, 0.06534426759935694, 0.9333593534791812]),
        atol=2.0e-6,
    )
    np.testing.assert_allclose(
        vapor.composition,
        np.asarray([0.10763138403472723, 0.3181426779517501, 0.5742259380135226]),
        atol=2.0e-6,
    )
    np.testing.assert_allclose(liquid.composition.sum(), 1.0)
    np.testing.assert_allclose(vapor.composition.sum(), 1.0)
    assert np.all(liquid.composition > 0.0)
    assert np.all(vapor.composition > 0.0)

    reconstructed = liquid.phase_fraction * liquid.composition + vapor.phase_fraction * vapor.composition
    np.testing.assert_allclose(reconstructed, feed, atol=1.0e-10)
    assert result.diagnostics["material_balance_error"] < 1.0e-10
    assert result.diagnostics["fugacity_residual_norm"] < 1.0e-6
    assert result.diagnostics["equilibrium_route"] == "neutral_vle"
    assert result.diagnostics["route_reason"] == "requested vapor-liquid path"
    assert result.diagnostics["stability_analysis"] == "neutral_tpd"
    assert result.diagnostics["neutral_fast_path"] is True
    assert result.diagnostics["neutral_fallback_used"] is False
    assert result.diagnostics["neutral_fallback_reason"] == ""
    assert isinstance(result.diagnostics["stability_stable"], bool)
    assert np.isfinite(result.diagnostics["min_tpd"])
    assert result.diagnostics["unstable_trial_count"] >= 1

    fugacity_residual = (
        np.log(vapor.composition)
        + vapor.ln_fugacity_coefficient
        - np.log(liquid.composition)
        - liquid.ln_fugacity_coefficient
    )
    np.testing.assert_allclose(fugacity_residual, np.zeros_like(feed), atol=1.0e-6)


def test_tp_flash_reports_no_split_when_rachford_rice_has_no_bracket() -> None:
    mix = _hydrocarbon_basis_mixture()

    result = mix.equilibrium(kind="tp_flash", T=300.0, P=1.0e5, z=[0.1, 0.3, 0.6])

    assert result.split_detected is False
    assert result.stable is False
    assert len(result.phases) == 1
    assert result.phases[0].label in {"liq", "vap"}
    assert "no two-phase Rachford-Rice bracket" in result.diagnostics["message"]
    assert result.diagnostics["point_solver_split_detected"] is False
    assert "no two-phase Rachford-Rice bracket" in result.diagnostics["point_solver_message"]
    assert result.diagnostics["stability_analysis"] == "neutral_tpd"
    assert result.diagnostics["stability_stable"] is False
    assert np.isfinite(result.diagnostics["min_tpd"])


def test_tp_flash_can_skip_stability_precheck_for_debug_workflows() -> None:
    mix = _hydrocarbon_basis_mixture()

    result = mix.equilibrium(
        kind="tp_flash",
        T=300.0,
        P=1.0e5,
        z=[0.1, 0.3, 0.6],
        options=epcsaft.EquilibriumOptions(stability_precheck=False),
    )

    assert result.split_detected is False
    assert result.stable is False
    assert result.diagnostics["stability_analysis"] == "not_run"
    assert result.diagnostics["stability_checked"] is False
    assert result.diagnostics["stability_stable"] is None
    assert "stability precheck skipped" in result.diagnostics["stability_message"]
    assert result.diagnostics["point_solver_split_detected"] is False
    assert "min_tpd" not in result.diagnostics


def test_tp_flash_phase_diagnostics_are_json_serializable_when_requested() -> None:
    mix = _hydrocarbon_basis_mixture()

    result = mix.equilibrium(
        kind="tp_flash",
        T=220.0,
        P=1.0e5,
        z=[0.1, 0.3, 0.6],
        options=epcsaft.EquilibriumOptions(include_phase_diagnostics=True),
    )

    payload = result.to_dict()
    json.dumps(payload, allow_nan=False)
    for phase, phase_payload in zip(result.phases, payload["phases"]):
        assert phase.diagnostics is not None
        assert "phase" in phase_payload["diagnostics"]
        assert "density" in phase_payload["diagnostics"]
        assert "fugacity_coefficient_terms" in phase_payload["diagnostics"]
