from __future__ import annotations

import json

import numpy as np
import pytest

import epcsaft
from epcsaft import ePCSAFTMixture


def _hydrocarbon_mixture() -> ePCSAFTMixture:
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


def _ionic_mixture() -> ePCSAFTMixture:
    params = {
        "m": np.asarray([1.2047, 1.0, 1.0]),
        "s": np.asarray([2.7927, 2.8232, 2.7560]),
        "e": np.asarray([353.95, 230.0, 170.0]),
        "z": np.asarray([0.0, 1.0, -1.0]),
        "dielc": np.asarray([78.09, 8.0, 8.0]),
    }
    return ePCSAFTMixture.from_params(params, species=["water", "Na+", "Cl-"])


def _assert_json_like(value):
    if isinstance(value, dict):
        for item in value.values():
            _assert_json_like(item)
    elif isinstance(value, list):
        for item in value:
            _assert_json_like(item)
    else:
        assert not isinstance(value, np.ndarray)


def test_equilibrium_public_exports_are_available() -> None:
    assert hasattr(epcsaft, "EquilibriumOptions")
    assert hasattr(epcsaft, "EquilibriumPhase")
    assert hasattr(epcsaft, "EquilibriumResult")
    assert hasattr(epcsaft, "StabilityTrial")
    assert hasattr(epcsaft, "StabilityResult")
    assert hasattr(epcsaft, "bubble_p")
    assert hasattr(epcsaft, "bubble_t")
    assert hasattr(epcsaft, "dew_p")
    assert hasattr(epcsaft, "dew_t")


def test_tp_flash_returns_structured_result_and_json_like_dict() -> None:
    mix = _hydrocarbon_mixture()

    result = mix.equilibrium(kind="tp_flash", T=220.0, P=1.0e5, z=np.asarray([0.1, 0.3, 0.6]))

    assert isinstance(result, epcsaft.EquilibriumResult)
    assert result.backend == "neutral_vle"
    assert result.problem_kind == "tp_flash"
    assert result.phase_labels == ["liq", "vap"]
    payload = result.to_dict()
    assert payload["phase_labels"] == ["liq", "vap"]
    json.dumps(payload, allow_nan=False)
    _assert_json_like(payload)


def test_explicit_flash_tp_matches_legacy_equilibrium_dispatch() -> None:
    mix = _hydrocarbon_mixture()
    feed = np.asarray([0.1, 0.3, 0.6])

    direct = mix.flash_tp(T=220.0, P=1.0e5, z=feed)
    legacy = mix.equilibrium(kind="tp_flash", T=220.0, P=1.0e5, z=feed)

    assert isinstance(direct, epcsaft.EquilibriumResult)
    assert direct.problem_kind == legacy.problem_kind
    assert direct.phase_labels == legacy.phase_labels
    np.testing.assert_allclose(direct.phases[0].composition, legacy.phases[0].composition)
    np.testing.assert_allclose(direct.phases[1].composition, legacy.phases[1].composition)


def test_explicit_stability_tp_matches_legacy_equilibrium_dispatch() -> None:
    feed = np.asarray([0.1, 0.3, 0.6])

    direct = _hydrocarbon_mixture().stability_tp(T=300.0, P=1.0e5, z=feed, parent_phase="liq", trial_phases=("liq",))
    legacy = _hydrocarbon_mixture().equilibrium(
        kind="stability",
        T=300.0,
        P=1.0e5,
        z=feed,
        parent_phase="liq",
        trial_phases=("liq",),
    )

    assert isinstance(direct, epcsaft.StabilityResult)
    assert direct.problem_kind == legacy.problem_kind
    assert direct.parent_phase == legacy.parent_phase
    assert direct.trial_phase == legacy.trial_phase
    assert direct.min_tpd == pytest.approx(legacy.min_tpd)


def test_explicit_lle_tp_matches_legacy_equilibrium_dispatch() -> None:
    mix = ePCSAFTMixture.from_params(
        {
            "m": np.asarray([1.5255, 2.5303]),
            "s": np.asarray([3.2300, 3.8499]),
            "e": np.asarray([188.90, 278.11]),
            "e_assoc": np.asarray([2899.5, 0.0]),
            "vol_a": np.asarray([0.035176, 0.0]),
            "assoc_scheme": ["2B", None],
            "k_ij": np.asarray([[0.0, 0.051], [0.051, 0.0]]),
        },
        species=["Methanol", "Cyclohexane"],
    )
    feed = np.asarray([0.5, 0.5])

    direct = mix.lle_tp(T=298.15, P=1.013e5, z=feed)
    legacy = mix.equilibrium(kind="lle_flash", T=298.15, P=1.013e5, z=feed)

    assert isinstance(direct, epcsaft.EquilibriumResult)
    assert direct.problem_kind == legacy.problem_kind
    assert direct.phase_labels == legacy.phase_labels
    assert direct.split_detected == legacy.split_detected


def test_explicit_chemical_equilibrium_matches_legacy_equilibrium_dispatch() -> None:
    mix = ePCSAFTMixture.from_params(
        {
            "m": np.asarray([1.0, 1.0]),
            "s": np.asarray([3.0, 3.0]),
            "e": np.asarray([200.0, 200.0]),
        },
        species=["A", "B"],
    )
    request = {
        "T": 298.15,
        "P": 1.0e5,
        "z": [0.5, 0.5],
        "balances": {"total": {"A": 1.0, "B": 1.0}},
        "totals": {"total": 1.0},
        "reactions": [epcsaft.ReactionDefinition({"A": -1.0, "B": 1.0}, np.log(3.0))],
        "options": epcsaft.ReactiveSpeciationOptions(tolerance=1.0e-10),
    }

    direct = mix.chemical_equilibrium(**request)
    legacy = mix.equilibrium(kind="chemical_equilibrium", **request)

    assert isinstance(direct, epcsaft.ReactiveSpeciationResult)
    assert direct.x == pytest.approx(legacy.x)
    assert direct.reaction_residuals == pytest.approx(legacy.reaction_residuals)


def test_equilibrium_phase_exposes_explicit_ln_fugacity_alias() -> None:
    mix = _hydrocarbon_mixture()

    result = mix.equilibrium(kind="tp_flash", T=220.0, P=1.0e5, z=np.asarray([0.1, 0.3, 0.6]))

    phase = result.phases[0]
    np.testing.assert_allclose(phase.ln_fugacity_coefficient, phase.fugacity_coefficient)
    payload = phase.to_dict()
    assert "ln_fugacity_coefficient" in payload
    assert "fugacity_coefficient" in payload
    np.testing.assert_allclose(payload["ln_fugacity_coefficient"], payload["fugacity_coefficient"])
    json.dumps(result.to_dict(), allow_nan=False)


@pytest.mark.parametrize(
    ("kwargs", "match"),
    [
        ({"kind": "bubble_point", "T": 220.0, "P": 1.0e5, "z": [0.1, 0.3, 0.6]}, "Only kind='tp_flash'"),
        ({"kind": "tp_flash", "P": 1.0e5, "z": [0.1, 0.3, 0.6]}, "T"),
        ({"kind": "tp_flash", "T": 220.0, "z": [0.1, 0.3, 0.6]}, "P"),
        ({"kind": "tp_flash", "T": 220.0, "P": 1.0e5}, "z"),
        ({"kind": "tp_flash", "T": 220.0, "P": 1.0e5, "z": [1.0]}, "length"),
        ({"kind": "tp_flash", "T": 220.0, "P": 1.0e5, "z": [0.1, 0.3, -0.4]}, "non-negative"),
    ],
)
def test_equilibrium_rejects_invalid_public_inputs(kwargs, match) -> None:
    mix = _hydrocarbon_mixture()

    with pytest.raises(epcsaft.InputError, match=match):
        mix.equilibrium(**kwargs)


def test_equilibrium_rejects_ionic_mixtures_for_v1() -> None:
    mix = _ionic_mixture()

    with pytest.raises(epcsaft.InputError, match="ion-containing"):
        mix.equilibrium(kind="tp_flash", T=298.15, P=1.0e5, z=[0.9998, 1.0e-4, 1.0e-4])
