from __future__ import annotations

import json

import numpy as np
import pytest

import epcsaft
from epcsaft import ePCSAFTMixture
from tests.equilibrium.core.test_vle import _assert_tp_flash_route_pending
from tests.helpers.numeric import assert_allclose


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


def test_tp_flash_requires_native_ipopt_route() -> None:
    mix = _hydrocarbon_mixture()

    with pytest.raises(epcsaft.InputError) as excinfo:
        mix.equilibrium(kind="tp_flash", T=220.0, P=1.0e5, z=np.asarray([0.1, 0.3, 0.6]))

    _assert_tp_flash_route_pending(excinfo)


def test_explicit_flash_tp_matches_equilibrium_route_pending_policy() -> None:
    mix = _hydrocarbon_mixture()
    feed = np.asarray([0.1, 0.3, 0.6])

    with pytest.raises(epcsaft.InputError) as direct_exc:
        mix.flash_tp(T=220.0, P=1.0e5, z=feed)
    with pytest.raises(epcsaft.InputError) as dispatched_exc:
        mix.equilibrium(kind="tp_flash", T=220.0, P=1.0e5, z=feed)

    _assert_tp_flash_route_pending(direct_exc)
    _assert_tp_flash_route_pending(dispatched_exc)


@pytest.mark.parametrize(
    ("kind", "kwargs", "route"),
    [
        ("bubble_p", {"T": 220.0, "x_liq": [0.2, 0.3, 0.5]}, "bubble_p"),
        ("bubble_t", {"P": 1.0e5, "z": [0.2, 0.3, 0.5]}, "bubble_t"),
        ("dew_p", {"T": 220.0, "z": [0.1, 0.3, 0.6]}, "dew_p"),
        ("dew_t", {"P": 1.0e5, "z": [0.1, 0.3, 0.6]}, "dew_t"),
    ],
)
def test_equilibrium_dispatch_rejects_neutral_bubble_dew_until_native_ipopt_route(
    kind: str, kwargs: dict[str, object], route: str
) -> None:
    mix = _hydrocarbon_mixture()

    with pytest.raises(epcsaft.InputError, match=rf"{route} requires a native Ipopt equilibrium NLP route"):
        mix.equilibrium(kind=kind, **kwargs)


def test_solve_equilibrium_accepts_typed_problem_objects() -> None:
    mix = _hydrocarbon_mixture()
    feed = np.asarray([0.1, 0.3, 0.6])

    with pytest.raises(epcsaft.InputError) as flash_exc:
        mix.solve_equilibrium(epcsaft.TPFlash(T=220.0, P=1.0e5, z=feed))
    stability = mix.solve_equilibrium(
        epcsaft.StabilityAnalysis(T=300.0, P=1.0e5, z=feed, parent_phase="liq", trial_phases=("liq",))
    )

    _assert_tp_flash_route_pending(flash_exc)
    assert isinstance(stability, epcsaft.StabilityResult)
    assert stability.problem_kind == "stability"

    with pytest.raises(epcsaft.InputError, match=r"dew_p requires a native Ipopt equilibrium NLP route"):
        mix.solve_equilibrium(epcsaft.DewPoint(T=260.0, y=feed))


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


def test_explicit_lle_tp_matches_equilibrium_route_pending_policy() -> None:
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

    with pytest.raises(epcsaft.InputError, match=r"lle_flash requires a native Ipopt equilibrium NLP route"):
        mix.lle_tp(T=298.15, P=1.013e5, z=feed)
    with pytest.raises(epcsaft.InputError, match=r"lle_flash requires a native Ipopt equilibrium NLP route"):
        mix.equilibrium(kind="lle_flash", T=298.15, P=1.013e5, z=feed)


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
        "reactions": [
            epcsaft.ReactionDefinition(
                {"A": -1.0, "B": 1.0},
                np.log(3.0),
                standard_state="ideal_mole_fraction",
            )
        ],
        "options": epcsaft.ReactiveSpeciationOptions(tolerance=1.0e-10),
    }

    direct = mix.chemical_equilibrium(**request)
    legacy = mix.equilibrium(kind="chemical_equilibrium", **request)

    assert isinstance(direct, epcsaft.ReactiveSpeciationResult)
    assert direct.x == pytest.approx(legacy.x)
    assert direct.reaction_residuals == pytest.approx(legacy.reaction_residuals)


def test_equilibrium_phase_exposes_explicit_ln_fugacity_alias() -> None:
    phase = epcsaft.EquilibriumPhase(
        "liq",
        composition=np.asarray([0.1, 0.3, 0.6]),
        density=10.0,
        temperature=220.0,
        pressure=1.0e5,
        phase_fraction=1.0,
        ln_fugacity_coefficient=np.asarray([0.0, 0.1, 0.2]),
    )
    assert_allclose(phase.ln_fugacity_coefficient, phase.fugacity_coefficient)
    payload = phase.to_dict()
    assert "ln_fugacity_coefficient" in payload
    assert "fugacity_coefficient" in payload
    assert_allclose(payload["ln_fugacity_coefficient"], payload["fugacity_coefficient"])
    json.dumps(payload, allow_nan=False)


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
