from __future__ import annotations

import json

import numpy as np
import pytest

import epcsaft
from epcsaft import ePCSAFTMixture
from tests.equilibrium.core.test_stability import _assert_stability_route_pending
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
    with pytest.raises(epcsaft.InputError) as stability_exc:
        mix.solve_equilibrium(
            epcsaft.StabilityAnalysis(T=300.0, P=1.0e5, z=feed, parent_phase="liq", trial_phases=("liq",))
        )

    _assert_tp_flash_route_pending(flash_exc)
    _assert_stability_route_pending(stability_exc)

    with pytest.raises(epcsaft.InputError, match=r"dew_p requires a native Ipopt equilibrium NLP route"):
        mix.solve_equilibrium(epcsaft.DewPoint(T=260.0, y=feed))


def test_equilibrium_phase_exposes_ln_and_coefficient_fugacity_fields() -> None:
    ln_phi = np.asarray([0.0, 0.1, 0.2])
    phase = epcsaft.EquilibriumPhase(
        "liq",
        composition=np.asarray([0.1, 0.3, 0.6]),
        density=10.0,
        temperature=220.0,
        pressure=1.0e5,
        phase_fraction=1.0,
        ln_fugacity_coefficient=ln_phi,
    )
    assert_allclose(phase.ln_fugacity_coefficient, ln_phi)
    assert_allclose(phase.fugacity_coefficient, np.exp(ln_phi))
    payload = phase.to_dict()
    assert "ln_fugacity_coefficient" in payload
    assert "fugacity_coefficient" in payload
    assert_allclose(payload["ln_fugacity_coefficient"], ln_phi)
    assert_allclose(payload["fugacity_coefficient"], np.exp(ln_phi))
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
