from __future__ import annotations

import numpy as np
import pytest

import epcsaft
from epcsaft import ePCSAFTMixture


def _methanol_cyclohexane_mixture() -> ePCSAFTMixture:
    # Gross/Sadowski 2002 methanol parameters and kij; Gross/Sadowski 2001 cyclohexane parameters.
    params = {
        "MW": np.asarray([32.042e-3, 84.147e-3]),
        "m": np.asarray([1.5255, 2.5303]),
        "s": np.asarray([3.2300, 3.8499]),
        "e": np.asarray([188.90, 278.11]),
        "e_assoc": np.asarray([2899.5, 0.0]),
        "vol_a": np.asarray([0.035176, 0.0]),
        "assoc_scheme": ["2B", None],
        "k_ij": np.asarray(
            [
                [0.0, 0.051],
                [0.051, 0.0],
            ]
        ),
        "z": np.asarray([0.0, 0.0]),
        "dielc": np.asarray([33.05, 2.02]),
    }
    return ePCSAFTMixture.from_params(params, species=["Methanol", "Cyclohexane"])


def _methanol_cyclohexane_lle_benchmark() -> tuple[np.ndarray, dict[str, object]]:
    methanol_poor = np.asarray([0.05, 0.95], dtype=float)
    methanol_rich = np.asarray([0.85, 0.15], dtype=float)
    feed = 0.5 * methanol_poor + 0.5 * methanol_rich
    initial_phases = {
        "liq1": methanol_poor,
        "liq2": methanol_rich,
        "phase_fraction": 0.5,
    }
    return feed, initial_phases


def _assert_json_like(value):
    if isinstance(value, dict):
        for item in value.values():
            _assert_json_like(item)
    elif isinstance(value, list):
        for item in value:
            _assert_json_like(item)
    else:
        assert not isinstance(value, np.ndarray)


def test_methanol_cyclohexane_lle_flash_closes_material_and_fugacity_balance() -> None:
    mix = _methanol_cyclohexane_mixture()
    feed, initial_phases = _methanol_cyclohexane_lle_benchmark()

    result = mix.equilibrium(
        kind="lle_flash",
        T=298.15,
        P=1.013e5,
        z=feed,
        backend="neutral_lle",
        initial_phases=initial_phases,
        options=epcsaft.EquilibriumOptions(max_iterations=240, tolerance=1.0e-10, damping=0.5),
    )

    assert result.split_detected is True
    assert result.stable is False
    assert result.backend == "neutral_lle"
    assert result.problem_kind == "lle_flash"
    assert result.phase_labels == ["liq1", "liq2"]
    liq1, liq2 = result.phases
    assert 0.0 < liq1.phase_fraction < 1.0
    assert 0.0 < liq2.phase_fraction < 1.0
    np.testing.assert_allclose(liq1.composition.sum(), 1.0)
    np.testing.assert_allclose(liq2.composition.sum(), 1.0)
    assert np.all(liq1.composition > 0.0)
    assert np.all(liq2.composition > 0.0)

    np.testing.assert_allclose(liq1.composition, [0.1175783826, 0.8824216174], atol=5.0e-8)
    np.testing.assert_allclose(liq2.composition, [0.7985874309, 0.2014125691], atol=5.0e-8)
    np.testing.assert_allclose(liq2.phase_fraction, 0.4881309848, atol=5.0e-8)
    assert result.diagnostics["phase_distance"] > 0.65

    reconstructed = liq1.phase_fraction * liq1.composition + liq2.phase_fraction * liq2.composition
    np.testing.assert_allclose(reconstructed, feed, atol=1.0e-10)
    assert result.diagnostics["material_balance_error"] < 1.0e-10
    assert result.diagnostics["fugacity_residual_norm"] < 1.0e-9

    fugacity_residual = (
        np.log(liq2.composition)
        + liq2.fugacity_coefficient
        - np.log(liq1.composition)
        - liq1.fugacity_coefficient
    )
    np.testing.assert_allclose(fugacity_residual, np.zeros_like(feed), atol=1.0e-9)

    payload = result.to_dict()
    assert payload["phase_labels"] == ["liq1", "liq2"]
    _assert_json_like(payload)


def test_lle_flash_reports_no_split_for_identical_initial_phases() -> None:
    mix = _methanol_cyclohexane_mixture()
    feed, _initial_phases = _methanol_cyclohexane_lle_benchmark()

    result = mix.equilibrium(
        kind="lle_flash",
        T=298.15,
        P=1.013e5,
        z=feed,
        initial_phases={"liq1": feed, "liq2": feed, "phase_fraction": 0.5},
    )

    assert result.split_detected is False
    assert result.stable is True
    assert result.phase_labels == ["liq"]
    assert "no V2 LLE split" in result.diagnostics["message"]


@pytest.mark.parametrize(
    ("kwargs", "match"),
    [
        ({"kind": "lle_flash", "T": 298.15, "P": 1.0e5, "z": [0.5, 0.5], "backend": "native"}, "backend"),
        ({"kind": "lle_flash", "P": 1.0e5, "z": [0.5, 0.5]}, "T"),
        ({"kind": "lle_flash", "T": 298.15, "z": [0.5, 0.5]}, "P"),
        ({"kind": "lle_flash", "T": 298.15, "P": 1.0e5}, "z"),
        ({"kind": "lle_flash", "T": 298.15, "P": 1.0e5, "z": [1.0]}, "length"),
        ({"kind": "lle_flash", "T": 298.15, "P": 1.0e5, "z": [0.5, -0.5]}, "non-negative"),
        (
            {
                "kind": "lle_flash",
                "T": 298.15,
                "P": 1.0e5,
                "z": [0.5, 0.5],
                "initial_phases": {"liq1": [0.5, 0.5], "phase_fraction": 0.5},
            },
            "initial_phases",
        ),
        (
            {
                "kind": "lle_flash",
                "T": 298.15,
                "P": 1.0e5,
                "z": [0.5, 0.5],
                "initial_phases": {"liq1": [0.5, 0.5], "liq2": [0.2, 0.8], "phase_fraction": 1.2},
            },
            "phase_fraction",
        ),
    ],
)
def test_lle_flash_rejects_invalid_public_inputs(kwargs, match) -> None:
    mix = _methanol_cyclohexane_mixture()

    with pytest.raises(epcsaft.InputError, match=match):
        mix.equilibrium(**kwargs)


def test_lle_flash_rejects_ionic_mixtures_for_v2() -> None:
    params = {
        "m": np.asarray([1.2047, 1.0, 1.0]),
        "s": np.asarray([2.7927, 2.8232, 2.7560]),
        "e": np.asarray([353.95, 230.0, 170.0]),
        "z": np.asarray([0.0, 1.0, -1.0]),
        "dielc": np.asarray([78.09, 8.0, 8.0]),
    }
    mix = ePCSAFTMixture.from_params(params, species=["water", "Na+", "Cl-"])

    with pytest.raises(epcsaft.InputError, match="ion-containing"):
        mix.equilibrium(kind="lle_flash", T=298.15, P=1.0e5, z=[0.9998, 1.0e-4, 1.0e-4])
