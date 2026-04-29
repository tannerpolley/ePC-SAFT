from __future__ import annotations

import numpy as np
import pytest

import epcsaft
from epcsaft import ePCSAFTMixture
from epcsaft.parameters import get_prop_dict


MW_WATER_BUTANOL = np.asarray([18.01528e-3, 74.1216e-3], dtype=float)


def _mass_to_mole_fraction(mass_fraction: list[float]) -> np.ndarray:
    moles = np.asarray(mass_fraction, dtype=float) / MW_WATER_BUTANOL
    return moles / np.sum(moles)


def _water_butanol_mixture() -> ePCSAFTMixture:
    temperature = 298.15
    species = ["H2O", "Butanol"]
    params = get_prop_dict("2022_Ascani", species, [0.5, 0.5], temperature)
    return ePCSAFTMixture.from_params(params, species=species)


def _water_butanol_case_study() -> tuple[np.ndarray, np.ndarray, np.ndarray, dict[str, object]]:
    aqueous = _mass_to_mole_fraction([1.0 - 0.07001321004, 0.07001321004])
    organic = _mass_to_mole_fraction([1.0 - 0.7926023778, 0.7926023778])
    feed = 0.5 * aqueous + 0.5 * organic
    initial_phases = {
        "liq1": aqueous,
        "liq2": organic,
        "phase_fraction": 0.5,
    }
    return feed, aqueous, organic, initial_phases


def _assert_json_like(value):
    if isinstance(value, dict):
        for item in value.values():
            _assert_json_like(item)
    elif isinstance(value, list):
        for item in value:
            _assert_json_like(item)
    else:
        assert not isinstance(value, np.ndarray)


@pytest.mark.xfail(
    reason=(
        "The current 2022_Ascani neutral water/1-butanol surface collapses "
        "to one liquid phase; keep this as the V2 split-acceptance target."
    ),
    strict=True,
)
def test_water_butanol_lle_flash_closes_material_and_fugacity_balance() -> None:
    mix = _water_butanol_mixture()
    feed, _aqueous_seed, _organic_seed, initial_phases = _water_butanol_case_study()

    result = mix.equilibrium(
        kind="lle_flash",
        T=298.15,
        P=1.0e5,
        z=feed,
        backend="neutral_lle",
        initial_phases=initial_phases,
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
    assert np.max(np.abs(liq1.composition - liq2.composition)) > 1.0e-4

    reconstructed = liq1.phase_fraction * liq1.composition + liq2.phase_fraction * liq2.composition
    np.testing.assert_allclose(reconstructed, feed, atol=1.0e-10)
    assert result.diagnostics["material_balance_error"] < 1.0e-10
    assert result.diagnostics["fugacity_residual_norm"] < 1.0e-6

    fugacity_residual = (
        np.log(liq2.composition)
        + liq2.fugacity_coefficient
        - np.log(liq1.composition)
        - liq1.fugacity_coefficient
    )
    np.testing.assert_allclose(fugacity_residual, np.zeros_like(feed), atol=1.0e-6)

    payload = result.to_dict()
    assert payload["phase_labels"] == ["liq1", "liq2"]
    _assert_json_like(payload)


def test_water_butanol_lle_case_study_reports_current_native_no_split_blocker() -> None:
    mix = _water_butanol_mixture()
    feed, _aqueous_seed, _organic_seed, initial_phases = _water_butanol_case_study()

    result = mix.equilibrium(
        kind="lle_flash",
        T=298.15,
        P=1.0e5,
        z=feed,
        initial_phases=initial_phases,
    )

    assert result.split_detected is False
    assert result.stable is True
    assert result.backend == "neutral_lle"
    assert result.problem_kind == "lle_flash"
    assert result.phase_labels == ["liq"]
    assert "no V2 LLE split" in result.diagnostics["message"]
    assert result.diagnostics["phase_distance"] < 5.0e-2
    assert result.diagnostics["fugacity_residual_norm"] < 2.0e-1


def test_lle_flash_reports_no_split_for_identical_initial_phases() -> None:
    mix = _water_butanol_mixture()
    feed, _aqueous_seed, _organic_seed, _initial_phases = _water_butanol_case_study()

    result = mix.equilibrium(
        kind="lle_flash",
        T=298.15,
        P=1.0e5,
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
    mix = _water_butanol_mixture()

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
