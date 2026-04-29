from __future__ import annotations

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


def _methanol_cyclohexane_mixture() -> ePCSAFTMixture:
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


def _assert_json_like(value):
    if isinstance(value, dict):
        for item in value.values():
            _assert_json_like(item)
    elif isinstance(value, list):
        for item in value:
            _assert_json_like(item)
    else:
        assert not isinstance(value, np.ndarray)


def test_stability_public_exports_are_available() -> None:
    assert hasattr(epcsaft, "StabilityTrial")
    assert hasattr(epcsaft, "StabilityResult")


def test_stability_returns_structured_result_and_json_like_dict() -> None:
    mix = _hydrocarbon_mixture()

    result = mix.equilibrium(
        kind="stability",
        T=300.0,
        P=1.0e5,
        z=[0.1, 0.3, 0.6],
        parent_phase="vap",
        trial_phases=("vap",),
    )

    assert isinstance(result, epcsaft.StabilityResult)
    assert result.backend == "neutral_tpd"
    assert result.problem_kind == "stability"
    assert result.parent_phase == "vap"
    assert result.trial_phase == "vap"
    np.testing.assert_allclose(result.trial_composition.sum(), 1.0)
    assert result.diagnostics["stability_analysis"] == "neutral_tpd"
    assert result.diagnostics["parent_phases"] == ["vap"]
    assert result.diagnostics["trial_phases"] == ["vap"]
    _assert_json_like(result.to_dict())


def test_explicit_stability_runs_even_when_equilibrium_precheck_is_disabled() -> None:
    mix = _hydrocarbon_mixture()

    result = mix.equilibrium(
        kind="stability",
        T=300.0,
        P=1.0e5,
        z=[0.1, 0.3, 0.6],
        parent_phase="vap",
        trial_phases=("vap",),
        options=epcsaft.EquilibriumOptions(stability_precheck=False),
    )

    assert result.diagnostics["stability_analysis"] == "neutral_tpd"
    assert result.diagnostics["trial_count"] > 0


def test_methanol_cyclohexane_liquid_parent_detects_liquid_instability() -> None:
    mix = _methanol_cyclohexane_mixture()

    result = mix.equilibrium(
        kind="stability",
        T=298.15,
        P=1.013e5,
        z=[0.45, 0.55],
        parent_phase="liq",
        trial_phases=("liq",),
    )

    assert result.stable is False
    assert result.min_tpd < -1.0e-4
    assert result.parent_phase == "liq"
    assert result.trial_phase == "liq"
    assert result.trial_composition[0] > 0.75
    assert result.trial_composition[1] < 0.25
    assert any(trial.unstable for trial in result.trials)


def test_hydrocarbon_vapor_parent_vapor_trial_is_stable() -> None:
    mix = _hydrocarbon_mixture()

    result = mix.equilibrium(
        kind="stability",
        T=300.0,
        P=1.0e5,
        z=[0.1, 0.3, 0.6],
        parent_phase="vap",
        trial_phases=("vap",),
    )

    assert result.stable is True
    assert result.min_tpd >= -1.0e-6
    assert all(not trial.unstable for trial in result.trials)


def test_stability_without_parent_phase_runs_liquid_and_vapor_parent_roots() -> None:
    mix = _methanol_cyclohexane_mixture()

    result = mix.equilibrium(
        kind="stability",
        T=298.15,
        P=1.013e5,
        z=[0.45, 0.55],
        trial_phases=("liq", "vap"),
    )

    parent_phases = {trial.parent_phase for trial in result.trials}
    assert parent_phases == {"liq", "vap"}
    assert result.diagnostics["parent_phases"] == ["liq", "vap"]
    assert result.diagnostics["trial_phases"] == ["liq", "vap"]


@pytest.mark.parametrize(
    ("kwargs", "match"),
    [
        ({"kind": "stability", "T": 300.0, "P": 1.0e5, "z": [0.1, 0.3, 0.6], "backend": "native"}, "backend"),
        ({"kind": "stability", "P": 1.0e5, "z": [0.1, 0.3, 0.6]}, "T"),
        ({"kind": "stability", "T": 300.0, "z": [0.1, 0.3, 0.6]}, "P"),
        ({"kind": "stability", "T": 300.0, "P": 1.0e5}, "z"),
        ({"kind": "stability", "T": 300.0, "P": 1.0e5, "z": [1.0]}, "length"),
        ({"kind": "stability", "T": 300.0, "P": 1.0e5, "z": [0.1, -0.3, 0.6]}, "non-negative"),
        (
            {"kind": "stability", "T": 300.0, "P": 1.0e5, "z": [0.1, 0.3, 0.6], "parent_phase": "solid"},
            "parent_phase",
        ),
        (
            {"kind": "stability", "T": 300.0, "P": 1.0e5, "z": [0.1, 0.3, 0.6], "trial_phases": ("liq", "solid")},
            "trial_phases",
        ),
    ],
)
def test_stability_rejects_invalid_public_inputs(kwargs, match) -> None:
    mix = _hydrocarbon_mixture()

    with pytest.raises(epcsaft.InputError, match=match):
        mix.equilibrium(**kwargs)


def test_stability_rejects_ionic_mixtures_for_v3() -> None:
    params = {
        "m": np.asarray([1.2047, 1.0, 1.0]),
        "s": np.asarray([2.7927, 2.8232, 2.7560]),
        "e": np.asarray([353.95, 230.0, 170.0]),
        "z": np.asarray([0.0, 1.0, -1.0]),
        "dielc": np.asarray([78.09, 8.0, 8.0]),
    }
    mix = ePCSAFTMixture.from_params(params, species=["water", "Na+", "Cl-"])

    with pytest.raises(epcsaft.InputError, match="ion-containing"):
        mix.equilibrium(kind="stability", T=298.15, P=1.0e5, z=[0.9998, 1.0e-4, 1.0e-4])
