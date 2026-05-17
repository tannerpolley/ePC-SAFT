from __future__ import annotations

import numpy as np
import pytest

import epcsaft
from epcsaft import _core, ePCSAFTMixture


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


def _assert_stability_native_ipopt_gate(excinfo: pytest.ExceptionInfo[epcsaft.InputError], route: str = "stability") -> None:
    message = str(excinfo.value)
    assert f"{route} requires a native Ipopt equilibrium stability NLP route" in message


def test_stability_uses_native_ipopt_route_after_validation() -> None:
    mix = _hydrocarbon_mixture()

    try:
        result = mix.equilibrium(
            kind="stability",
            T=300.0,
            P=1.0e5,
            z=[0.1, 0.3, 0.6],
            parent_phase="vap",
            trial_phases=("vap",),
        )
    except epcsaft.InputError as exc:
        assert "stability requires a native Ipopt equilibrium stability NLP route" in str(exc)
        return

    assert isinstance(result, epcsaft.StabilityResult)
    assert result.backend == "native_equilibrium_nlp"
    assert result.problem_kind == "neutral_stability"
    assert result.diagnostics["derivative_backend"] == "cppad_implicit"


def test_stability_builds_one_native_route_request_before_ipopt_gate(monkeypatch: pytest.MonkeyPatch) -> None:
    mix = _hydrocarbon_mixture()
    calls: list[dict[str, object]] = []

    def fake_route(
        _native,
        temperature,
        pressure,
        feed_composition,
        parent_phase,
        trial_phase,
        max_iterations,
        tolerance,
        timeout_seconds,
        stability_tolerance,
    ):
        calls.append(
            {
                "temperature": temperature,
                "pressure": pressure,
                "feed_composition": feed_composition,
                "parent_phase": parent_phase,
                "trial_phase": trial_phase,
                "max_iterations": max_iterations,
                "tolerance": tolerance,
                "timeout_seconds": timeout_seconds,
                "stability_tolerance": stability_tolerance,
            }
        )
        return {
            "backend": "ipopt",
            "compiled": False,
            "ran": False,
            "accepted": False,
            "status": "ipopt_dependency_required",
        }

    monkeypatch.setattr(_core, "_native_neutral_stability_tpd_route_result", fake_route)

    with pytest.raises(epcsaft.InputError, match=r"stability requires a native Ipopt equilibrium stability NLP route"):
        mix.equilibrium(
            kind="stability",
            T=300.0,
            P=1.0e5,
            z=[0.1, 0.3, 0.6],
            parent_phase="vap",
            trial_phases=("vap",),
            options=epcsaft.EquilibriumOptions(max_iterations=19, tolerance=3.0e-8, timeout_seconds=4.5),
        )

    assert len(calls) == 1
    call = calls[0]
    assert call["temperature"] == pytest.approx(300.0)
    assert call["pressure"] == pytest.approx(1.0e5)
    assert call["feed_composition"] == pytest.approx([0.1, 0.3, 0.6])
    assert call["parent_phase"] == "vap"
    assert call["trial_phase"] == "vap"
    assert call["max_iterations"] == 19
    assert call["tolerance"] == pytest.approx(3.0e-8)
    assert call["timeout_seconds"] == pytest.approx(4.5)
    assert call["stability_tolerance"] == pytest.approx(3.0e-8)


def test_stability_converts_accepted_native_route_payload(monkeypatch: pytest.MonkeyPatch) -> None:
    mix = _hydrocarbon_mixture()

    def fake_route(
        _native,
        temperature,
        pressure,
        feed_composition,
        parent_phase,
        trial_phase,
        max_iterations,
        tolerance,
        timeout_seconds,
        stability_tolerance,
    ):
        assert temperature == pytest.approx(300.0)
        assert pressure == pytest.approx(1.0e5)
        assert feed_composition == pytest.approx([0.1, 0.3, 0.6])
        assert parent_phase == "vap"
        assert trial_phase == "vap"
        assert max_iterations > 0
        assert tolerance > 0.0
        assert timeout_seconds == 0.0
        assert stability_tolerance > 0.0
        return {
            "backend": "ipopt",
            "compiled": True,
            "ran": True,
            "solver_accepted": True,
            "accepted": True,
            "stable": True,
            "status": "accepted",
            "solver_status": "Solve_Succeeded",
            "application_status": "Solve_Succeeded",
            "parent_phase": "vap",
            "trial_phase": "vap",
            "seed_name": "canonical_shifted_feed",
            "min_tpd": 2.0e-6,
            "objective": 2.0e-6,
            "trial_composition": [0.12, 0.28, 0.60],
            "constraints": [0.0],
            "derivative_backend": "cppad_implicit",
        }

    monkeypatch.setattr(_core, "_native_neutral_stability_tpd_route_result", fake_route)

    result = mix.equilibrium(
        kind="stability",
        T=300.0,
        P=1.0e5,
        z=[0.1, 0.3, 0.6],
        parent_phase="vap",
        trial_phases=("vap",),
    )

    assert isinstance(result, epcsaft.StabilityResult)
    assert result.backend == "native_equilibrium_nlp"
    assert result.problem_kind == "neutral_stability"
    assert result.stable is True
    assert result.min_tpd == pytest.approx(2.0e-6)
    assert result.parent_phase == "vap"
    assert result.trial_phase == "vap"
    assert result.trial_composition == pytest.approx([0.12, 0.28, 0.60])
    assert len(result.trials) == 1
    assert result.trials[0].converged is True
    assert result.trials[0].unstable is False
    assert result.diagnostics["route_count"] == 1
    assert result.diagnostics["derivative_backend"] == "cppad_implicit"


@pytest.mark.parametrize(
    ("kwargs", "match"),
    [
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
