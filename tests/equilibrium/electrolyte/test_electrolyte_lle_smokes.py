from __future__ import annotations

import numpy as np
import pytest

import epcsaft
from epcsaft import _core, ePCSAFTMixture


def _ascani_water_butanol_nacl_mixture(feed=None) -> ePCSAFTMixture:
    species = ["H2O", "Butanol", "Na+", "Cl-"]
    if feed is None:
        feed = np.asarray([0.55, 0.40, 0.025, 0.025], dtype=float)
    return ePCSAFTMixture.from_dataset("2022_Ascani", species, feed, 298.15)


def _assert_electrolyte_lle_native_ipopt_gate(excinfo: pytest.ExceptionInfo[epcsaft.InputError]) -> None:
    message = str(excinfo.value)
    assert "electrolyte_lle requires a native Ipopt equilibrium NLP route" in message


def test_electrolyte_lle_builds_native_route_before_ipopt_gate(monkeypatch) -> None:
    feed = np.asarray([0.55, 0.40, 0.025, 0.025], dtype=float)
    mix = _ascani_water_butanol_nacl_mixture(feed)
    calls: list[tuple[object, ...]] = []

    def fake_route(_native, *route_args):
        calls.append(route_args)
        return {
            "backend": "ipopt",
            "compiled": False,
            "ran": False,
            "accepted": False,
            "status": "ipopt_dependency_required",
            "postsolve": {"accepted": False},
        }

    monkeypatch.setattr(_core, "_native_electrolyte_lle_eos_route_result", fake_route)

    with pytest.raises(epcsaft.InputError) as excinfo:
        mix.equilibrium(
            kind="electrolyte_lle",
            T=298.15,
            P=1.013e5,
            z=feed,
            options=epcsaft.EquilibriumOptions(
                solver_backend="ipopt",
                max_iterations=70,
                tolerance=1.0e-7,
                timeout_seconds=9.0,
                hessian_mode="exact",
                ipopt_iteration_history_limit=4,
            ),
        )

    _assert_electrolyte_lle_native_ipopt_gate(excinfo)
    assert len(calls) == 1
    (
        temperature,
        target_pressure,
        feed_amounts,
        max_iterations,
        tolerance,
        timeout_seconds,
        hessian_mode,
        iteration_history_limit,
        material_tolerance,
        pressure_tolerance,
        charge_tolerance,
        chemical_potential_tolerance,
        phase_distance_tolerance,
    ) = calls[0]
    assert temperature == pytest.approx(298.15)
    assert target_pressure == pytest.approx(1.013e5)
    assert feed_amounts == pytest.approx(feed.tolist())
    assert max_iterations == 70
    assert tolerance == pytest.approx(1.0e-7)
    assert timeout_seconds == pytest.approx(9.0)
    assert hessian_mode == "exact"
    assert iteration_history_limit == 4
    assert material_tolerance == pytest.approx(1.0e-7)
    assert pressure_tolerance == pytest.approx(1.013e-2)
    assert charge_tolerance == pytest.approx(1.0e-8)
    assert chemical_potential_tolerance == pytest.approx(1.0e-7)
    assert phase_distance_tolerance > 0.0


def test_electrolyte_lle_molality_feed_requires_native_ipopt_route() -> None:
    mix = _ascani_water_butanol_nacl_mixture()

    result = mix.equilibrium(
        kind="electrolyte_lle",
        T=298.15,
        P=1.013e5,
        solvent_feed={"H2O": 0.58, "Butanol": 0.42},
        salt_molality={"NaCl": 1.0},
    )

    assert result.backend == "native_equilibrium_nlp"
    assert result.problem_kind == "electrolyte_lle"
    assert result.diagnostics["route_status"] == "accepted"
    assert result.diagnostics["solver_backend"] == "ipopt"
    assert result.diagnostics["charge_balance_norm"] <= 1.0e-8


def test_electrolyte_lle_rejects_non_neutral_direct_feed() -> None:
    mix = _ascani_water_butanol_nacl_mixture()

    with pytest.raises(epcsaft.InputError, match="charge neutral"):
        mix.equilibrium(kind="electrolyte_lle", T=298.15, P=1.013e5, z=[0.55, 0.40, 0.04, 0.01])
