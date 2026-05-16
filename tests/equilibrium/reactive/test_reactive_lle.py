from __future__ import annotations

import math

import numpy as np
import pytest

import epcsaft
from epcsaft import ePCSAFTMixture
from tests.helpers.numeric import assert_allclose


def _reactive_lle_mixture() -> ePCSAFTMixture:
    params = {
        "MW": np.asarray([32.042e-3, 84.147e-3]),
        "m": np.asarray([1.5255, 2.5303]),
        "s": np.asarray([3.2300, 3.8499]),
        "e": np.asarray([188.90, 278.11]),
        "e_assoc": np.asarray([2899.5, 0.0]),
        "vol_a": np.asarray([0.035176, 0.0]),
        "assoc_scheme": ["2B", None],
        "k_ij": np.asarray([[0.0, 0.051], [0.051, 0.0]]),
        "z": np.asarray([0.0, 0.0]),
        "dielc": np.asarray([33.05, 2.02]),
    }
    return ePCSAFTMixture.from_params(params, species=["Methanol", "Cyclohexane"])


def _lle_seed() -> tuple[np.ndarray, dict[str, object]]:
    phase_a = np.asarray([0.05, 0.95], dtype=float)
    phase_b = np.asarray([0.85, 0.15], dtype=float)
    feed = 0.5 * phase_a + 0.5 * phase_b
    return feed, {"liq1": phase_a, "liq2": phase_b, "phase_fraction": 0.5}


def _reaction_for_feed(feed: np.ndarray) -> epcsaft.ReactionDefinition:
    return epcsaft.ReactionDefinition.from_literature_constant(
        {"Methanol": -1.0, "Cyclohexane": 1.0},
        log_equilibrium_constant=math.log(float(feed[1] / feed[0])),
        name="methanol_to_cyclohexane",
        standard_state="ideal_mole_fraction",
        source="generic reactive LLE regression fixture",
    )


def _phase_result(feed: np.ndarray) -> epcsaft.EquilibriumResult:
    phase_a = epcsaft.EquilibriumPhase(
        "liq1",
        [0.05, 0.95],
        density=850.0,
        temperature=298.15,
        pressure=1.013e5,
        phase_fraction=0.5,
        ln_fugacity_coefficient=[0.1, 0.2],
    )
    phase_b = epcsaft.EquilibriumPhase(
        "liq2",
        [0.85, 0.15],
        density=780.0,
        temperature=298.15,
        pressure=1.013e5,
        phase_fraction=0.5,
        ln_fugacity_coefficient=[0.1, 0.2],
    )
    return epcsaft.EquilibriumResult(
        backend="unit_test_phase_route",
        problem_kind="lle_flash",
        phases=(phase_a, phase_b),
        stable=False,
        split_detected=True,
        diagnostics={
            "equilibrium_route": "lle_flash",
            "fugacity_residual_norm": 1.0e-12,
            "material_balance_error": 1.0e-12,
            "phase_distance": 0.8,
            "feed_composition": feed.tolist(),
        },
    )


def test_explicit_reactive_staged_equilibrium_routes_reaction_coordinates_into_neutral_lle_split(monkeypatch) -> None:
    mix = _reactive_lle_mixture()
    feed, initial_phases = _lle_seed()
    monkeypatch.setattr(
        mix,
        "lle_tp",
        lambda *, T, P, z, options=None, initial_phases=None: _phase_result(np.asarray(z, dtype=float)),
    )

    result = mix.reactive_staged_equilibrium(
        T=298.15,
        P=1.013e5,
        z=[0.5, 0.5],
        balances={"total": {"Methanol": 1.0, "Cyclohexane": 1.0}},
        totals={"total": 1.0},
        reactions=[_reaction_for_feed(feed)],
        phase_kind="lle_flash",
        phase_options=epcsaft.EquilibriumOptions(max_iterations=240, tolerance=1.0e-10, damping=0.5),
        phase_kwargs={"initial_phases": initial_phases},
    )

    assert result.success is True
    assert result.diagnostics["reactive_phase_method"] == "chemical_equilibrium_then_phase_equilibrium"
    assert result.diagnostics["coupling_level"] == "staged_not_full_simultaneous_nlp"
    assert result.diagnostics["full_simultaneous_reactive_nlp"] is False
    assert result.diagnostics["derivative_policy"]["numerical_derivative_backend_available"] is False
    assert result.diagnostics["phase_kind"] == "lle_flash"
    assert result.diagnostics["reaction_coordinates"]["reaction_count"] == 1
    assert result.diagnostics["reaction_coordinates"]["named_reactions"] == ["methanol_to_cyclohexane"]
    assert result.diagnostics["reaction_equilibrium_residuals"]["methanol_to_cyclohexane"] == pytest.approx(
        0.0, abs=1.0e-8
    )
    assert result.diagnostics["element_balance_residuals"]["total"] == pytest.approx(0.0, abs=1.0e-10)
    assert result.diagnostics["nonnegativity"]["status"] == "pass"
    assert result.diagnostics["phase_split"]["status"] == "split_detected"
    assert result.diagnostics["phase_split"]["phase_labels"] == ["liq1", "liq2"]
    assert result.diagnostics["fugacity_equality"]["fugacity_residual_norm"] < 1.0e-8
    assert result.diagnostics["material_balance_error"] < 1.0e-8
    assert_allclose([result.z["Methanol"], result.z["Cyclohexane"]], feed, atol=1.0e-10)


def test_explicit_reactive_staged_equilibrium_routes_generic_lle() -> None:
    mix = _reactive_lle_mixture()
    feed, initial_phases = _lle_seed()
    mix.lle_tp = lambda *, T, P, z, options=None, initial_phases=None: _phase_result(np.asarray(z, dtype=float))
    result = mix.reactive_staged_equilibrium(
        T=298.15,
        P=1.013e5,
        z=[0.5, 0.5],
        balances={"total": {"Methanol": 1.0, "Cyclohexane": 1.0}},
        totals={"total": 1.0},
        reactions=[_reaction_for_feed(feed)],
        phase_kind="lle_flash",
        phase_options=epcsaft.EquilibriumOptions(max_iterations=240, tolerance=1.0e-10, damping=0.5),
        phase_kwargs={"initial_phases": initial_phases},
    )

    assert result.success is True
    assert result.phase.split_detected is True
    assert result.diagnostics["ascani_benchmark_attempt"]["status"] == "not_applicable_to_neutral_route"
