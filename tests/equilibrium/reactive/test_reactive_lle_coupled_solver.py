from __future__ import annotations

import numpy as np
import pytest

import epcsaft
from epcsaft import ePCSAFTMixture
from tests.helpers.numeric import assert_allclose


def _neutral_reactive_lle_fixture() -> tuple[ePCSAFTMixture, np.ndarray, dict[str, object], epcsaft.ReactionDefinition]:
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
    mix = ePCSAFTMixture.from_params(params, species=["Methanol", "Cyclohexane"])
    liq1 = np.asarray([0.11757838279937723, 0.8824216172006228])
    liq2 = np.asarray([0.7985874308392054, 0.20141256916079467])
    beta2 = 0.48813098468607985
    feed = (1.0 - beta2) * liq1 + beta2 * liq2
    reaction = epcsaft.ReactionDefinition.from_literature_constant(
        {"Methanol": -1.0, "Cyclohexane": 1.0},
        log_equilibrium_constant=-0.079259405371,
        name="methanol_to_cyclohexane",
        standard_state="mole_fraction_activity",
        source="repo-contained model-consistent reactive LLE fixture",
    )
    return mix, feed, {"liq1": liq1, "liq2": liq2, "phase_fraction": beta2}, reaction


def _assert_coupled_reactive_phase_diagnostics(diagnostics: dict[str, object]) -> None:
    assert diagnostics["reactive_workflow_class"] == "coupled_native"
    assert diagnostics["reactive_phase_method"] == "native_coupled_reactive_phase_equilibrium"
    assert diagnostics["staged_route_used"] is False
    assert diagnostics["solver_backend"] == "ceres"
    assert diagnostics["selected_solver_backend"] == "ceres"
    assert diagnostics["solver_method"] == "ceres_trust_region_coupled_reactive_phase_equilibrium"
    assert diagnostics["jacobian_backend"] == "cppad_implicit"
    assert diagnostics["derivative_backend"] == "cppad_implicit"
    assert diagnostics["jacobian_available"] is True
    assert diagnostics["derivative_available"] is True
    assert diagnostics["solved_state_sensitivity_backend"] == "cppad_implicit"
    assert diagnostics["reaction_and_phase_residuals_share_state"] is True


def test_neutral_reactive_lle_public_route_uses_coupled_native_solver() -> None:
    mix, feed, initial_phases, reaction = _neutral_reactive_lle_fixture()

    result = mix.equilibrium(
        kind="reactive_lle",
        T=298.15,
        P=1.013e5,
        z=feed,
        balances={"total": {"Methanol": 1.0, "Cyclohexane": 1.0}},
        totals={"total": 1.0},
        reactions=[reaction],
        initial_phases=initial_phases,
        phase_options=epcsaft.EquilibriumOptions(max_iterations=80, tolerance=1.0e-8, min_composition=1.0e-12),
    )
    diagnostics = result.diagnostics

    assert result.problem_kind == "reactive_phase_equilibrium"
    assert result.split_detected is True
    _assert_coupled_reactive_phase_diagnostics(diagnostics)
    assert diagnostics["phase_kind"] == "lle_flash"
    assert diagnostics["reaction_residual_norm"] <= 1.0e-8
    assert diagnostics["phase_equilibrium_residual_norm"] <= 1.0e-8
    assert diagnostics["material_balance_norm"] <= 1.0e-10
    assert diagnostics["element_balance_norm"] <= 1.0e-10
    assert diagnostics["phase_charge_balance_norm"] <= 1.0e-10
    assert diagnostics["phase_distance"] > 0.6
    assert set(diagnostics["phase_compositions"]) == {"liq1", "liq2"}
    assert diagnostics["phase_fraction_sum"] == pytest.approx(1.0, abs=1.0e-12)
    assert set(diagnostics["reaction_extents"]) == {"methanol_to_cyclohexane"}
    assert "chemical_equilibrium_then_phase_equilibrium" not in diagnostics.values()

    reconstructed = np.zeros_like(feed)
    for phase in result.phases:
        assert phase.composition.sum() == pytest.approx(1.0, abs=1.0e-12)
        reconstructed += phase.phase_fraction * phase.composition
    assert_allclose(reconstructed, feed, atol=1.0e-10)
