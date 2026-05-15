from __future__ import annotations

import numpy as np
import pytest

import epcsaft
from epcsaft import ePCSAFTMixture
from tests.equilibrium.reactive.test_reactive_lle_coupled_solver import _assert_coupled_reactive_phase_diagnostics


def _reactive_electrolyte_lle_fixture() -> tuple[
    ePCSAFTMixture,
    np.ndarray,
    dict[str, object],
    epcsaft.ReactionDefinition,
]:
    species = ["H2O", "Butanol", "Na+", "Cl-"]
    aq = np.asarray([0.798324680201737, 0.016320352824141723, 0.09267748348706063, 0.09267748348706063])
    org = np.asarray([0.37006036048879404, 0.6214918588210971, 0.004223890345054407, 0.004223890345054407])
    beta_org = 0.613766575013417
    feed = (1.0 - beta_org) * aq + beta_org * org
    mix = ePCSAFTMixture.from_dataset("2022_Ascani", species, feed, 298.15)
    reaction = epcsaft.ReactionDefinition.from_literature_constant(
        {"H2O": -1.0, "Butanol": 1.0},
        log_equilibrium_constant=-1.265500953237746,
        name="water_to_butanol",
        standard_state="mole_fraction_activity",
        source="repo-contained model-consistent reactive electrolyte LLE fixture",
    )
    return mix, feed, {"aq": aq, "org": org, "phase_fraction": beta_org}, reaction


def test_reactive_electrolyte_lle_public_route_uses_coupled_native_solver() -> None:
    mix, feed, initial_phases, reaction = _reactive_electrolyte_lle_fixture()

    result = mix.equilibrium(
        kind="reactive_electrolyte_lle",
        T=298.15,
        P=1.013e5,
        z=feed,
        balances={
            "solvent_total": {"H2O": 1.0, "Butanol": 1.0},
            "sodium": {"Na+": 1.0},
            "chloride": {"Cl-": 1.0},
        },
        totals={
            "solvent_total": float(feed[0] + feed[1]),
            "sodium": float(feed[2]),
            "chloride": float(feed[3]),
        },
        reactions=[reaction],
        initial_phases=initial_phases,
        phase_options=epcsaft.EquilibriumOptions(max_iterations=80, tolerance=1.0e-8, min_composition=1.0e-12),
    )
    diagnostics = result.diagnostics

    assert result.problem_kind == "reactive_phase_equilibrium"
    assert result.split_detected is True
    _assert_coupled_reactive_phase_diagnostics(diagnostics)
    assert diagnostics["phase_kind"] == "electrolyte_lle"
    assert diagnostics["reaction_residual_norm"] <= 1.0e-8
    assert diagnostics["phase_equilibrium_residual_norm"] <= 2.0e-8
    assert diagnostics["ionic_equilibrium_residual_norm"] <= 1.0e-8
    assert diagnostics["material_balance_norm"] <= 1.0e-8
    assert diagnostics["element_balance_norm"] <= 1.0e-8
    assert diagnostics["phase_charge_balance_norm"] <= 1.0e-8
    assert diagnostics["phase_distance"] > 0.5
    assert set(diagnostics["element_balance_residuals"]) == {"solvent_total", "sodium", "chloride"}
    assert set(diagnostics["reaction_extents"]) == {"water_to_butanol"}

    charges = np.asarray(mix.parameters["z"], dtype=float)
    reconstructed = np.zeros_like(feed)
    for phase in result.phases:
        assert phase.composition.sum() == pytest.approx(1.0, abs=1.0e-12)
        assert float(np.dot(phase.composition, charges)) == pytest.approx(0.0, abs=1.0e-8)
        reconstructed += phase.phase_fraction * phase.composition
    np.testing.assert_allclose(
        [
            reconstructed[0] + reconstructed[1],
            reconstructed[2],
            reconstructed[3],
        ],
        [feed[0] + feed[1], feed[2], feed[3]],
        atol=1.0e-8,
    )
