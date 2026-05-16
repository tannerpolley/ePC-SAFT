from __future__ import annotations

import math

import numpy as np
import pytest

import epcsaft
from epcsaft import ePCSAFTMixture
from tests.equilibrium.electrolyte.test_electrolyte_lle_solver_contracts import _assert_ceres_production_diagnostics
from tests.helpers.numeric import assert_allclose


def _salting_out_fixture() -> tuple[ePCSAFTMixture, np.ndarray, dict[str, object]]:
    species = ["H2O", "Butanol", "Na+", "Cl-"]
    aq = np.asarray([0.798324680201737, 0.016320352824141723, 0.09267748348706063, 0.09267748348706063])
    org = np.asarray([0.37006036048879404, 0.6214918588210971, 0.004223890345054407, 0.004223890345054407])
    beta_org = 0.613766575013417
    feed = (1.0 - beta_org) * aq + beta_org * org
    mix = ePCSAFTMixture.from_dataset("2022_Ascani", species, feed, 298.15)
    initial_phases = {"aq": aq, "org": org, "phase_fraction": beta_org}
    return mix, feed, initial_phases


def test_quaternary_salting_out_lle_benchmark_uses_production_solver() -> None:
    mix, feed, initial_phases = _salting_out_fixture()

    result = mix.equilibrium(
        kind="electrolyte_lle",
        T=298.15,
        P=1.013e5,
        z=feed,
        initial_phases=initial_phases,
        options=epcsaft.EquilibriumOptions(max_iterations=80, tolerance=1.0e-8, min_composition=1.0e-12),
    )
    diagnostics = result.diagnostics

    assert result.split_detected is True
    assert len(result.phases) == 2
    _assert_ceres_production_diagnostics(diagnostics)
    assert diagnostics["neutral_fugacity_residual_norm"] <= 1.0e-8
    assert diagnostics["ionic_equilibrium_residual_norm"] <= 1.0e-8
    assert diagnostics["material_balance_norm"] <= 1.0e-10
    assert diagnostics["phase_charge_balance_norm"] <= 1.0e-8
    assert diagnostics["gibbs_delta"] < 0.0
    assert diagnostics["phase_distance"] > 0.1

    phases = {phase.label: phase for phase in result.phases}
    charges = np.asarray(mix.parameters["z"], dtype=float)
    reconstructed = np.zeros_like(feed)
    for phase in result.phases:
        assert math.isfinite(float(phase.density))
        assert phase.composition.sum() == pytest.approx(1.0, abs=1.0e-12)
        assert float(np.dot(phase.composition, charges)) == pytest.approx(0.0, abs=1.0e-8)
        reconstructed += phase.phase_fraction * phase.composition
    assert_allclose(reconstructed, feed, atol=1.0e-10)
    assert phases["aq"].composition[0] > phases["org"].composition[0]
    assert phases["org"].composition[1] > phases["aq"].composition[1]
    assert phases["aq"].composition[2] > phases["org"].composition[2]
