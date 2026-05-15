from __future__ import annotations

import math

import numpy as np
import pytest

import epcsaft
from tests.equilibrium.electrolyte.test_electrolyte_lle_solver_contracts import (
    _assert_ceres_production_diagnostics,
    _case2_feed,
    _case2_mixture,
)


def test_distributed_ion_lle_production_solver_reports_residual_proof() -> None:
    feed = _case2_feed()
    mix = _case2_mixture(feed)

    result = mix.equilibrium(
        kind="electrolyte_lle",
        T=298.15,
        P=1.0e5,
        z=feed,
        options=epcsaft.EquilibriumOptions(max_iterations=180, tolerance=1.0e-8, damping=0.5),
    )
    diagnostics = result.diagnostics

    assert result.split_detected is True
    assert len(result.phases) == 2
    _assert_ceres_production_diagnostics(diagnostics)
    assert diagnostics["acceptance_gate"] == "ceres_residual_solve"
    assert diagnostics["neutral_fugacity_residual_norm"] <= 1.0e-8
    assert diagnostics["ionic_equilibrium_residual_norm"] <= 1.0e-8
    assert diagnostics["material_balance_norm"] <= 1.0e-10
    assert diagnostics["phase_charge_balance_norm"] <= 1.0e-8
    assert diagnostics["scaled_solver_residual_norm"] <= 1.0e-8
    assert diagnostics["unscaled_solver_residual_norm"] <= 1.0e-8
    assert diagnostics["gibbs_delta"] < 0.0
    assert diagnostics["phase_distance"] > 0.1

    charges = np.asarray(mix.parameters["z"], dtype=float)
    reconstructed = np.zeros_like(feed)
    for phase in result.phases:
        assert math.isfinite(float(phase.density))
        assert phase.density > 0.0
        assert phase.composition.sum() == pytest.approx(1.0, abs=1.0e-12)
        assert float(np.dot(phase.composition, charges)) == pytest.approx(0.0, abs=1.0e-8)
        reconstructed += phase.phase_fraction * phase.composition
    np.testing.assert_allclose(reconstructed, feed, atol=1.0e-10)

    salt_pairs = diagnostics["salt_pairs"]
    assert [pair["label"] for pair in salt_pairs] == ["NaCl", "KCl"]
