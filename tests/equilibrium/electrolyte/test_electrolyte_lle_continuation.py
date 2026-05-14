from __future__ import annotations

import os

import numpy as np
import pytest

pytestmark = pytest.mark.skipif(
    os.environ.get("EPCSAFT_RUN_HUBACH_LLE") != "1",
    reason="Hubach continuation solve is an opt-in hard-case regression.",
)

import epcsaft
from epcsaft import initial_phases_from_result
from tests.equilibrium.electrolyte.test_hubach_electrolyte_lle import (
    P_PA,
    T_K,
    _hubach_mixture,
    _row0_feed,
    _row0_initial_phases,
)


def test_initial_phases_from_result_round_trips_hubach_split() -> None:
    feed = _row0_feed()
    mix = _hubach_mixture(feed)
    result = mix.equilibrium(
        kind="electrolyte_lle",
        T=T_K,
        P=P_PA,
        z=feed,
        initial_phases=_row0_initial_phases(),
        options=epcsaft.EquilibriumOptions(max_iterations=180, tolerance=1.0e-8),
    )

    seed = initial_phases_from_result(result)

    assert set(seed) == {"aq", "org", "phase_fraction"}
    assert seed["phase_fraction"] == result.phases[1].phase_fraction


def test_equilibrium_curve_uses_previous_hubach_split_as_seed() -> None:
    feed0 = _row0_feed()
    feed1 = feed0.copy()
    feed1[0] -= 0.002
    feed1[1] += 0.0015
    feed1[2] += 0.0005
    feed1 = feed1 / float(np.sum(feed1))
    mix = _hubach_mixture(feed0)

    results = mix.equilibrium_curve(
        [{"z": feed0}, {"z": feed1}],
        kind="electrolyte_lle",
        T=T_K,
        P=P_PA,
        initial_phases=_row0_initial_phases(),
        options=epcsaft.EquilibriumOptions(max_iterations=180, tolerance=1.0e-8),
    )

    assert len(results) == 2
    assert all(result.split_detected for result in results)
    assert results[1].diagnostics["solver_seed_name"] == "initial_phases"
