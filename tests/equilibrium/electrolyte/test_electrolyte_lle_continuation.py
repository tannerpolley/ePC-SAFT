from __future__ import annotations

import os

import numpy as np
import pytest

pytestmark = pytest.mark.skipif(
    os.environ.get("EPCSAFT_RUN_HUBACH_LLE") != "1",
    reason="Hubach continuation solve is an opt-in hard-case regression.",
)

import epcsaft
from tests.equilibrium.electrolyte.test_electrolyte_lle_smokes import _assert_electrolyte_lle_route_pending
from tests.equilibrium.electrolyte.test_hubach_electrolyte_lle import (
    P_PA,
    T_K,
    _hubach_mixture,
    _row0_feed,
    _row0_initial_phases,
)


def test_hubach_initial_phase_seed_requires_native_ipopt_route() -> None:
    feed = _row0_feed()
    mix = _hubach_mixture(feed)
    with pytest.raises(epcsaft.InputError) as excinfo:
        mix.equilibrium(
            kind="electrolyte_lle",
            T=T_K,
            P=P_PA,
            z=feed,
            initial_phases=_row0_initial_phases(),
            options=epcsaft.EquilibriumOptions(max_iterations=180, tolerance=1.0e-8),
        )

    _assert_electrolyte_lle_route_pending(excinfo)


def test_equilibrium_curve_requires_native_ipopt_route_for_hubach_lle() -> None:
    feed0 = _row0_feed()
    feed1 = feed0.copy()
    feed1[0] -= 0.002
    feed1[1] += 0.0015
    feed1[2] += 0.0005
    feed1 = feed1 / float(np.sum(feed1))
    mix = _hubach_mixture(feed0)

    with pytest.raises(epcsaft.InputError) as excinfo:
        mix.equilibrium_curve(
            [{"z": feed0}, {"z": feed1}],
            kind="electrolyte_lle",
            T=T_K,
            P=P_PA,
            initial_phases=_row0_initial_phases(),
            options=epcsaft.EquilibriumOptions(max_iterations=180, tolerance=1.0e-8),
        )

    _assert_electrolyte_lle_route_pending(excinfo)
